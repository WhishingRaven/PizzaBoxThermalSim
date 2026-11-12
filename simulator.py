"""Pizza box contest: CAD thin-wall thermal shells + volumetric lumped pizza/air.

This is an independently specified evaluator for the Fall 2026 assignment,
not the TA's undisclosed reference model. STEP and DXF are mandatory inputs.
STEP coordinates are mm, DXF coordinates are read in its declared units.
"""
from __future__ import annotations
import argparse
import json
import math
import time
from pathlib import Path
from dataclasses import dataclass, asdict
import numpy as np
from scipy import ndimage, sparse
from scipy.sparse.linalg import cg, LinearOperator
from scipy.integrate import solve_ivp
import cadquery as cq
import ezdxf
from shapely.geometry import Polygon
from shapely.ops import unary_union
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.gp import gp_Pnt
from OCP.TopAbs import TopAbs_IN, TopAbs_ON
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEFAULTS = {
    'pizza_diameter_m': 0.300, 'pizza_thickness_m': 0.015,
    'pizza_mass_kg': 0.750, 'pizza_cp_J_kgK': 3000.0,
    'box_k_W_mK': 0.08, 'box_density_kg_m3': 200.0,
    'box_cp_J_kgK': 1400.0,
    'air_density_kg_m3': 1.2, 'air_cp_J_kgK': 1005.0,
    'T_pizza_initial_C': 90.0, 'T_environment_C': 22.0,
    'h_internal_W_m2K': 3.0, 'h_external_W_m2K': 10.0,
    'h_contact_W_m2K': 10.0, 'duration_s': 1200.0,
    'water_fraction_initial': 0.4, 'RH_initial': 0.5,
    'water_activity': 1.0, 'km_m_s': 1e-3,
    'Rv_J_kgK': 461.526,
    'min_box_mass_g': 83.0, 'max_box_mass_g': 249.0,
    'min_pizza_temperature_C': 68.0,
    'nominal_board_thickness_m': 0.003,
}

FACE = [
    # neighbor dz,dy,dx; Q1 face node ids; exposed area axis
    ((-1,0,0), (0,1,2,3), 0),
    (( 1,0,0), (4,5,6,7), 0),
    ((0,-1,0), (0,1,4,5), 1),
    ((0, 1,0), (2,3,6,7), 1),
    ((0,0,-1), (0,2,4,6), 2),
    ((0,0, 1), (1,3,5,7), 2),
]
CORNERS=np.array([(0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1)],dtype=int)
QUAD_MASS = np.array([[4,2,2,1],[2,4,1,2],[2,1,4,2],[1,2,2,4]],dtype=float)/36

class GeometryError(ValueError):
    pass


def extract_flat(dxf_path:Path, declared_units='mm', thickness_m=None):
    """Read closed CUT polylines (outline / holes), report board area and fold info."""
    doc=ezdxf.readfile(str(dxf_path))
    unit=doc.header.get('$INSUNITS',0)
    scales={1:0.0254,2:0.3048,4:0.001,5:0.01,6:1.0}
    if unit not in scales:
        scales={'mm':0.001,'cm':0.01,'m':1.0,'in':0.0254}[declared_units]
    scale=scales.get(unit,0.001)
    msp=doc.modelspace()
    cut_geoms=[]; folds=[]; unsupported=[]
    for ent in msp:
        layer=ent.dxf.layer.upper()
        typ=ent.dxftype()
        if layer not in ('CUT','FOLD'):
            continue
        if typ in ('LWPOLYLINE','POLYLINE'):
            if typ=='LWPOLYLINE':
                coords=[(float(x)*scale,float(y)*scale) for x,y,*_ in ent.get_points('xy')]
                isclosed=ent.closed
                if any(abs(p[4])>1e-10 for p in ent.get_points('xyseb')):
                    unsupported.append('bulge arcs in CUT/FOLD LWPOLYLINE')
            else:
                coords=[(float(v.dxf.location.x)*scale,float(v.dxf.location.y)*scale) for v in ent.vertices]
                isclosed=ent.is_closed
            if layer=='CUT':
                if isclosed and len(coords)>=3:
                    poly=Polygon(coords)
                    if not poly.is_valid:
                        raise GeometryError('Invalid/self-crossing CUT polygon.')
                    cut_geoms.append(poly)
                else:
                    unsupported.append('open CUT polyline')
            else:
                folds.append(coords)
        elif typ=='LINE':
            a=ent.dxf.start; b=ent.dxf.end
            seg=[(a.x*scale,a.y*scale),(b.x*scale,b.y*scale)]
            if layer=='FOLD':folds.append(seg)
            else: unsupported.append('CUT line not joined to closed polygon')
        else:
            unsupported.append(f'{typ} on {layer} layer')
    if unsupported:
        raise GeometryError('Unsupported DXF entities (simplify to closed CUT LWPOLYLINE and straight FOLD lines): '+', '.join(sorted(set(unsupported))))
    if not cut_geoms:
        raise GeometryError('DXF requires at least one closed polyline on CUT layer.')
    # Multiple nested contours are interpreted as through-holes using even-odd parity.
    # Disjoint solids produce MultiPolygon and fail the one-sheet criterion.
    cut_geoms.sort(key=lambda p:p.area, reverse=True)
    shape=cut_geoms[0]
    for poly in cut_geoms[1:]:
        shape=shape.symmetric_difference(poly)
    if shape.geom_type!='Polygon':
        raise GeometryError('Flat pattern must represent exactly one connected sheet (one polygon).')
    if len(folds)<1:
        raise GeometryError('No FOLD lines found. Use DXF layer FOLD.')
    mass_g=shape.area*(DEFAULTS['nominal_board_thickness_m'] if thickness_m is None else thickness_m)*DEFAULTS['box_density_kg_m3']*1000.0
    return {'area_m2':float(shape.area), 'mass_g':float(mass_g),
            'fold_count':len(folds), 'hole_count':len(shape.interiors),
            'unit':unit, 'connected_sheet':True,
            'warning':'2D topology only: folding kinematics, interpenetration and seals need separate verification.'}


def inspect_step(path:Path):
    model=cq.importers.importStep(str(path))
    solids=model.solids().vals()
    if len(solids)!=1:
        raise GeometryError(f'STEP must contain exactly one solid; found {len(solids)}')
    shape=solids[0]
    if not shape.isValid():
        raise GeometryError('STEP B-rep solid is invalid.')
    bb=shape.BoundingBox()
    bbox=np.array([[bb.xmin,bb.xmax],[bb.ymin,bb.ymax],[bb.zmin,bb.zmax]],float)*1e-3
    # We normalize CAD coordinates from mm to meters; for OCC point queries, keep mm.
    vol_mm3=shape.Volume()
    mass_g=vol_mm3*1e-9*DEFAULTS['box_density_kg_m3']*1e3
    # Flat exterior bottom requirement (exact B-rep planar face test).
    planar_bottom_area_mm2=0.0
    for face in shape.Faces():
        fbb=face.BoundingBox()
        if abs(fbb.zmin-bb.zmin)<1e-4 and abs(fbb.zmax-bb.zmin)<1e-4:
            planar_bottom_area_mm2+=face.Area()
    if planar_bottom_area_mm2 < math.pi*(150.0**2)*0.95:
        raise GeometryError(f'Exterior base must have a flat horizontal face supporting the pizza; '
                            f'bottom planar area only {planar_bottom_area_mm2:.1f} mm2.')
    if min(np.diff(bbox,axis=1).ravel())<=0:
        raise GeometryError('STEP has zero dimension.')
    return shape,bbox,float(mass_g)


def voxelize(shape,bbox,h_m,max_cells=2_000_000):
    dims=bbox[:,1]-bbox[:,0]
    nxyz=np.ceil(dims/h_m).astype(int)
    nx,ny,nz=(int(nxyz[0]),int(nxyz[1]),int(nxyz[2]))
    count=nx*ny*nz
    if count>max_cells:
        raise GeometryError(f'Mesh too fine: {count:,} candidate voxels (limit {max_cells:,}). Increase --voxel-mm.')
    dx,dy,dz=dims/nxyz
    # OCC classifier in CAD mm coordinates. Note intentional z,y,x array order.
    xs=(bbox[0,0]+(np.arange(nx)+0.5)*dx)*1000
    ys=(bbox[1,0]+(np.arange(ny)+0.5)*dy)*1000
    zs=(bbox[2,0]+(np.arange(nz)+0.5)*dz)*1000
    classifier=BRepClass3d_SolidClassifier(shape.wrapped)
    occ=np.zeros((nz,ny,nx),dtype=bool)
    tol_mm=0.00001
    for iz,z in enumerate(zs):
        for iy,y in enumerate(ys):
            for ix,x in enumerate(xs):
                classifier.Perform(gp_Pnt(float(x),float(y),float(z)),tol_mm)
                occ[iz,iy,ix]=classifier.State() in (TopAbs_IN,TopAbs_ON)
    if not np.any(occ):
        raise GeometryError('No board material resolved in voxel mesh. Decrease --voxel-mm.')
    return occ,(float(dx),float(dy),float(dz)),nxyz


def label_air(occ):
    void=np.pad(~occ,1,constant_values=True)
    structure=ndimage.generate_binary_structure(3,1)
    lab,nlab=ndimage.label(void,structure)
    ext=int(lab[0,0,0])
    values,counts=np.unique(lab,return_counts=True)
    cavity_ids=[int(v) for v,c in zip(values,counts) if v>0 and v!=ext and c>5]
    if len(cavity_ids)!=1:
        raise GeometryError(f'Expected one sealed air cavity; found {len(cavity_ids)} resolved cavities (IDs {cavity_ids}). ' 
                            'Check lid closure, STEP solid and voxel resolution.')
    return lab,ext,cavity_ids[0]


def validate_pizza_fit(occ,labels,cavity_id,bbox,spacing,params):
    """Require flat full-disk support and a 15 mm-clear pizza volume above it."""
    nz,ny,nx=occ.shape;dx,dy,dz=spacing
    x=bbox[0,0]+(np.arange(nx)+.5)*dx
    y=bbox[1,0]+(np.arange(ny)+.5)*dy
    cx=(bbox[0,0]+bbox[0,1])/2;cy=(bbox[1,0]+bbox[1,1])/2
    rad=params['pizza_diameter_m']/2
    # Use strictly interior cell centers to avoid an aliased pizza boundary.
    in_disk=(x[None,:]-cx)**2+(y[:,None]-cy)**2 < max(rad-0.25*max(dx,dy),0)**2
    if not in_disk.any():raise GeometryError('Could not rasterize pizza area.')
    columns=labels[1:-1,1:-1,1:-1][:,in_disk]==cavity_id
    has_cavity=columns.any(axis=0)
    if not has_cavity.all():
        raise GeometryError(f'Pizza disk obstructed in {np.sum(~has_cavity)} cross-section cells.')
    first=columns.argmax(axis=0)
    if int(first.max())!=int(first.min()):
        raise GeometryError('Pizza support base varies in height across the pizza disk.')
    base_idx=int(first[0])
    need=max(1,int(math.ceil(params['pizza_thickness_m']/dz)))
    if base_idx+need>nz or not np.all(columns[base_idx:base_idx+need]):
        raise GeometryError('15 mm-thick pizza collides with lid or internal ribs.')
    if base_idx==0 or not np.all(occ[base_idx-1][in_disk]):
        raise GeometryError('Pizza disk is not directly and fully supported by solid box floor.')
    return {'pizza_floor_inner_z_mm':float((bbox[2,0]+base_idx*dz)*1000),
            'free_pizza_clearance_mm':float(need*dz*1000),
            'pizza_disk_grid_points':int(np.sum(in_disk))}

def build_mesh(occ,nxyz):
    xyz=np.argwhere(occ).astype(np.int64)
    nx,ny,nz=(int(v) for v in nxyz)
    corners=xyz[:,None,:]+CORNERS[None,:,:]
    globals_=(corners[:,:,0]*(ny+1)*(nx+1)+corners[:,:,1]*(nx+1)+corners[:,:,2]).astype(np.int64)
    uniq, inv=np.unique(globals_,return_inverse=True)
    elem=inv.reshape((-1,8)).astype(np.int32)
    z_idx=uniq//((ny+1)*(nx+1)); y_idx=(uniq//(nx+1))%(ny+1); x_idx=uniq%(nx+1)
    return xyz,elem,np.column_stack((x_idx,y_idx,z_idx)),uniq


def local_hex_K(dx,dy,dz,k):
    # Q1 hexahedral stiffness ∫ grad(N_i).grad(N_j) dV, 2x2x2 Gauss.
    out=np.zeros((8,8),dtype=float)
    eta=1/math.sqrt(3)
    for zz in (-eta,eta):
        for yy in (-eta,eta):
            for xx in (-eta,eta):
                grad=np.zeros((8,3))
                for a,(iz,iy,ix) in enumerate(CORNERS):
                    sx=-1 if ix==0 else 1
                    sy=-1 if iy==0 else 1
                    sz=-1 if iz==0 else 1
                    grad[a,0] = sx*(1+sy*yy)*(1+sz*zz)/(4*dx)
                    grad[a,1] = sy*(1+sx*xx)*(1+sz*zz)/(4*dy)
                    grad[a,2] = sz*(1+sx*xx)*(1+sy*yy)/(4*dz)
                out+=k*grad@grad.T*(dx*dy*dz/8)
    return out


def append_blocks(blocks, nodes, matrix):
    if len(nodes)==0:return
    q=nodes.shape[1]
    blocks.append((np.broadcast_to(nodes[:,:,None],(len(nodes),q,q)).ravel().copy(),
                   np.broadcast_to(nodes[:,None,:],(len(nodes),q,q)).ravel().copy(),
                   np.broadcast_to(matrix[None,:,:],(len(nodes),q,q)).ravel().copy()))


def append_pair(blocks, nodes, lump_id, coefficient):
    # off-diagonal -h A / (#face nodes), for each node
    count=nodes.size
    if count==0:return
    r=nodes.ravel(); c=np.full(count,lump_id,dtype=np.int64); v=np.full(count,-coefficient,dtype=float)
    blocks.append((np.r_[r,c],np.r_[c,r],np.r_[v,v]))


def assemble_system(occ, labels, exterior_id, cavity_id, bbox, steps, elem, idx, spacing, params):
    dx,dy,dz=spacing
    nx,ny,nz=(int(v) for v in steps)
    nbox=int(elem.max())+1; ia=nbox; ip=nbox+1; n=nbox+2
    density=params['box_density_kg_m3']; cb=params['box_cp_J_kgK']
    conduction=local_hex_K(dx,dy,dz,params['box_k_W_mK'])
    blocks=[]
    append_blocks(blocks,elem,conduction)
    capacities=np.zeros(n)
    np.add.at(capacities,elem.ravel(),density*cb*dx*dy*dz/8)
    heat_source=np.zeros(n)
    outer_faces=[]; inner_air_faces=[]; pizza_base_faces=[]; bottom_faces=[]
    contacts=[]
    face_areas=[dx*dy,dx*dz,dy*dz]
    loc0=(bbox[0,0],bbox[1,0],bbox[2,0])
    cx=(bbox[0,0]+bbox[0,1])/2
    cy=(bbox[1,0]+bbox[1,1])/2
    radius=params['pizza_diameter_m']/2
    contact_area=0
    air_area=0
    convection_air_faces=[]
    for side,(d,face_ids,axis) in enumerate(FACE):
        zi,yi,xi=idx[:,0],idx[:,1],idx[:,2]
        neigh=labels[zi+1+d[0],yi+1+d[1],xi+1+d[2]]
        isexterior=(neigh==exterior_id)
        isinterior=(neigh==cavity_id)
        if not (np.any(isexterior) or np.any(isinterior)):continue
        area=face_areas[axis]
        qnodes=elem[:,face_ids]
        if np.any(isexterior):
            exnodes=qnodes[isexterior]
            # Exterior bottom is exactly at the ambient/floor prescribed temperature.
            if side==0:
                bottom_faces.append(exnodes)
            else:
                outer_faces.append((exnodes,area))
                append_blocks(blocks,exnodes,params['h_external_W_m2K']*area*QUAD_MASS)
                np.add.at(heat_source,exnodes.ravel(),params['h_external_W_m2K']*area/4*params['T_environment_C'])
        if np.any(isinterior):
            inner_cells=idx[isinterior]
            innodes=qnodes[isinterior]
            # The upward-facing internal base surface is in contact with the disk of pizza.
            if side==1:
                xc=loc0[0]+(inner_cells[:,2]+0.5)*dx
                yc=loc0[1]+(inner_cells[:,1]+0.5)*dy
                iscontact=(xc-cx)**2+(yc-cy)**2 < radius**2
            else:
                iscontact=np.zeros(len(innodes),dtype=bool)
            if np.any(iscontact):
                p_nodes=innodes[iscontact]
                pizza_base_faces.append((p_nodes,area))
                contact_area+=len(p_nodes)*area
                h=params['h_contact_W_m2K']
                append_blocks(blocks,p_nodes,h*area*QUAD_MASS)
                append_pair(blocks,p_nodes,ip,h*area/4)
            if np.any(~iscontact):
                a_nodes=innodes[~iscontact]
                inner_air_faces.append((a_nodes,area))
                h=params['h_internal_W_m2K']
                append_blocks(blocks,a_nodes,h*area*QUAD_MASS)
                append_pair(blocks,a_nodes,ia,h*area/4)
                air_area+=len(a_nodes)*area
    # Add the two lumped capacities and conductances.
    pizza_area=math.pi*radius**2+math.pi*(2*radius)*params['pizza_thickness_m']
    hpa=params['h_internal_W_m2K']*pizza_area
    cavity_count=int(np.count_nonzero(labels==cavity_id))
    cavity_volume=cavity_count*dx*dy*dz
    pizza_volume=math.pi*radius**2*params['pizza_thickness_m']
    free_air_volume=cavity_volume-pizza_volume
    if free_air_volume<=0:
        raise GeometryError(f'No free cavity air volume after pizza insertion: {free_air_volume:.5g} m^3')
    capacities[ia]=params['air_density_kg_m3']*params['air_cp_J_kgK']*free_air_volume
    capacities[ip]=params['pizza_mass_kg']*params['pizza_cp_J_kgK']
    hpc=params['h_contact_W_m2K']*contact_area
    hac=params['h_internal_W_m2K']*air_area
    blocks.append((np.array([ia,ip,ia,ip,ip,ia]),np.array([ia,ip,ip,ia,ip,ia]),np.array([hpa,hpa,-hpa,-hpa,hpc,hac])))
    # NOTE lumped diagonal for air/box-contact are already supplied per interface;
    # pizza/air additional hpa, pizza/contact hpc, air/box hac.
    if not bottom_faces or not pizza_base_faces or not inner_air_faces:
        raise GeometryError('Failed to identify exterior floor, pizza contact or inner-air surfaces.')
    if abs(contact_area-math.pi*radius**2)/(math.pi*radius**2)>0.12:
        raise GeometryError(f'Pizza bottom contact area {contact_area:.5f} m² differs >12% from pizza disk area. '
                            'Pizza may not fit, sit centered, or touch the base; or mesh too coarse.')
    if contact_area <=0.065:
        raise GeometryError('Insufficient pizza contact surface.')
    # For 4-node Q1 hex, block means positive semidefinite, global lump coupling SPD after floor restriction.
    rr=np.concatenate([b[0] for b in blocks]); cc=np.concatenate([b[1] for b in blocks]); vv=np.concatenate([b[2] for b in blocks])
    H=sparse.coo_matrix((vv,(rr,cc)),shape=(n,n)).tocsr()
    H.sum_duplicates(); del rr,cc,vv,blocks
    fixed=np.unique(np.concatenate(bottom_faces).ravel())
    free_mask=np.ones(n,dtype=bool);free_mask[fixed]=False
    free=np.flatnonzero(free_mask)
    # Validate the pizza disk is inside the void at a small height above the floor.
    base_z_cells=np.min(idx[:,0])
    cavity_2d=labels[base_z_cells+2] if base_z_cells+2 < labels.shape[0] else None
    return {'H':H,'capacity':capacities,'source':heat_source,'fixed':fixed,'free':free,'air_idx':ia,'pizza_idx':ip,
            'contact_area':contact_area,'pizza_area':pizza_area,'cavity_volume':cavity_volume,
            'free_air_volume':free_air_volume,'n_dof':n,'n_elements':len(elem),
            'inner_faces':inner_air_faces,'base_faces':pizza_base_faces,'outside_faces':outer_faces,
            'h_pizza_air':hpa,'h_pizza_contact':hpc,'h_air_box':hac}


def solve_thermal(sys,params,dt,output_dir):
    dur=params['duration_s']; nt=math.ceil(dur/dt);dt=dur/nt
    t=np.linspace(0,dur,nt+1)
    H=sys['H']; C=sys['capacity']; free=sys['free']; fixed=sys['fixed']; ia=sys['air_idx']; ip=sys['pizza_idx']
    # Implicit Euler. All coefficients constant, assemble operator once.
    A=(H+sparse.diags(C/dt,format='csr')).tocsr()
    A_uu=A[free,:][:,free].tocsr()
    f=sys['source'][free]-np.asarray(A[free,:][:,fixed].sum(axis=1)).ravel()*params['T_environment_C']
    diagonal=A_uu.diagonal()
    if np.any(diagonal<=0):raise GeometryError('FEM matrix is not positive definite (non-positive diagonal).')
    invdiag=1/diagonal
    prec=LinearOperator(A_uu.shape,matvec=lambda v:invdiag*v)
    initial=np.full(len(C),params['T_environment_C'])
    initial[ip]=params['T_pizza_initial_C']
    u=initial.copy()
    hist_p=np.empty(nt+1);hist_a=np.empty(nt+1);hist_p[0]=u[ip];hist_a[0]=u[ia]
    # Moisture step uses exactly these wall temperature histories after temperature integration.
    innerface=np.concatenate([q for q,a in sys['inner_faces']],axis=0)
    innerarea=np.concatenate([np.full(len(q),a) for q,a in sys['inner_faces']])
    # Contact faces are excluded from exposed air wall. Partition is exact.
    twall=np.empty((nt+1,len(innerface)),dtype=np.float32)
    twall[0]=np.mean(u[innerface],axis=1)
    fluxes=np.empty((nt+1,5))
    facebase=np.concatenate([q for q,a in sys['base_faces']],axis=0)
    area_base=np.concatenate([np.full(len(q),a) for q,a in sys['base_faces']])
    faceout=np.concatenate([q for q,a in sys['outside_faces']],axis=0)
    area_out=np.concatenate([np.full(len(q),a) for q,a in sys['outside_faces']])
    reaction_H_rows=H[fixed,:]
    def flux(v):
        p,a=v[ip],v[ia]
        qpa=sys['h_pizza_air']*(p-a)
        qpbase=params['h_contact_W_m2K']*np.sum(area_base*(p-v[facebase].mean(axis=1)))
        qab=params['h_internal_W_m2K']*np.sum(innerarea*(a-v[innerface].mean(axis=1)))
        qout=params['h_external_W_m2K']*np.sum(area_out*(v[faceout].mean(axis=1)-params['T_environment_C']))
        qfloor=-np.sum(reaction_H_rows.dot(v)-sys['source'][fixed])
        return qpa, qpbase, qab,qout,qfloor
    fluxes[0]=flux(u)
    cg_failures=0;cg_iterations=[]
    for k in range(1,nt+1):
        rhs=f+C[free]/dt*u[free]
        counter=[0]
        def itercount(x):counter[0]+=1
        sol,info=cg(A_uu,rhs,x0=u[free],M=prec,rtol=2e-8,atol=1e-10,maxiter=700,callback=itercount)
        if info!=0:
            cg_failures+=1
            if info<0:raise RuntimeError(f'FEM CG breakdown at step {k}')
            raise RuntimeError(f'FEM CG did not converge at step {k} ({info} iterations); choose coarser mesh or smaller dt')
        cg_iterations.append(counter[0])
        u[free]=sol
        u[fixed]=params['T_environment_C']
        hist_p[k]=u[ip];hist_a[k]=u[ia]
        twall[k]=np.mean(u[innerface],axis=1)
        fluxes[k]=flux(u)
    np.savetxt(output_dir/'temperature.csv',np.column_stack((t,hist_p,hist_a)),delimiter=',',
               header='time_s,pizza_C,air_C',comments='')
    np.savetxt(output_dir/'heat_flux.csv',np.column_stack((t,fluxes)),delimiter=',',
               header='time_s,pizza_to_air_W,pizza_to_base_W,air_to_box_W,box_to_ambient_W,box_to_floor_W',comments='')
    return {'time':t,'Tp':hist_p,'Ta':hist_a,'wallT':twall,'wallareas':innerarea,'flux':fluxes,
            'cg_iterations_avg':float(np.mean(cg_iterations)), 'cg_failures':cg_failures,
            'last_box_temp':u[:ia].copy()}


def p_sat_Pa(T_C):
    """Buck (1981) saturation pressure over liquid water, within water's relevant range."""
    c=np.asarray(T_C)
    return 611.21*np.exp((18.678-c/234.5)*c/(257.14+c))


def solve_moisture(thermal,sys,params,outdir):
    # Internal wall temperature histories are fixed from thermal Step 1.
    # Condensation only if p_v>psat(Twall) and no artificial negative wall water.
    ts=thermal['time']; wp=thermal['Tp'];ta=thermal['Ta']
    wall=thermal['wallT']; area=thermal['wallareas'];Va=sys['free_air_volume']
    km=params['km_m_s'];Rv=params['Rv_J_kgK']
    pizza_area=sys['pizza_area']
    initial_water=params['pizza_mass_kg']*params['water_fraction_initial']
    pv0=params['RH_initial']*float(p_sat_Pa(params['T_environment_C']))
    vapor0=pv0*Va/(Rv*(params['T_environment_C']+273.15))
    def rates(t,mv):
        Tp=float(np.interp(t,ts,wp));Ta=float(np.interp(t,ts,ta))
        lo=min(np.searchsorted(ts,t,side='right')-1,len(ts)-2);lo=max(0,lo)
        r=np.clip((t-ts[lo])/(ts[lo+1]-ts[lo]),0,1)
        Tw=wall[lo]+r*(wall[lo+1]-wall[lo])
        pv=max(mv,0)*Rv*(Ta+273.15)/Va
        Tfp=(Tp+Ta)/2+273.15
        evaporation=max(0.,km*pizza_area*(params['water_activity']*float(p_sat_Pa(Tp))-pv)/(Rv*Tfp))
        Tfw=(Tw+Ta)/2+273.15
        cond=km*area*np.maximum(pv-p_sat_Pa(Tw),0.)/(Rv*Tfw)
        return evaporation,float(cond.sum()),pv
    # Substep integration via solve_ivp BDF, robust for stiff small cavity vapor reservoir.
    def ode(t,y):
        evap,cond,_=rates(t,float(y[0]))
        return [evap-cond,-evap,cond]
    sol=solve_ivp(ode,(0.,float(ts[-1])),[vapor0,initial_water,0.],t_eval=ts,
                  rtol=1e-7,atol=1e-10,method='BDF',max_step=10)
    if not sol.success:
        raise RuntimeError('Moisture integration failed: '+sol.message)
    mv,mwp,wb=sol.y
    # Signed rates were constrained to nonnegative; numerical integration should conserve water.
    residual=mv+mwp+wb-(vapor0+initial_water)
    if np.max(np.abs(residual))>1e-6:
        raise RuntimeError('Moisture balance exceeded tolerance.')
    retained=100*mwp/initial_water
    pvs=np.array([rates(tt,m)[2] for tt,m in zip(ts,mv)])
    np.savetxt(outdir/'moisture.csv',np.column_stack((ts,mwp*1000,mv*1000,wb*1000,retained,pvs)),
               delimiter=',',header='time_s,pizza_water_g,vapor_water_g,wall_condensate_g,retention_pct,vapor_pressure_Pa',comments='')
    return {'initial_pizza_water_g':float(initial_water*1e3), 'initial_air_vapor_g':float(vapor0*1e3),
            'final_pizza_water_g':float(mwp[-1]*1e3), 'final_wall_condensate_g':float(wb[-1]*1e3),
            'final_air_vapor_g':float(mv[-1]*1e3),'retention_pct':float(retained[-1]),
            'max_balance_error_g':float(max(abs(residual))*1000),'retention_history':retained}


def plot_results(thermal,moisture,outdir):
    ts=thermal['time']/60
    fig,ax=plt.subplots(figsize=(7.6,4.5),dpi=160)
    ax.plot(ts,thermal['Tp'],lw=2,label='Pizza')
    ax.plot(ts,thermal['Ta'],lw=1.4,label='Enclosed air')
    ax.axhline(68,color='gray',lw=1.1,ls='--',label='68 C threshold')
    ax.set(xlabel='Time (min)',ylabel='Temperature (C)',title='Thermal FEM: pizza / air')
    ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(outdir/'temperature.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.6,4.5),dpi=160)
    ax.plot(ts,moisture['retention_history'],lw=2)
    ax.set(xlabel='Time (min)',ylabel='Pizza moisture retention (%)',title='Sequential moisture model')
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(outdir/'moisture.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.6,4.5),dpi=160)
    for i,name in enumerate(['Pizza -> air','Pizza -> base','Air -> box','Box -> outside','Box -> floor']):
        ax.plot(ts,thermal['flux'][:,i],label=name,lw=1.4)
    ax.set(xlabel='Time (min)',ylabel='Heat flow (W)',title='Heat-flow pathways')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(outdir/'heat_flux.png');plt.close(fig)


def run(step_path,dxf_path,outdir,voxel_mm=3.0,dt=10.,dxf_unit='mm',max_cells=2_000_000,strict=False,wall_model='shell',thickness_mm=3.):
    if wall_model=='shell':
        from thinwall import run_shell
        return run_shell(step_path,dxf_path,outdir,voxel_mm,dt,thickness_mm,dxf_unit,max_cells)
    if wall_model!='voxel':raise ValueError('Unknown wall model')
    if thickness_mm!=3.:raise ValueError('Archived voxel model has fixed 3 mm material; use shell for explicit thickness')
    start=time.perf_counter()
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    params=DEFAULTS.copy()
    flat=extract_flat(Path(dxf_path),dxf_unit)
    shape,bbox,mass_step=inspect_step(Path(step_path))
    mass_flat=flat['mass_g']
    mismatch=abs(mass_flat-mass_step)/max(mass_flat,1e-10)
    if mismatch>0.12:
        raise GeometryError(f'DXF sheet mass ({mass_flat:.2f} g) and STEP solid mass ({mass_step:.2f} g) differ by {mismatch:.1%}; must be <=12%. '
                            'Include overlap tabs/layers in both representations.')
    if not params['min_box_mass_g']<=mass_flat<=params['max_box_mass_g']:
        raise GeometryError(f'Box mass {mass_flat:.1f} g not in allowed 83-249 g range')
    occ,spacing,nxyz=voxelize(shape,bbox,voxel_mm*1e-3,max_cells=max_cells)
    lab,ext,cavity=label_air(occ)
    pizza_fit=validate_pizza_fit(occ,lab,cavity,bbox,spacing,params)
    xyz,elem,idx,uniq=build_mesh(occ,nxyz)
    # Check 6-connected board voxels: element bridge-only contacts are not acceptable.
    cc,ncc=ndimage.label(occ,ndimage.generate_binary_structure(3,1))
    if ncc>1:raise GeometryError(f'Voxelized board disconnected: {ncc} components. Reduce voxel size.')
    sys=assemble_system(occ,lab,ext,cavity,bbox,nxyz,elem,xyz,spacing,params)
    thermo=solve_thermal(sys,params,dt,outdir)
    np.savez_compressed(outdir/'thermal_field.npz',node_index_xyz=idx,
                        temperature_C=thermo['last_box_temp'],element_nodes=elem,
                        origin_m=bbox[:,0],cell_size_m=np.array(spacing))
    moisture=solve_moisture(thermo,sys,params,outdir)
    plot_results(thermo,moisture,outdir)
    tp=float(thermo['Tp'][-1]);passed=tp>=params['min_pizza_temperature_C']
    result={
        'model':'3D conforming Cartesian Q1 hexahedral FEM on voxelized STEP board + lumped pizza/air; prescribed sequential vapor balance',
        'official_metrics':{'mass_g':mass_flat,'pizza_temperature_20min_C':tp,
                            'pizza_moisture_retention_20min_pct':moisture['retention_pct'],
                            'temperature_threshold_pass':bool(passed)},
        'geometry':{'step':str(step_path),'dxf':str(dxf_path),'flat':flat,'pizza_fit':pizza_fit,
                    'step_solid_mass_g':mass_step,'mass_mismatch_fraction':mismatch,
                    'bbox_mm':(bbox*1000).tolist(),'voxel_mm_actual':(np.array(spacing)*1000).tolist(),
                    'board_voxel_mass_g':float(np.sum(occ)*np.prod(spacing)*params['box_density_kg_m3']*1000),
                    'volume_mesh_error_pct':float((np.sum(occ)*np.prod(spacing)*params['box_density_kg_m3']*1000/mass_step-1)*100),
                    'cavity_volume_m3':sys['cavity_volume'],'free_air_volume_m3':sys['free_air_volume'],
                    'pizza_base_contact_m2':sys['contact_area']},
        'numerics':{'duration_s':params['duration_s'],'dt_s':float(thermo['time'][1]),
                    'requested_voxel_mm':float(voxel_mm),
                    'voxels_xyz':[int(v) for v in nxyz],'box_elements':sys['n_elements'],
                    'degrees_freedom':sys['n_dof'], 'average_CG_iterations':thermo['cg_iterations_avg'],
                    'geometry_aliasing':'The imported STEP is voxelized; check grid convergence before ranking close results.',
                    'time_integration':'backward Euler (1st order)',
                    'water_balance_error_g':moisture['max_balance_error_g']},
        'moisture':{k:v for k,v in moisture.items() if k!='retention_history'},
        'fluxes_20min_W':dict(zip(['pizza_to_air','pizza_to_base','air_to_box','box_to_ambient','box_to_floor'],
                                [float(v) for v in thermo['flux'][-1]])),
        'assumptions':[ 'One single closed STEP solid with one voxel-resolved sealed cavity.',
            'STEP and DXF geometry in standard CAD units; center inferred from exterior box bounding box.',
            'Pizza rests flat at center on inner base: its underside is coupled via h_contact.',
            'Pizza top + cylindrical side exposed to internal air; pizza is spatially uniform.',
            'Air is well mixed and uniform temperature; no spatial CFD or air-flow heat equation.',
            '3D board conduction discretized by uniform Cartesian Q1 brick FEM; stair-stepped surface geometry.',
            'Thermal radiation, edge leakage, vapor latent heat and evaporation cooling are neglected.',
            'Moisture transport uses Buck water saturation vapor pressure and film temperatures (mean air/surface).',
            'Condensation is unidirectional onto inner walls, without re-evaporation (no negative wall liquid).',
            'Flat DXF validates cut-outline area and fold mark presence, not physical foldability.' ],
        'officiality':'INDEPENDENT ESTIMATE, NOT AN OFFICIAL TA SCORE',
        'runtime_s':float(time.perf_counter()-start),
        'status':'temperature_pass' if passed else 'temperature_fail',
    }
    with open(outdir/'result.json','w',encoding='utf-8') as f:
        json.dump(result,f,indent=2,ensure_ascii=False)
    return result


def main():
    parser=argparse.ArgumentParser(description='Fall 2026 pizza box contest independent 3D FEM scorer')
    parser.add_argument('--step',required=True,help='Folded closed STEP solid file (.step/.stp), mm geometry')
    parser.add_argument('--dxf',required=True,help='One-sheet flat pattern DXF (CUT and FOLD layers)')
    parser.add_argument('--out',default='results')
    parser.add_argument('--mesh-mm','--voxel-mm',dest='voxel_mm',type=float,default=3.,help='Surface element size in mm; --voxel-mm is a compatibility alias')
    parser.add_argument('--wall-model',choices=['shell','voxel'],default='shell',help='Default shell; voxel is the archived comparison model')
    parser.add_argument('--wall-thickness-mm',type=float,default=3.)
    parser.add_argument('--dt',type=float,default=10.,help='Thermal time step (s), default 10')
    parser.add_argument('--dxf-unit',default='mm',choices=['mm','cm','m','in'],help='Only for DXF with no units defined')
    parser.add_argument('--max-cells',type=int,default=2_000_000)
    args=parser.parse_args()
    try:
        result=run(args.step,args.dxf,args.out,args.voxel_mm,args.dt,args.dxf_unit,args.max_cells,wall_model=args.wall_model,thickness_mm=args.wall_thickness_mm)
        print(json.dumps({'official_metrics':result['official_metrics'],'numerics':result['numerics'],
                          'runtime_s':result['runtime_s'],'report':str(Path(args.out)/'result.json')},ensure_ascii=False,indent=2))
    except Exception as exc:
        print(f'ERROR: {type(exc).__name__}: {exc}',flush=True)
        raise

if __name__=='__main__':main()
