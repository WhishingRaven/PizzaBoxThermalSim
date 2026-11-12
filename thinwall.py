"""CAD midsurface thermal shells with analytical half-wall resistance.

One P1 surface temperature carries rho*c*t storage and k*t in-plane conduction.
Each side adds t/(2*k) in series with its film/contact resistance. This reduced
through-thickness model is shared by Python and constant-temperature SHELL131.
"""
from pathlib import Path
import csv
import json
import math
import tempfile
import time

import cadquery as cq
import gmsh
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg, LinearOperator
from OCP.BRepOffsetAPI import BRepOffsetAPI_MakeOffsetShape
from OCP.BRepOffset import BRepOffset_Skin
from OCP.GeomAbs import GeomAbs_Intersection


def half_wall_film(h, thickness, conductivity):
    if not np.all(np.isfinite([h,thickness,conductivity])) or thickness <= 0 or conductivity <= 0 or h <= 0:
        raise ValueError('Thickness, conductivity and film coefficient must be positive')
    return 1.0 / (1.0 / h + thickness / (2.0 * conductivity))


def _offset(solid, distance):
    from simulator import GeometryError
    builder = BRepOffsetAPI_MakeOffsetShape()
    try:
        builder.PerformByJoin(solid.wrapped, distance, 1e-6, BRepOffset_Skin,
                              False, False, GeomAbs_Intersection, True)
        result = cq.Shape.cast(builder.Shape())
        if not builder.IsDone() or not result.isValid() or len(result.Solids()) != 1:
            raise ValueError('offset is not one valid solid')
        return result.Solids()[0]
    except Exception as exc:
        raise GeometryError('Cannot form a nonintersecting constant-thickness midsurface; '
                            'check panel junctions, small features and wall thickness') from exc


def build_shell(step, dxf, size_mm, params, dxf_unit='mm', max_elements=500000):
    from simulator import inspect_step, extract_flat, GeometryError
    t = params['nominal_board_thickness_m']
    if not np.isfinite(size_mm) or size_mm <= 0 or not np.isfinite(t) or t <= 0:
        raise ValueError('Surface mesh size and wall thickness must be finite and positive')
    shape, bbox, step_mass = inspect_step(Path(step))
    flat = extract_flat(Path(dxf), dxf_unit, thickness_m=t)
    mismatch = abs(flat['mass_g'] - step_mass) / flat['mass_g']
    if mismatch > .12 or not params['min_box_mass_g'] <= flat['mass_g'] <= params['max_box_mass_g']:
        raise GeometryError('DXF/STEP mass must agree within 12% and sheet mass must be 83-249 g')
    shells = shape.Shells()
    if len(shells) != 2:
        raise GeometryError('Thin-wall sealed-air model requires outside + one closed cavity; '
                            'openings/vents require a ventilation model and are not silently sealed')
    enclosures = [cq.Solid.makeSolid(s) for s in shells]
    cavity = min(enclosures, key=lambda s: abs(s.Volume()))
    if cavity.Volume() < 0:
        cavity = cq.Solid(cavity.wrapped.Reversed())
    exterior = _offset(cavity, t * 1000)
    rebuilt = exterior.cut(cavity)
    error = (rebuilt.cut(shape).Volume() + shape.cut(rebuilt).Volume()) / shape.Volume()
    if error > 1e-4:
        raise GeometryError(f'STEP is not a uniform {t*1000:g} mm wall ({error:.2%} mismatch). '
                            'Overlapping tabs/multiple layers need an explicit layered surface input; '
                            'they cannot be discarded or inferred from mesh layers.')
    midsolid = _offset(cavity, t * 500)
    mid_area = midsolid.Area() * 1e-6
    step_volume = shape.Volume() * 1e-9
    volume_error = (mid_area * t / step_volume - 1) * 100
    if abs(volume_error) > .5:
        raise GeometryError('Midsurface thin-shell approximation differs >0.5% from STEP volume; '
                            'curvature/junction thickness is too large for this thin-wall model')
    cb = cavity.BoundingBox()
    cx, cy = bbox[:2].mean(axis=1) * 1000
    radius = params['pizza_diameter_m'] * 500
    pizza = cq.Workplane('XY').center(cx, cy).circle(radius).extrude(
        params['pizza_thickness_m']*1000).translate((0, 0, cb.zmin)).val()
    if abs(cavity.intersect(pizza).Volume() / pizza.Volume() - 1) > 1e-6:
        raise GeometryError('Centered pizza is obstructed, lacks flat full-disk support, or touches the lid')
    cavity_volume = cavity.Volume()*1e-9
    air_volume = cavity_volume - pizza.Volume()*1e-9
    if air_volume <= 0:
        raise GeometryError('No free air volume after pizza insertion')
    # Gmsh shares vertices across panel edges and respects the exact disk boundary.
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.option.setNumber('General.NumThreads', 1)
        gmsh.model.add('thin_wall')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'midsurface.step'
            cq.exporters.export(midsolid, str(path))
            volumes = gmsh.model.occ.importShapes(str(path))
        gmsh.model.occ.synchronize()
        floor_z = midsolid.BoundingBox().zmin
        disk = gmsh.model.occ.addDisk(cx, cy, floor_z, radius, radius)
        _, mapping = gmsh.model.occ.fragment(volumes, [(2, disk)])
        gmsh.model.occ.synchronize()
        contact_tags = {tag for dim, tag in mapping[1] if dim == 2}
        boundary = {abs(tag) for dim, tag in gmsh.model.getBoundary(gmsh.model.getEntities(3)) if dim == 2}
        if not contact_tags <= boundary:
            raise GeometryError('Pizza contact disk is not wholly contained in the floor')
        gmsh.option.setNumber('Mesh.MeshSizeMin', size_mm)
        gmsh.option.setNumber('Mesh.MeshSizeMax', size_mm)
        gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
        gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
        gmsh.model.mesh.generate(2)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        tags = np.asarray(tags)[order]
        xyz = np.asarray(coords).reshape(-1, 3)[order]*.001
        triangles, areas, scales, floor_flags, contact_flags = [], [], [], [], []
        for tag in sorted(boundary):
            types, _, blocks = gmsh.model.mesh.getElements(2, tag)
            if list(types) != [2]:
                raise GeometryError('Surface mesher did not return linear triangles')
            tri = np.searchsorted(tags, np.asarray(blocks[0]).reshape(-1, 3))
            points = xyz[tri]
            raw_area = np.linalg.norm(np.cross(points[:,1]-points[:,0], points[:,2]-points[:,0]), axis=1)/2
            exact_area = gmsh.model.occ.getMass(2, tag)*1e-6
            scale = exact_area/raw_area.sum()
            bounds = gmsh.model.getBoundingBox(2, tag)
            is_floor = abs(bounds[2]-floor_z) < 1e-4 and abs(bounds[5]-floor_z) < 1e-4
            triangles.extend(tri.tolist())
            areas.extend((raw_area*scale).tolist())
            scales.extend([scale]*len(tri))
            floor_flags.extend([is_floor]*len(tri))
            contact_flags.extend([tag in contact_tags]*len(tri))
    finally:
        gmsh.finalize()
    tri = np.asarray(triangles, dtype=np.int64)
    if len(tri) > max_elements:
        raise GeometryError(f'Surface mesh has {len(tri)} elements; limit is {max_elements}')
    used, tri = np.unique(tri, return_inverse=True)
    xyz = xyz[used]; tri = tri.reshape(-1, 3)
    area = np.asarray(areas)
    floor = np.asarray(floor_flags); contact = np.asarray(contact_flags)
    # Closed-manifold edge multiplicity prevents gaps/duplicated panels/junctions.
    edges = np.sort(np.concatenate((tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]])), axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    if not np.all(counts == 2):
        raise GeometryError('Midsurface has an open or nonmanifold edge; no thermal bridges are invented')
    if abs(area.sum()/mid_area-1) > 1e-8 or abs(area[contact].sum()/(math.pi*(radius*.001)**2)-1)>1e-8:
        raise GeometryError('Surface area or pizza contact partition lost physical CAD area')
    return dict(xyz=xyz, triangles=tri, area=area, area_scale=np.asarray(scales), floor=floor,
                contact=contact, bbox=bbox, flat=flat, step_mass_g=step_mass,
                area_m2=mid_area, thickness_m=t, volume_m3=mid_area*t,
                mass_g=mid_area*t*params['box_density_kg_m3']*1000, step_volume_m3=step_volume,
                shell_volume_error_pct=volume_error, cavity_volume_m3=cavity_volume,
                free_air_volume_m3=air_volume, mesh_mm=size_mm,
                pizza_fit={'pizza_floor_inner_z_mm': cb.zmin, 'full_disk_support':True},
                interior=[(q,float(a)) for q,a,c in zip(tri,area,contact) if not c])


def coefficients(mesh, p):
    t=mesh['thickness_m'];k=p['box_k_W_mK']
    return dict(air=half_wall_film(p['h_internal_W_m2K'],t,k),
                contact=half_wall_film(p['h_contact_W_m2K'],t,k),
                exterior=half_wall_film(p['h_external_W_m2K'],t,k),floor=2*k/t)


def assemble(mesh,p):
    tri=mesh['triangles'];xyz=mesh['xyz'];a=mesh['area'];n=len(xyz);ia=n;ip=n+1
    pts=xyz[tri];normal=np.cross(pts[:,1]-pts[:,0],pts[:,2]-pts[:,0])
    norm2=(normal*normal).sum(axis=1)
    grad=np.stack([np.cross(normal,pts[:,2]-pts[:,1]),np.cross(normal,pts[:,0]-pts[:,2]),
                   np.cross(normal,pts[:,1]-pts[:,0])],axis=1)/norm2[:,None,None]
    local=p['box_k_W_mK']*mesh['thickness_m']*a[:,None,None]*np.einsum('eik,ejk->eij',grad,grad)
    rows=np.repeat(tri,3,axis=1).ravel();cols=np.tile(tri,(1,3)).ravel()
    blocks=[(rows,cols,local.ravel())]
    cap=np.zeros(n+2);np.add.at(cap,tri.ravel(),np.repeat(p['box_density_kg_m3']*p['box_cp_J_kgK']*mesh['thickness_m']*a/3,3))
    source=np.zeros(n+2);weights={};h=coefficients(mesh,p)
    masks={'air':~mesh['contact'],'contact':mesh['contact'],'exterior':~mesh['floor'],'floor':mesh['floor']}
    for label,mask in masks.items():
        q=tri[mask];ar=a[mask];w=np.zeros(n)
        np.add.at(w,q.ravel(),np.repeat(ar/3,3));weights[label]=w
        mass=h[label]*ar[:,None,None]*(np.ones((3,3))+np.eye(3))[None,:,:]/12
        blocks.append((np.repeat(q,3,axis=1).ravel(),np.tile(q,(1,3)).ravel(),mass.ravel()))
        if label in ('air','contact'):
            lump=ia if label=='air' else ip;nodes=np.flatnonzero(w);v=-h[label]*w[nodes]
            blocks.append((np.r_[nodes,np.full(len(nodes),lump),lump],np.r_[np.full(len(nodes),lump),nodes,lump],np.r_[v,v,h[label]*ar.sum()]))
        else: source[:n]+=h[label]*w*p['T_environment_C']
    pizza_area=math.pi*(p['pizza_diameter_m']/2)**2+math.pi*p['pizza_diameter_m']*p['pizza_thickness_m']
    hpa=p['h_internal_W_m2K']*pizza_area
    blocks.append((np.array([ia,ip,ia,ip]),np.array([ia,ip,ip,ia]),np.array([hpa,hpa,-hpa,-hpa])))
    cap[ia]=p['air_density_kg_m3']*p['air_cp_J_kgK']*mesh['free_air_volume_m3']
    cap[ip]=p['pizza_mass_kg']*p['pizza_cp_J_kgK']
    H=sparse.coo_matrix((np.concatenate([b[2] for b in blocks]),(np.concatenate([b[0] for b in blocks]),np.concatenate([b[1] for b in blocks]))),shape=(n+2,n+2)).tocsr()
    return dict(H=H,capacity=cap,source=source,weights=weights,h=h,hpa=hpa,pizza_area=pizza_area,air_idx=ia,pizza_idx=ip)


def flux(u,sys,p):
    ia=sys['air_idx'];ip=sys['pizza_idx'];wall=u[:ia];a=u[ia];pizza=u[ip];w=sys['weights'];h=sys['h']
    return np.array([sys['hpa']*(pizza-a),h['contact']*(w['contact'].sum()*pizza-w['contact']@wall),
                     h['air']*(w['air'].sum()*a-w['air']@wall),
                     h['exterior']*(w['exterior']@wall-w['exterior'].sum()*p['T_environment_C']),
                     h['floor']*(w['floor']@wall-w['floor'].sum()*p['T_environment_C'])])


def sensors(mesh,count=64):
    if count < 6:raise ValueError('At least six wall samples are required')
    from ansys_bridge import cluster_sensors
    return cluster_sensors(mesh,max_samples=count)


def solve(mesh,p,dt,selected):
    sys=assemble(mesh,p);C=sys['capacity'];n=len(mesh['xyz']);N=math.ceil(p['duration_s']/dt);dt=p['duration_s']/N
    times=np.linspace(0,p['duration_s'],N+1);A=sys['H']+sparse.diags(C/dt)
    diagonal=A.diagonal()
    pre=LinearOperator(A.shape,matvec=lambda x:x/diagonal)
    u=np.full(n+2,p['T_environment_C']);u[-1]=p['T_pizza_initial_C'];initial=u.copy()
    sample=np.array([s['nodes'] for s in selected])-1;area=np.array([s['area_m2'] for s in selected])
    hist=np.empty((N+1,2));wall=np.empty((N+1,len(sample)));flows=np.empty((N+1,5))
    factor=mesh['thickness_m']/(2*p['box_k_W_mK'])*sys['h']['air']
    iterations=[]
    for j in range(N+1):
        if j:
            it=[0]
            def count(_):it[0]+=1
            u,info=cg(A,sys['source']+C/dt*u,x0=u,M=pre,rtol=1e-10,atol=1e-11,maxiter=2000,callback=count)
            if info:raise RuntimeError(f'Thin-wall linear solve failed at step {j}: {info}')
            iterations.append(it[0])
        hist[j]=u[-2:][::-1];center=u[sample].mean(axis=1)
        wall[j]=center+factor*(u[-2]-center)
        flows[j]=flux(u,sys,p)
    balance=(C@(u-initial)+dt*flows[1:,3:].sum())
    if abs(balance)>1e-5*p['pizza_mass_kg']*p['pizza_cp_J_kgK']*(p['T_pizza_initial_C']-p['T_environment_C']):
        raise RuntimeError(f'Thermal energy balance failed: {balance} J')
    return dict(time=times,Tp=hist[:,0],Ta=hist[:,1],wallT=wall,wallareas=area,flux=flows,
                last_box_temp=u[:n],cg_iterations_avg=float(np.mean(iterations)),energy_balance_error_J=float(balance))


def geometry_metrics(mesh):
    return {k:float(mesh[k]) for k in ('area_m2','thickness_m','volume_m3','mass_g','step_volume_m3','step_mass_g','shell_volume_error_pct','cavity_volume_m3','free_air_volume_m3')}


def report(mesh,p,thermal,out,step,dxf,runtime,solver):
    from simulator import solve_moisture,plot_results
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    sys={'free_air_volume':mesh['free_air_volume_m3'],'pizza_area':math.pi*(p['pizza_diameter_m']/2)**2+math.pi*p['pizza_diameter_m']*p['pizza_thickness_m']}
    moisture=solve_moisture(thermal,sys,p,out)
    plot_results(thermal,moisture,out)
    times=thermal['time'];dt=float(times[1]-times[0]);names=['pizza_to_air','pizza_to_base','air_to_box','box_to_ambient','box_to_floor']
    np.savetxt(out/'temperature.csv',np.column_stack((times,thermal['Tp'],thermal['Ta'])),delimiter=',',header='time_s,pizza_C,air_C',comments='')
    np.savetxt(out/'heat_flux.csv',np.column_stack((times,thermal['flux'])),delimiter=',',header='time_s,'+','.join(n+'_W' for n in names),comments='')
    energies=dict(zip(names,(dt*thermal['flux'][1:].sum(axis=0)).tolist()))
    energies['cardboard_to_environment']=energies['box_to_ambient']+energies['box_to_floor']
    temp=float(thermal['Tp'][-1]);passed=temp>=p['min_pizza_temperature_C']
    result=dict(model_version='thinwall-v1',model='CAD midsurface P1 thermal shell + analytical half-wall interfaces + volumetric lumped pizza/air',solver=solver,geometry_mesh_mode='shell',
        official_metrics={'mass_g':mesh['flat']['mass_g'],'pizza_temperature_20min_C':temp,'pizza_moisture_retention_20min_pct':moisture['retention_pct'],'temperature_threshold_pass':passed},
        geometry=dict(step=str(step),dxf=str(dxf),flat=mesh['flat'],**geometry_metrics(mesh)),
        material={'conductivity_W_mK':p['box_k_W_mK'],'density_kg_m3':p['box_density_kg_m3'],'specific_heat_J_kgK':p['box_cp_J_kgK']},
        runtime_scope='CAD import + meshing + thermal solve' if solver.startswith('Python') else 'MAPDL solver process wall time; excludes local mesh preparation',
        numerics={'mesh_mm':mesh['mesh_mm'],'dt_s':dt,'steps':len(times)-1,'elements':len(mesh['triangles']),'nodes':len(mesh['xyz']),
                  'wall_storage':'one midsurface temperature; rho*c*t per CAD area','normal_resistance_m2K_W':mesh['thickness_m']/p['box_k_W_mK'],
                  'surface_quadrature':'CAD face areas normalized exactly; mesh-independent mass and film integrals','time_integration':'backward Euler'},
        fluxes_20min_W=dict(zip(names,thermal['flux'][-1].tolist())),heat_transferred_J=energies,
        moisture={k:v for k,v in moisture.items() if k!='retention_history'},runtime_s=runtime,
        energy_balance_error_J=thermal.get('energy_balance_error_J'),status='temperature_pass' if passed else 'temperature_fail',
        assumptions=['Air retains exact 3D cavity volume and volumetric heat capacity, represented by the existing well-mixed lumped model, not CFD.',
                     'Pizza retains its volumetric mass and heat capacity; no surface-only conversion of air or pizza.',
                     'Uniform-thickness watertight STEP; open vents and unresolved overlapping layers are rejected.',
                     'Finite through-wall resistance with lumped midwall storage; sub-thickness transient gradients are reduced analytically.',
                     'Radiation, leakage, latent cooling and folding kinematics remain outside the original specified model.'],
        officiality='INDEPENDENT ESTIMATE, NOT AN OFFICIAL TA SCORE')
    (out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def run_shell(step,dxf,out,size_mm=3,dt=10,thickness_mm=3,dxf_unit='mm',max_elements=500000):
    from simulator import DEFAULTS
    start=time.perf_counter();p=DEFAULTS.copy();p['nominal_board_thickness_m']=thickness_mm*.001
    if not np.isfinite(dt) or dt<=0:raise ValueError('Time step must be finite and positive')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    mesh=build_shell(step,dxf,size_mm,p,dxf_unit,max_elements)
    selected=sensors(mesh)
    thermal=solve(mesh,p,dt,selected)
    np.savez_compressed(out/'thermal_field.npz',xyz_m=mesh['xyz'],triangles=mesh['triangles'],temperature_C=thermal['last_box_temp'],physical_area_m2=mesh['area'],thickness_m=mesh['thickness_m'])
    return report(mesh,p,thermal,out,step,dxf,time.perf_counter()-start,'Python thin-wall FEM')


def write_shell_apdl(mesh,out,p,dt,selected):
    sys=assemble(mesh,p);xyz=mesh['xyz'];tri=mesh['triangles'];n=len(xyz);ia=n+1;ip=n+2
    N=math.ceil(p['duration_s']/dt);dt=p['duration_s']/N
    factors,mat=np.unique(np.round(mesh['area_scale'],12),return_inverse=True)
    h=sys['h'];weights=sys['weights'];T0=p['T_environment_C']
    with (Path(out)/'run.inp').open('w',encoding='ascii') as f:
        def wr(*v):f.write(','.join(map(str,v))+'\n')
        f.write('/BATCH\n/CLEAR,START\n/FILNAME,PIZZA\n/PREP7\n')
        f.write('! Constant-temperature thermal shell; each face film includes t/(2k).\n')
        wr('ET',1,'SHELL131');wr('KEYOPT',1,3,2);wr('KEYOPT',1,15,1)
        wr('ET',2,'SURF152');wr('KEYOPT',2,4,1);wr('KEYOPT',2,5,1);wr('KEYOPT',2,8,2)
        wr('ET',3,'MASS71');wr('KEYOPT',3,3,1)
        wr('ET',4,'LINK33');wr('MP','KXX',10000,1)
        for j,scale in enumerate(factors,1):
            f.write('! CAD-area quadrature correction, not a change of physical thickness.\n')
            wr('MP','KXX',j,p['box_k_W_mK']*scale);wr('MP','DENS',j,p['box_density_kg_m3']*scale);wr('MP','C',j,p['box_cp_J_kgK'])
            wr('SECTYPE',j,'SHELL');wr('SECDATA',mesh['thickness_m'],j)
        for j,pos in enumerate(xyz,1):wr('N',j,*pos)
        wr('N',ia,-1,0,0);wr('N',ip,-2,0,0)
        wr('R',1,0);wr('R',2,sys['capacity'][-2]);wr('R',3,sys['capacity'][-1]);wr('R',4,sys['hpa'])
        for j,(q,m,scale) in enumerate(zip(tri,mat,mesh['area_scale']),1):
            wr('TYPE',1);wr('MAT',int(m)+1);wr('SECNUM',int(m)+1)
            ids=(q+1).tolist();wr('EN',j,*ids,ids[-1])
            he=h['floor'] if mesh['floor'][j-1] else h['exterior']
            wr('SFE',j,1,'CONV',0,he*scale);wr('SFE',j,1,'CONV',2,T0)
        eid=len(tri)
        for q,c,scale in zip(tri,mesh['contact'],mesh['area_scale']):
            wr('TYPE',2);wr('MAT',1);wr('REAL',1)
            ids=(q+1).tolist();eid+=1;wr('EN',eid,*ids,ids[-1],ip if c else ia)
            wr('SFE',eid,1,'CONV',0,(h['contact'] if c else h['air'])*scale)
        wr('TYPE',3);wr('REAL',2);eid+=1;wr('EN',eid,ia)
        wr('REAL',3);eid+=1;wr('EN',eid,ip)
        wr('TYPE',4);wr('MAT',10000);wr('REAL',4);eid+=1;wr('EN',eid,ia,ip)
        f.write('ALLSEL,ALL\nFINISH\n/SOLU\nANTYPE,TRANS\nTRNOPT,FULL\nAUTOTS,OFF\nKBC,1\nOUTRES,NSOL,ALL\n')
        wr('TUNIF',T0);wr('IC',ip,'TEMP',p['T_pizza_initial_C']);wr('TIME',p['duration_s']);wr('NSUBST',N,N,N)
        wr('DELTIM',dt,dt,dt);f.write('SOLVE\nFINISH\n/POST1\n')
        for name in ('TV','VV','WA','WC','WE','WF','CW'):wr('*DIM',name,'ARRAY',n)
        for j,value in enumerate(sys['capacity'][:n],1):
            f.write(f'CW({j})={value:.15g}\n')
        for name,label in [('WA','air'),('WC','contact'),('WE','exterior'),('WF','floor')]:
            for j,value in enumerate(weights[label],1):
                if value:f.write(f'{name}({j})={value:.15g}\n')
        # Fixed-width headers are literal formats; APDL character parameters truncate to 8 characters.
        f.write('*CFOPEN,history,csv\n*VWRITE\n(\'time_s,pizza_C,air_C\')\n*CFCLOS\n')
        f.write('*CFOPEN,wall_samples,csv\n*VWRITE\n(\'time_s,sensor_id,wall_C\')\n*CFCLOS\n')
        f.write('*CFOPEN,heat_flux_raw,csv\n*VWRITE\n(\'time_s,pa,pb,ab,out,floor,stored\')\n*CFCLOS\n')
        wr('*DO','j',1,N);f.write('SET,1,j\n');wr('*GET','clock','ACTIVE',0,'SET','TIME')
        wr('*GET','pzt','NODE',ip,'TEMP');wr('*GET','art','NODE',ia,'TEMP')
        f.write('*CFOPEN,history,csv,,APPEND\n*VWRITE,clock,pzt,art\n(F16.6,\',\',F16.9,\',\',F16.9)\n*CFCLOS\n')
        f.write('*CFOPEN,wall_samples,csv,,APPEND\n')
        factor=mesh['thickness_m']/(2*p['box_k_W_mK'])*h['air']
        for sample in selected:
            for k,node in enumerate(sample['nodes'],1):wr('*GET',f't{k}','NODE',node,'TEMP')
            f.write('tw=(t1+t2+t3)/3\n');f.write(f'tw=tw+{factor:.15g}*(art-tw)\n')
            wr('*VWRITE','clock',sample['id'],'tw');f.write('(F16.6,\',\',F8.0,\',\',F16.9)\n')
        f.write('*CFCLOS\n');wr('*VLEN',n);wr('*VGET','TV(1)','NODE',1,'TEMP')
        for var,weight,label,lump,sign in [('qpb','WC','contact','pzt',-1),('qab','WA','air','art',-1),('qout','WE','exterior',str(T0),1),('qfloor','WF','floor',str(T0),1)]:
            wr('*VOPER','VV(1)','TV(1)','MULT',f'{weight}(1)');wr('*VSCFUN','sumt','SUM','VV(1)')
            total=weights[label].sum()
            f.write(f'{var}={sign*h[label]:.15g}*(sumt-{total:.15g}*{lump})\n')
        wr('*VOPER','VV(1)','TV(1)','MULT','CW(1)');wr('*VSCFUN','sumenergy','SUM','VV(1)')
        f.write(f'energy=sumenergy-{sys["capacity"][:n].sum()*T0:.15g}+{sys["capacity"][-2]:.15g}*(art-{T0})+{sys["capacity"][-1]:.15g}*(pzt-{p["T_pizza_initial_C"]})\n')
        f.write(f'qpa={sys["hpa"]:.15g}*(pzt-art)\n*CFOPEN,heat_flux_raw,csv,,APPEND\n')
        wr('*VWRITE','clock','qpa','qpb','qab','qout','qfloor','energy');f.write('(F12.4,6(\',\',E17.9))\n*CFCLOS\n*ENDDO\nFINISH\n/EXIT,NOSAVE\n')
    return {'path':str(Path(out)/'run.inp'),'steps':N,'dt_s':dt,'mesh_mm':mesh['mesh_mm'],'elements':len(tri),'nodes':n,'total_ansys_elements':eid,'wall_thickness_m':mesh['thickness_m']}


def prepare_shell(step,dxf,out,size_mm,dt,samples=64,thickness_mm=3):
    from simulator import DEFAULTS
    start=time.perf_counter();p=DEFAULTS.copy();p['nominal_board_thickness_m']=thickness_mm*.001
    if dt<=0 or not np.isfinite(dt):raise ValueError('Time step must be finite and positive')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    mesh=build_shell(step,dxf,size_mm,p);selected=sensors(mesh,samples)
    deck=write_shell_apdl(mesh,out,p,dt,selected)
    np.savez_compressed(out/'shell_mesh.npz',**{k:mesh[k] for k in ['xyz','triangles','area','area_scale','floor','contact']})
    meta=dict(kind='pizza-contest-ansys-mapdl',model_version='thinwall-v1',mesh_mode='shell',params=p,step=str(step),dxf=str(dxf),mesh=deck,
              geometry=geometry_metrics(mesh),flat=mesh['flat'],mass_g=mesh['flat']['mass_g'],sensors=selected,
              preparation_runtime_s=time.perf_counter()-start,status='prepared_not_run')
    (out/'metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    with (out/'wall_patches.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['sensor_id','area_m2','node1','node2','node3'])
        for s in selected:w.writerow([s['id'],s['area_m2'],*s['nodes']])
    return meta


def finalize_shell(out,meta):
    from ansys_bridge import parse_result_csv
    out=Path(out)
    import re
    log=(out/'ansys.out').read_text(errors='replace')
    if re.search(r'\*\*\*\s+(ERROR|FATAL)\s+\*\*\*',log) or 'RUN COMPLETED' not in log:
        raise RuntimeError('ANSYS shell log contains errors or lacks RUN COMPLETED')
    p=meta['params'];h=parse_result_csv(out/'history.csv',3);w=parse_result_csv(out/'wall_samples.csv',3);flows=parse_result_csv(out/'heat_flux_raw.csv',7)
    N=meta['mesh']['steps'];dt=meta['mesh']['dt_s'];ns=len(meta['sensors'])
    if len(h)!=N or len(w)!=N*ns or len(flows)!=N:raise RuntimeError('Incomplete shell ANSYS outputs')
    expected=np.arange(1,N+1)*dt
    w=w.reshape(N,ns,3)
    if not np.allclose(h[:,0],expected) or not np.allclose(flows[:,0],expected) or not np.allclose(w[:,:,0],expected[:,None]) or not np.all(w[:,:,1]==np.arange(1,ns+1)[None,:]):
        raise RuntimeError('ANSYS shell time grid or sensor IDs do not match metadata')
    mesh=dict(np.load(out/'shell_mesh.npz'));mesh.update(meta['geometry']);mesh['flat']=meta['flat'];mesh['mesh_mm']=meta['mesh']['mesh_mm']
    T0=p['T_environment_C'];init=np.full(len(mesh['xyz'])+2,T0);init[-1]=p['T_pizza_initial_C']
    initial_flux=flux(init,assemble(mesh,p),p)
    thermal=dict(time=np.r_[0,h[:,0]],Tp=np.r_[p['T_pizza_initial_C'],h[:,1]],Ta=np.r_[T0,h[:,2]],wallT=np.vstack((np.full(ns,T0),w[:,:,2])),
                 wallareas=np.array([s['area_m2'] for s in meta['sensors']]),flux=np.vstack((initial_flux,flows[:,1:6])),energy_balance_error_J=float(flows[-1,6]+dt*flows[:,4:6].sum()))
    if abs(thermal['energy_balance_error_J'])>0.001*p['pizza_mass_kg']*p['pizza_cp_J_kgK']*(p['T_pizza_initial_C']-T0):
        raise RuntimeError('ANSYS shell energy balance error exceeds 0.1%')
    runtime=float((out/'runtime_s.txt').read_text())
    if not np.isfinite(runtime) or runtime<0:raise RuntimeError('Invalid ANSYS runtime')
    return report(mesh,p,thermal,out,meta['step'],meta['dxf'],runtime,'ANSYS SHELL131 thin-wall FEM')
