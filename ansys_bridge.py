"""Windows ANSYS MAPDL batch + Tailscale SSH backend for pizza-box contest.

The actual mathematical model is kept identical between solvers: 3D transient
solid conduction, internal/external convection, thermal contact conductance,
and two lumped thermal masses.  MAPDL is a separate licensed solver; none of
these methods imply ANSYS is installed on the computer generating the deck.

Python dependencies: numpy, scipy, cadquery, ezdxf, shapely, matplotlib.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import zipfile
import numpy as np
from scipy.spatial import cKDTree

from simulator import (DEFAULTS, FACE, GeometryError, extract_flat, inspect_step,
                       voxelize, label_air, validate_pizza_fit, build_mesh,
                       p_sat_Pa)

# HEX8 node mapping: local indexing in simulator is (z, y, x) with x fastest.
ANSYS_HEX = np.array([0,1,3,2,4,5,7,6])
# Faces have cyclic ordering for valid SURF152 quadrilaterals.
ANSYS_FACE = ((0,1,3,2),(4,5,7,6),(0,1,5,4),(2,3,7,6),(0,2,6,4),(1,3,7,5))


def validate_common(step: Path, dxf: Path):
    flat=extract_flat(dxf)
    shape,bbox,mass_step=inspect_step(step)
    mismatch=abs(mass_step-flat['mass_g']) / flat['mass_g']
    if mismatch>0.12:
        raise GeometryError(f'DXF/STEP mass mismatch {mismatch:.1%}; maximum 12%.')
    if not DEFAULTS['min_box_mass_g'] <= flat['mass_g'] <= DEFAULTS['max_box_mass_g']:
        raise GeometryError(f'Box mass {flat["mass_g"]:.1f} g is outside 83-249 g')
    shells=shape.Shells()
    if len(shells)!=2:
        raise GeometryError(f'STEP needs exactly two closed shells (outside + one sealed cavity); got {len(shells)}')
    volumes=[]
    for shell in shells:
        s=__import__('cadquery').Solid.makeSolid(shell)
        if not s.isValid():
            raise GeometryError('A cavity/surface shell cannot form a valid closed solid.')
        volumes.append(s.Volume()*1e-9)
    cavity_volume=min(volumes)
    pizza_volume=math.pi*(DEFAULTS['pizza_diameter_m']/2)**2*DEFAULTS['pizza_thickness_m']
    air_volume=cavity_volume-pizza_volume
    if air_volume<=0:raise GeometryError('Insufficient internal cavity for pizza and air.')
    return shape,bbox,flat,mass_step,cavity_volume,air_volume


def make_voxel_mesh(shape,bbox,voxel_mm,cavity_volume,params):
    occ,spacing,nxyz=voxelize(shape,bbox,voxel_mm*0.001,max_cells=2_500_000)
    lab,ext,cavity=label_air(occ)
    fit=validate_pizza_fit(occ,lab,cavity,bbox,spacing,params)
    cell_idx,elements,node_idx,_=build_mesh(occ,nxyz)
    from scipy import ndimage
    _,ncomponents=ndimage.label(occ,ndimage.generate_binary_structure(3,1))
    if ncomponents!=1:raise GeometryError(f'Voxel FEM disconnected into {ncomponents} components.')
    dx,dy,dz=spacing
    xyz=bbox[:,0] + node_idx*np.array([dx,dy,dz])[None,:]
    bottom=[]; exterior=[]; interior=[]; contact=[]
    cx=(bbox[0,0]+bbox[0,1])/2;cy=(bbox[1,0]+bbox[1,1])/2
    r=params['pizza_diameter_m']/2
    face_area=np.array([dx*dy,dx*dy,dx*dz,dx*dz,dy*dz,dy*dz])
    for side,(direction,_,_) in enumerate(FACE):
        iz,iy,ix=cell_idx.T
        region=lab[iz+1+direction[0],iy+1+direction[1],ix+1+direction[2]]
        qnodes=elements[:,ANSYS_FACE[side]]
        exids=np.flatnonzero(region==ext)
        inids=np.flatnonzero(region==cavity)
        area=face_area[side]
        if len(exids):
            if side==0:bottom.extend(qnodes[exids].tolist())
            else:
                for nodes in qnodes[exids]: exterior.append((nodes,area))
        if len(inids):
            nodes=qnodes[inids]
            if side==1:
                xc=bbox[0,0]+(cell_idx[inids,2]+.5)*dx
                yc=bbox[1,0]+(cell_idx[inids,1]+.5)*dy
                touches=(xc-cx)**2+(yc-cy)**2<r**2
            else: touches=np.zeros(len(inids),dtype=bool)
            for q in nodes[touches]:contact.append((q,area))
            for q in nodes[~touches]:interior.append((q,area))
    if not bottom or not contact or not interior:
        raise GeometryError('No exterior floor or pizza-floor contact or internal air faces')
    ca=sum(a for _,a in contact)
    if abs(ca/(math.pi*r*r)-1)>.12:
        raise GeometryError(f'Contact area {ca:.6f} m² differs from pizza disk by more than 12%')
    mesh={
        'xyz':xyz,'solid_elements':elements[:,ANSYS_HEX], 'solid_type':'HEX8',
        'bottom':np.unique(np.array(bottom,dtype=np.int64)).tolist(),
        'exterior':exterior,'interior':interior,'contact':contact,
        'cavity_volume_m3':cavity_volume,'contact_area_m2':ca,
        'pizza_fit':fit,'mesh_mm':voxel_mm,
        'voxel_board_volume_m3':float(occ.sum()*dx*dy*dz),
        'warnings':['Voxel staircase approximates true STEP boundary; perform mesh convergence.',
                    'SURF152 film effects use a lumped air node and area-averaged conditions.'],
    }
    return mesh


def make_gmsh_tetra_mesh(step:Path,bbox,cavity_volume,params,mesh_mm):
    """Optional exact STEP-boundary tetrahedral meshing through Gmsh/OpenCascade.

    Gmsh is optional and installed only on the machine preparing this input.
    Mesh generation is NOT required on the Windows ANSYS host.
    Uses SOLID87 TET10 and linear triangular SURF152 overlays on tetra faces.
    """
    try:import gmsh
    except ImportError as exc:
        raise RuntimeError('Exact STEP tessellation requires pip install gmsh and the Gmsh runtime. '
                           'Or use --mesh-mode voxel.') from exc
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.model.add('pizza_box_step')
        tags=gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        volumes=gmsh.model.getEntities(3)
        if len(volumes)!=1:
            raise GeometryError(f'Gmsh STEP import found {len(volumes)} volumes, expected exactly one')
        gmsh.option.setNumber('Mesh.MeshSizeMin',mesh_mm*.65)
        gmsh.option.setNumber('Mesh.MeshSizeMax',mesh_mm)
        gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',10)
        gmsh.option.setNumber('Mesh.Algorithm3D',1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        nodetags, nodecoords,_=gmsh.model.mesh.getNodes()
        nodetags=np.asarray(nodetags,dtype=np.int64)
        xyz=np.asarray(nodecoords,dtype=np.float64).reshape(-1,3)*.001
        # Gmsh node tags are not guaranteed continuous.
        perm=np.argsort(nodetags)
        sortedtags=nodetags[perm]
        def remap(tags):
            loc=np.searchsorted(sortedtags,tags)
            if np.any(loc>=len(sortedtags)) or not np.all(sortedtags[np.minimum(loc,len(sortedtags)-1)]==tags):
                raise GeometryError('Gmsh node tag missing from mesh node coordinates')
            return perm[loc]
        element_types,_,conn=gmsh.model.mesh.getElements(3)
        if 11 not in element_types:
            raise GeometryError('No 10-node tetrahedra found after mesh.setOrder(2)')
        tet=np.asarray(conn[list(element_types).index(11)],dtype=np.int64).reshape(-1,10)
        # Gmsh TET10 last two mid-edge indices differ from ANSYS SOLID87.
        # Gmsh: (0-1,1-2,2-0,0-3,2-3,1-3); ANSYS: (I-J,J-K,K-I,I-L,J-L,K-L).
        tet=tet[:,[0,1,2,3,4,5,6,7,9,8]]
        tet=remap(tet)
        ftypes,_,fconn=gmsh.model.mesh.getElements(2)
        allfaces=[]
        for typ,block in zip(ftypes,fconn):
            if typ==9:allfaces.append(remap(np.asarray(block,dtype=np.int64).reshape(-1,6)[:,:3]))
            if typ==2:allfaces.append(remap(np.asarray(block,dtype=np.int64).reshape(-1,3)))
        if not allfaces:raise GeometryError('Gmsh found no surface triangle mesh')
        tris=np.vstack(allfaces).astype(np.int64)
        # Surface components: exterior shell and one internal cavity shell.
        # Each triangle shares corner nodes along edges with its component.
        allcorner=np.unique(tris)
        compact=np.searchsorted(allcorner,tris)
        rr=np.repeat(compact,3,axis=1).ravel()
        cc=np.concatenate((compact[:,1],compact[:,2],compact[:,0]),axis=1).ravel() if False else np.column_stack((compact[:,1],compact[:,2],compact[:,0])).ravel()
        g=csr_matrix((np.ones(len(rr),dtype=np.uint8),(rr,cc)),shape=(len(allcorner),len(allcorner)))
        ncc,labels=connected_components(g,directed=False)
        tri_comp=labels[compact[:,0]]
        if ncc!=2:
            raise GeometryError(f'Gmsh surface has {ncc} components, expected outside + one sealed cavity')
        # Bounding box of exterior shell matches STEP solid outside bounding box.
        comp_box=[]
        for k in range(ncc):
            pts=xyz[np.unique(tris[tri_comp==k].ravel())]
            dist=float(np.max(abs(pts.min(axis=0)-bbox[:,0])) + np.max(abs(pts.max(axis=0)-bbox[:,1])))
            comp_box.append(dist)
        extcomp=int(np.argmin(comp_box)); innercomp=1-extcomp
        if comp_box[extcomp]>mesh_mm*0.001*1.5:
            raise GeometryError('Could not confidently classify outside boundary surface')
        a=xyz[tris[:,0]]; b=xyz[tris[:,1]]; c=xyz[tris[:,2]]
        area=.5*np.linalg.norm(np.cross(b-a,c-a),axis=1)
        centers=(a+b+c)/3
        bottom=[];exterior=[];interior=[];contact=[]
        zmin=bbox[2,0];zbase=float(np.min(centers[tri_comp==innercomp,2]))
        r=params['pizza_diameter_m']/2
        cx=(bbox[0,0]+bbox[0,1])/2;cy=(bbox[1,0]+bbox[1,1])/2
        for j,q in enumerate(tris):
            areaj=float(area[j]); center=centers[j]
            if tri_comp[j]==extcomp:
                if np.max(np.abs(xyz[q,2]-zmin))<=.00015:bottom.extend(q.tolist())
                else:exterior.append((q,areaj))
            else:
                floor=(np.max(abs(xyz[q,2]-zbase))<.00015)
                if floor and (center[0]-cx)**2+(center[1]-cy)**2<r*r:
                    contact.append((q,areaj))
                else:interior.append((q,areaj))
        ca=sum(ar for _,ar in contact)
        if len(bottom)==0 or not contact or not interior:
            raise GeometryError('Incomplete tetra boundary face classification')
        if abs(ca/(math.pi*r*r)-1)>.10:
            raise GeometryError(f'Gmsh pizza contact area mismatch: {ca:.5f} m²')
        bbox_inner_max=np.max(centers[tri_comp==innercomp,2])
        if bbox_inner_max-zbase<params['pizza_thickness_m']:
            raise GeometryError('Pizza does not fit inside cavity height')
        return {'xyz':xyz,'solid_elements':tet,'solid_type':'TET10',
                'bottom':np.unique(bottom).tolist(), 'exterior':exterior,
                'interior':interior,'contact':contact,
                'cavity_volume_m3':cavity_volume,'contact_area_m2':ca,
                'pizza_fit':{'inner_floor_z_mm':zbase*1000,'inner_ceiling_z_mm':bbox_inner_max*1000},
                'mesh_mm':mesh_mm,'voxel_board_volume_m3':None,
                'warnings':['Quadratic TET10 STEP-conforming volume mesh; linearized convection surfaces.',
                            'Pizza-floor contact uses triangle centroid selection; verify mesh convergence.']}
    finally:
        gmsh.finalize()


def cluster_sensors(mesh,max_samples=64,seed=1234):
    """Deterministically sample wall patches; preserve total internal convection area."""
    walls=mesh['interior']; xyz=mesh['xyz']
    ctr=np.array([xyz[np.asarray(nodes)].mean(axis=0) for nodes,_ in walls])
    areas=np.array([v for _,v in walls])
    n=min(max_samples,len(walls))
    # Surface area weighted deterministic sampling with feature coverage.
    rng=np.random.default_rng(seed)
    # Split into 6 axis-based groups to ensure thin side walls are represented.
    sampled=[]
    centers=ctr
    for ax in range(3):
        for sign in (-1,1):
            norm=(centers-centers.mean(axis=0))[:,ax]
            pool=np.where(norm*sign>0)[0]
            if len(pool):
                sampled.append(int(pool[np.argmax(areas[pool])]))
    sampled=list(dict.fromkeys(sampled))
    p=areas/areas.sum()
    while len(sampled)<n:
        choices=rng.choice(len(walls),size=min(len(walls),max(2*n,80)),replace=True,p=p)
        if sampled:
            dist=cKDTree(centers[sampled]).query(centers[choices])[0]
            ix=int(choices[np.argmax(dist)])
        else:ix=int(choices[0])
        if ix not in sampled:sampled.append(ix)
        elif len(sampled)>=len(walls):break
    sampled=np.array(sampled,dtype=int)
    chosen=centers[sampled]
    region=cKDTree(chosen).query(centers)[1]
    patchareas=np.bincount(region,weights=areas,minlength=len(sampled))
    sensors=[{'id':i+1,'nodes':[int(k)+1 for k in np.asarray(walls[int(idx)][0])],
              'area_m2':float(patchareas[i])}
             for i,idx in enumerate(sampled)]
    assert np.isclose(sum(s['area_m2'] for s in sensors),areas.sum())
    return sensors


def write_apdl(mesh,outdir:Path, air_volume:float, params:dict, dt:float, sensors:list):
    """Emit a single native MAPDL APDL deck. All geometry units are SI."""
    xyz=mesh['xyz'];els=mesh['solid_elements']
    nnode=len(xyz);nelem=len(els)
    airnode=nnode+1;pizzanode=nnode+2
    x0=float(np.min(xyz[:,0]));yc=float(np.mean([xyz[:,1].min(),xyz[:,1].max()]));
    zc=float(np.mean([xyz[:,2].min(),xyz[:,2].max()]))
    N=max(1,math.ceil(params['duration_s']/dt)); dt=params['duration_s']/N
    pizza_area=math.pi*(params['pizza_diameter_m']/2)**2 + math.pi*params['pizza_diameter_m']*params['pizza_thickness_m']
    cap_air=params['air_density_kg_m3']*params['air_cp_J_kgK']*air_volume
    cap_pizza=params['pizza_mass_kg']*params['pizza_cp_J_kgK']
    hpa=params['h_internal_W_m2K']*pizza_area
    path=outdir/'run.inp'
    with open(path,'w',encoding='ascii',newline='\n') as f:
        def wr(*parts):
            f.write(','.join(str(k) for k in parts)+'\n')
        f.write('! Pizza contest - generated by ansys_bridge.py\n')
        f.write('! SI geometry/material, temperatures in degrees Celsius\n')
        f.write('/BATCH\n/CLEAR,START\n/FILNAME,PIZZA\n/PREP7\n')
        solidtype='SOLID70' if mesh['solid_type']=='HEX8' else 'SOLID87'
        wr('ET',1,solidtype)
        wr('MP','KXX',1,params['box_k_W_mK'])
        wr('MP','C',1,params['box_cp_J_kgK'])
        wr('MP','DENS',1,params['box_density_kg_m3'])
        # Surface types: exterior known 22C; internal unknown bulk extra node (air/pizza).
        wr('ET',2,'SURF152');wr('KEYOPT',2,4,1);wr('KEYOPT',2,5,0);wr('KEYOPT',2,8,2)
        wr('ET',3,'SURF152');wr('KEYOPT',3,4,1);wr('KEYOPT',3,5,1);wr('KEYOPT',3,8,2)
        wr('ET',4,'MASS71');wr('KEYOPT',4,3,1)
        wr('ET',5,'LINK33');wr('MP','KXX',5,1)
        wr('R',1,0)
        wr('R',4,cap_air);wr('R',5,cap_pizza)
        wr('R',6,hpa)  # Link length = 1 m, k = 1 W/(m.K); R=conductance in m².
        for j,pos in enumerate(xyz,start=1):
            wr('N',j,f'{pos[0]:.11g}',f'{pos[1]:.11g}',f'{pos[2]:.11g}')
        wr('N',airnode,f'{x0-1:.11g}',f'{yc:.11g}',f'{zc:.11g}')
        wr('N',pizzanode,f'{x0-2:.11g}',f'{yc:.11g}',f'{zc:.11g}')
        f.write('TYPE,1\nMAT,1\nREAL,1\n')
        for j,nd in enumerate(els,start=1):
            if mesh['solid_type']=='HEX8':wr('EN',j,*(int(n)+1 for n in nd))
            else:
                wr('EN',j,*(int(n)+1 for n in nd[:8]))
                wr('EMORE',*(int(n)+1 for n in nd[8:]))
        eid=nelem
        f.write('TYPE,4\nREAL,4\n');eid+=1;wr('EN',eid,airnode)
        f.write('TYPE,4\nREAL,5\n');eid+=1;wr('EN',eid,pizzanode)
        f.write('TYPE,5\nMAT,5\nREAL,6\n');eid+=1;wr('EN',eid,airnode,pizzanode)
        # SOLID70 and SOLID87 use quadratic/linear-volume conduction; surfaces get one extra node.
        f.write('TYPE,2\nMAT,1\nREAL,1\n')
        for nd,area in mesh['exterior']:
            eid+=1
            ids=[int(n)+1 for n in nd]
            if len(ids)==3:ids.append(ids[-1])
            wr('EN',eid,*ids)
            wr('SFE',eid,1,'CONV',0,params['h_external_W_m2K'])
            wr('SFE',eid,1,'CONV',2,params['T_environment_C'])
        f.write('TYPE,3\nMAT,1\nREAL,1\n')
        for faces,node,h in ((mesh['interior'],airnode,params['h_internal_W_m2K']),
                             (mesh['contact'],pizzanode,params['h_contact_W_m2K'])):
            for nd,area in faces:
                eid+=1
                ids=[int(n)+1 for n in nd]
                if len(ids)==3:ids.append(ids[-1])
                wr('EN',eid,*ids,node)
                wr('SFE',eid,1,'CONV',0,h)
        f.write('ALLSEL,ALL\n')
        for node in mesh['bottom']:wr('D',int(node)+1,'TEMP',params['T_environment_C'])
        f.write('FINISH\n/SOLU\nANTYPE,TRANS\nTRNOPT,FULL\n')
        f.write('AUTOTS,OFF\nKBC,1\nOUTRES,NSOL,ALL\n')
        wr('TUNIF',params['T_environment_C'])
        wr('IC',pizzanode,'TEMP',params['T_pizza_initial_C'])
        wr('TIME',params['duration_s'])
        wr('NSUBST',N,N,N)
        wr('DELTIM',dt,dt,dt)
        f.write('SOLVE\nFINISH\n/POST1\n')
        wr('*DIM','SN1','ARRAY',len(sensors))
        wr('*DIM','SN2','ARRAY',len(sensors))
        wr('*DIM','SN3','ARRAY',len(sensors))
        wr('*DIM','SN4','ARRAY',len(sensors))
        for item in sensors:
            i=item['id'];nodes=item['nodes']
            if len(nodes)==3:nodes.append(nodes[-1])
            for k,nn in enumerate(nodes,start=1):f.write(f'SN{k}({i})={int(nn)}\n')
        f.write('*CFOPEN,history,csv\n')
        f.write("*VWRITE,'time_s,pizza_C,air_C'\n(A)\n")
        wr('*DO','j',1,int(N))
        f.write('SET,1,j\n')
        wr('*GET','clock','ACTIVE',0,'SET','TIME')
        wr('*GET','pzt','NODE',pizzanode,'TEMP')
        wr('*GET','art','NODE',airnode,'TEMP')
        f.write('*VWRITE,clock,pzt,art\n(F16.6,\',\',F16.9,\',\',F16.9)\n')
        f.write('*ENDDO\n*CFCLOS\n')
        f.write('*CFOPEN,wall_samples,csv\n')
        f.write("*VWRITE,'time_s,sensor_id,wall_C'\n(A)\n")
        wr('*DO','j',1,int(N))
        f.write('SET,1,j\n')
        wr('*GET','clock','ACTIVE',0,'SET','TIME')
        wr('*DO','k',1,len(sensors))
        f.write('n1=SN1(k)\nn2=SN2(k)\nn3=SN3(k)\nn4=SN4(k)\n')
        f.write('*GET,tt1,NODE,n1,TEMP\n*GET,tt2,NODE,n2,TEMP\n')
        f.write('*GET,tt3,NODE,n3,TEMP\n*GET,tt4,NODE,n4,TEMP\n')
        f.write('meanT=(tt1+tt2+tt3)/3\n' if mesh['solid_type']=='TET10' else 'meanT=(tt1+tt2+tt3+tt4)/4\n')
        f.write('*VWRITE,clock,k,meanT\n(F16.6,\',\',F6.0,\',\',F16.9)\n')
        f.write('*ENDDO\n*ENDDO\n*CFCLOS\n')
        f.write('FINISH\n/EXIT,NOSAVE\n')
    return {'path':str(path),'steps':N,'dt_s':dt,'mesh_mm':mesh['mesh_mm'],'elements':nelem,
            'nodes':nnode,'surface_elements':eid-nelem-3,'pizza_node':pizzanode,
            'air_node':airnode,'pizza_air_conductance_W_K':hpa}


PS1 = r'''param([string]$AnsysExe = $env:ANSYS_EXE, [int]$Nproc = 4)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
if (-not $AnsysExe) {
  $base = Join-Path ${env:ProgramFiles} 'ANSYS Inc'
  if (-not (Test-Path $base)) { throw "ANSYS directory not found: $base. Set ANSYS_EXE." }
  $roots = @($base, (Join-Path $base 'ANSYS Student'))
  $candidates = foreach ($root in $roots) {
    Get-ChildItem -LiteralPath $root -Directory -ErrorAction SilentlyContinue | Where-Object { $_.Name -match '^v\d+' }
  }
  $candidates = $candidates | Sort-Object Name -Descending
  foreach ($v in $candidates) {
    $exes = Get-ChildItem -Path (Join-Path $v.FullName 'ANSYS\bin\winx64') -Filter 'ansys*.exe' -ErrorAction SilentlyContinue
    if ($exes) { $AnsysExe = ($exes | Sort-Object Name -Descending | Select-Object -First 1).FullName; break }
  }
}
if (-not $AnsysExe -or -not (Test-Path -LiteralPath $AnsysExe)) { throw "ANSYS executable unavailable: $AnsysExe" }
Write-Host "Running $AnsysExe"
$env:ANS_CONSEC = 'YES'
& $AnsysExe -b -np $Nproc -j PIZZA -i run.inp -o ansys.out
$exitCode = $LASTEXITCODE
if ($exitCode -ne 0) { throw "ANSYS process exited $exitCode. Inspect ansys.out" }
if (-not (Test-Path -LiteralPath history.csv)) { throw 'Missing history.csv - ANSYS did not finish post-processing.' }
if (-not (Test-Path -LiteralPath wall_samples.csv)) { throw 'Missing wall_samples.csv - ANSYS did not finish post-processing.' }
Write-Host 'ANSYS_RESULT_READY'
'''


def prepare(step,dxf,out,voxel_mm=1.8,dt=10.0,samples=64,mesh_mode='voxel'):
    step=Path(step).resolve();dxf=Path(dxf).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    params=DEFAULTS.copy()
    shape,bbox,flat,step_mass,Vc,Va=validate_common(step,dxf)
    if mesh_mode=='voxel':mesh=make_voxel_mesh(shape,bbox,voxel_mm,Vc,params)
    elif mesh_mode=='step':mesh=make_gmsh_tetra_mesh(step,bbox,Vc,params,voxel_mm)
    else:raise ValueError('Unknown mesh_mode')
    selected=cluster_sensors(mesh,max_samples=samples)
    deck=write_apdl(mesh,out,Va,params,dt,selected)
    meta={'kind':'pizza-contest-ansys-mapdl','status':'prepared_not_run',
          'officiality':'INDEPENDENT ESTIMATE - NOT TA RESULT',
          'mesh_mode':mesh_mode,'params':params,'dxf':str(dxf),'step':str(step),
          'mass_g':flat['mass_g'],'step_mass_g':step_mass,
          'cavity_volume_m3':Vc,'free_air_volume_m3':Va,
          'pizza_area_m2':math.pi*(params['pizza_diameter_m']/2)**2+math.pi*params['pizza_diameter_m']*params['pizza_thickness_m'],
          'contact_area_m2':mesh['contact_area_m2'],
          'inner_area_m2':sum(a for _,a in mesh['interior']),
          'mesh':deck,'mesh_warnings':mesh['warnings'],
          'pizza_fit':mesh['pizza_fit'],'sensors':selected,'flat':flat}
    (out/'metadata.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    (out/'run_ansys.ps1').write_text(PS1,encoding='utf-8')
    with open(out/'wall_patches.csv','w',newline='') as f:
        w=csv.writer(f);w.writerow(['sensor_id','area_m2','node1','node2','node3','node4'])
        for s in selected:w.writerow([s['id'],s['area_m2'],*s['nodes']])
    return meta


def parse_result_csv(path,ncols):
    data=np.genfromtxt(path,delimiter=',',skip_header=1)
    if data.size==0:raise RuntimeError(f'Empty ANSYS result: {path}')
    if data.ndim==1:data=data.reshape(1,-1)
    if data.shape[1]!=ncols or not np.all(np.isfinite(data)):
        raise RuntimeError(f'Invalid column count or nonfinite ANSYS results in {path}')
    return data


def finalize(folder):
    out=Path(folder).resolve()
    meta=json.loads((out/'metadata.json').read_text(encoding='utf-8'))
    h=parse_result_csv(out/'history.csv',3)
    w=parse_result_csv(out/'wall_samples.csv',3)
    nt=meta['mesh']['steps'];dt=meta['mesh']['dt_s']; ns=len(meta['sensors'])
    if len(h)!=nt or len(w)!=nt*ns:
        raise RuntimeError(f'Incomplete MAPDL results: expected {nt} history x {ns} sensors, '
                           f'got {len(h)} history and {len(w)} wall rows')
    if abs(h[-1,0]-1200)>dt*.01 or abs(h[0,0]-dt)>dt*.01:
        raise RuntimeError('ANSYS time stamps do not match the 20-minute run.')
    w2=w.reshape(nt,ns,3)
    ids=np.arange(1,ns+1)
    if not np.all(np.rint(w2[:,:,1]).astype(int)==ids[None,:]):
        raise RuntimeError('ANSYS wall sensor ordering mismatch')
    if not np.allclose(w2[:,:,0],h[:,0,None],atol=1e-2):
        raise RuntimeError('ANSYS wall time stamps mismatch')
    # Append initial 22C wall + air and 90C pizza, then solve moisture with BDF.
    from scipy.integrate import solve_ivp
    p=meta['params']; ts=np.r_[0,h[:,0]]
    tp=np.r_[p['T_pizza_initial_C'],h[:,1]]
    ta=np.r_[p['T_environment_C'],h[:,2]]
    wall=np.vstack((np.full(ns,p['T_environment_C']),w2[:,:,2]))
    area=np.array([s['area_m2'] for s in meta['sensors']])
    vp=meta['free_air_volume_m3'];rv=p['Rv_J_kgK'];km=p['km_m_s']
    a_p=meta['pizza_area_m2']
    mw0=p['pizza_mass_kg']*p['water_fraction_initial']
    mv0=p['RH_initial']*float(p_sat_Pa(p['T_environment_C']))*vp/(rv*(p['T_environment_C']+273.15))
    def moisture(t,y):
        j=max(0,min(len(ts)-2,int(np.searchsorted(ts,t,side='right'))-1))
        alpha=float(np.clip((t-ts[j])/(ts[j+1]-ts[j]),0.,1.))
        Tpizza=(1-alpha)*tp[j]+alpha*tp[j+1]
        Tair=(1-alpha)*ta[j]+alpha*ta[j+1]
        Tw=(1-alpha)*wall[j]+alpha*wall[j+1]
        pv=max(0,float(y[0]))*rv*(Tair+273.15)/vp
        evap=max(0.,km*a_p*(p['water_activity']*float(p_sat_Pa(Tpizza))-pv)/(rv*((Tpizza+Tair)/2+273.15)))
        cond=km*area*np.maximum(pv-p_sat_Pa(Tw),0.)/(rv*((Tw+Tair)/2+273.15))
        return [evap-float(np.sum(cond)),-evap,float(np.sum(cond))]
    res=solve_ivp(moisture,(0.,1200.),[mv0,mw0,0.],method='BDF',t_eval=ts,
                  rtol=1e-7,atol=1e-10,max_step=10)
    if not res.success:raise RuntimeError(res.message)
    vapor,water,condensate=res.y
    if abs(vapor[-1]+water[-1]+condensate[-1]-(mv0+mw0))>1e-7:
        raise RuntimeError('Water mass conservation check failed')
    retention=100.*water/mw0
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5))
    ax.plot(ts/60,tp,label='Pizza');ax.plot(ts/60,ta,label='Air')
    ax.axhline(68,color='grey',linestyle='--',label='68 C threshold')
    ax.set(xlabel='Time (minutes)',ylabel='Temperature (C)',title='ANSYS MAPDL transient FEM')
    ax.legend();ax.grid(alpha=.3);fig.tight_layout();fig.savefig(out/'temperature.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5));ax.plot(ts/60,retention)
    ax.set(xlabel='Time (minutes)',ylabel='Pizza water retained (%)',title='Sequential moisture')
    ax.grid(alpha=.3);fig.tight_layout();fig.savefig(out/'moisture.png',dpi=160);plt.close(fig)
    with open(out/'moisture.csv','w',newline='') as f:
        wri=csv.writer(f);wri.writerow(['time_s','pizza_water_g','air_vapor_g','wall_water_g','retention_pct'])
        wri.writerows(zip(ts,water*1000,vapor*1000,condensate*1000,retention))
    official={'mass_g':meta['mass_g'],'pizza_temperature_20min_C':float(tp[-1]),
              'pizza_moisture_retention_20min_pct':float(retention[-1]),
              'temperature_threshold_pass':bool(tp[-1]>=68)}
    result={'official_metrics':official,'solver':'ANSYS MAPDL full transient 3D FEM + lumped masses + sequential vapor model',
            'geometry_mesh_mode':meta['mesh_mode'], 'geometry':{'step':meta['step'],'dxf':meta['dxf'],
            'mass_step_g':meta['step_mass_g'],'box_air_volume_m3':vp},
            'numerics':meta['mesh'],'warnings':meta['mesh_warnings'],
            'moisture':{'initial_pizza_water_g':mw0*1000,'final_pizza_water_g':float(water[-1]*1000),
                        'wall_condensate_g':float(condensate[-1]*1000),
                        'air_vapor_g':float(vapor[-1]*1000)},
            'status':'temperature_pass' if tp[-1]>=68 else 'temperature_fail',
            'officiality':'INDEPENDENT ESTIMATE, NOT AN OFFICIAL TA SCORE',
            'validation':'Deck generated and postprocessed; requires independent ANSYS solver check and mesh/time convergence.'}
    (out/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    return result


def remote_run(jobdir,host,remote_base='pizza_box_contest/jobs',nproc=4,ansys_exe='',dry_run=False):
    """Uses system OpenSSH/scp; never stores passwords. Windows OpenSSH required."""
    jobdir=Path(jobdir).resolve();host=str(host).strip()
    if not re.fullmatch(r'[a-zA-Z0-9_.@-]+',host):
        raise ValueError('Invalid --host; use sshuser@100.x.y.z or a host alias')
    if not re.fullmatch(r'[a-zA-Z0-9_./-]+',remote_base):
        raise ValueError('Invalid --remote-base path')
    if not (jobdir/'run.inp').exists():raise FileNotFoundError(jobdir/'run.inp')
    name=re.sub(r'[^a-zA-Z0-9_-]','_',jobdir.name)
    remotedir=f'{remote_base.rstrip("/")}/{name}'
    required=['run.inp','metadata.json','run_ansys.ps1','wall_patches.csv']
    for file in required:
        if not (jobdir/file).exists():raise FileNotFoundError(jobdir/file)
    # Windows SSH's login working directory is usually %USERPROFILE%.
    remote_create=f'powershell -NoProfile -Command "New-Item -ItemType Directory -Force -Path {remotedir} | Out-Null"'
    psargs=f'-Nproc {int(nproc)}'
    if ansys_exe:
        if '"' in ansys_exe or '\n' in ansys_exe:raise ValueError('Invalid ANSYS executable path')
        psargs+=f' -AnsysExe "{ansys_exe}"'
    remote_exec=f'powershell -NoProfile -ExecutionPolicy Bypass -File {remotedir}/run_ansys.ps1 {psargs}'
    commands=[['ssh',host,remote_create],
              ['scp',*(str(jobdir/x) for x in required),f'{host}:{remotedir}/'],
              ['ssh',host,remote_exec],
              ['scp',f'{host}:{remotedir}/history.csv',str(jobdir/'history.csv')],
              ['scp',f'{host}:{remotedir}/wall_samples.csv',str(jobdir/'wall_samples.csv')],
              ['scp',f'{host}:{remotedir}/ansys.out',str(jobdir/'ansys.out')]]
    if dry_run:return commands
    for i,cmd in enumerate(commands):
        print(f'[{i+1}/{len(commands)}] {shlex.join(cmd)}',flush=True)
        try:subprocess.run(cmd,check=True,timeout=3600)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f'Remote step {i+1} failed (exit {e.returncode}). '
                               'Verify OpenSSH, Tailscale, directory permissions, ANSYS license, and binary path.') from e
    return finalize(jobdir)


def compare(folder_paths,outdir,allow_mixed=False):
    """User-defined equal-weight rank aggregation; official score weights are unspecified."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    paths=[];records=[]
    for folder in folder_paths:
        p=Path(folder)
        if p.is_dir():p=p/'result.json'
        d=json.loads(p.read_text(encoding='utf-8'))
        m=d['official_metrics']
        paths.append(str(p))
        records.append({'name':p.parent.name,'path':str(p),'mass_g':float(m['mass_g']),
                        'temperature_C':float(m['pizza_temperature_20min_C']),
                        'retention_pct':float(m['pizza_moisture_retention_20min_pct']),
                        'qualified':bool(m['temperature_threshold_pass']),
                        'solver':d.get('solver',d.get('model','unknown'))})
    if not records:raise ValueError('Need at least one result.json to rank')
    profiles=[]
    for rawpath in paths:
        d=json.loads(Path(rawpath).read_text(encoding='utf-8'))
        if 'geometry_mesh_mode' in d:
            config=d.get('numerics',{})
            signature=('ansys',d['geometry_mesh_mode'],config.get('dt_s'),config.get('mesh_mm'))
        else:
            num=d.get('numerics',{})
            geom=d.get('geometry',{})
            signature=('quick','voxel',round(float(num.get('dt_s',-1)),5),
                       round(float(num.get('requested_voxel_mm',np.mean(geom.get('voxel_mm_actual',[0])))),2))
        profiles.append(signature)
    if len({str(p) for p in profiles})>1 and not allow_mixed:
        raise ValueError('Incompatible evaluation profiles (solver, mesh mode or step size). '
                         'Evaluate all competitors using the same settings, or explicitly add --allow-mixed.')
    from scipy.stats import rankdata
    ok=[i for i,r in enumerate(records) if r['qualified']]
    if ok:
        ma=rankdata([records[i]['mass_g'] for i in ok],method='min')
        te=rankdata([-records[i]['temperature_C'] for i in ok],method='min')
        mo=rankdata([records[i]['retention_pct'] for i in ok],method='min')
        for j,i in enumerate(ok):
            records[i]['mass_rank']=int(ma[j]);records[i]['temperature_rank']=int(te[j])
            records[i]['moisture_rank']=int(mo[j]);records[i]['rank_sum']=int(ma[j]+te[j]+mo[j])
        for r in records:
            if not r['qualified']:r['rank_sum']=None
    sorted_records=sorted(records,key=lambda d:(not d['qualified'],float('inf') if d.get('rank_sum') is None else d['rank_sum'],
                                                -d['temperature_C']))
    for k,r in enumerate(sorted_records,1):r['internal_placement']=k if r['qualified'] else None
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    with open(outdir/'ranking.csv','w',newline='') as f:
        keys=['name','qualified','mass_g','temperature_C','retention_pct','mass_rank','temperature_rank',
              'moisture_rank','rank_sum','internal_placement','solver','path']
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(sorted_records)
    (outdir/'ranking.json').write_text(json.dumps({'method':'Equal-weight sum of 3 ordinal ranks (INTERNAL ONLY)',
                                                    'official_weights_known':False,'profiles':profiles,
                                                    'mixed_profiles':len({str(p) for p in profiles})>1,'ranking':sorted_records},
                                                    indent=2,ensure_ascii=False),encoding='utf-8')
    if records:
        fig,ax=plt.subplots(figsize=(7,5));
        for r in records:
            ax.scatter(r['mass_g'],r['temperature_C'],s=45)
            ax.annotate(r['name'],(r['mass_g'],r['temperature_C']),xytext=(4,4),textcoords='offset points')
        ax.axhline(68,ls='--',color='grey');ax.grid(alpha=.25)
        ax.set(xlabel='Box mass (g)',ylabel='Pizza temperature at 20 min (C)',title='Pizza Box Contest')
        fig.tight_layout();fig.savefig(outdir/'ranking.png',dpi=160);plt.close(fig)
    return sorted_records


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub=ap.add_subparsers(dest='action',required=True)
    p=sub.add_parser('prepare',help='STEP + DXF => ANSYS .inp and PowerShell runner')
    p.add_argument('--step',required=True);p.add_argument('--dxf',required=True);p.add_argument('--out',required=True)
    p.add_argument('--mesh-mode',choices=['voxel','step'],default='voxel',help='step requires gmsh; exact STEP-conforming mesh')
    p.add_argument('--mesh-mm',type=float,default=1.8);p.add_argument('--dt',type=float,default=10)
    p.add_argument('--sensors',type=int,default=64)
    p=sub.add_parser('ssh',help='Upload generated job, execute ANSYS on Windows through SSH, retrieve & score')
    p.add_argument('--job',required=True);p.add_argument('--host',required=True)
    p.add_argument('--remote-base',default='pizza_box_contest/jobs');p.add_argument('--nproc',type=int,default=4)
    p.add_argument('--ansys-exe',default='');p.add_argument('--dry-run',action='store_true')
    p=sub.add_parser('finalize',help='Compute moisture/metrics after collecting ANSYS CSV files')
    p.add_argument('--job',required=True)
    p=sub.add_parser('rank',help='Compare result.json files; internal equal-rank sum')
    p.add_argument('results',nargs='+');p.add_argument('--out',required=True)
    p.add_argument('--allow-mixed',action='store_true',help='Permit incompatible fidelity/mesh settings (not recommended)')
    a=ap.parse_args()
    if a.action=='prepare':
        result=prepare(a.step,a.dxf,a.out,a.mesh_mm,a.dt,a.sensors,a.mesh_mode)
        print(json.dumps({'status':'prepared_not_run','mesh':result['mesh'],'mass_g':result['mass_g'],
                          'mesh_mode':result['mesh_mode']},indent=2))
    elif a.action=='ssh':
        result=remote_run(a.job,a.host,a.remote_base,a.nproc,a.ansys_exe,a.dry_run)
        print(json.dumps(result,indent=2,ensure_ascii=False))
    elif a.action=='finalize':print(json.dumps(finalize(a.job),indent=2,ensure_ascii=False))
    elif a.action=='rank':print(json.dumps(compare(a.results,a.out,a.allow_mixed),indent=2,ensure_ascii=False))

if __name__=='__main__':main()
