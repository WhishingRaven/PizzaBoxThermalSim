"""Reproduce shell convergence and baseline/curvature comparison; optional live SSH."""
from pathlib import Path
import argparse
import csv
import json
import math
import shutil

import cadquery as cq
import numpy as np

from simulator import DEFAULTS
from thinwall import run_shell, prepare_shell

ROOT=Path(__file__).resolve().parent
WORK=ROOT/'jobs/thinwall_study'
REPORT=ROOT/'validation/thinwall'
SIZES=[24.,12.,6.,math.sqrt(21),3.5,3.,2.5,2.,1.5]


def make_rounded(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    def prism(w,d,r,z,h):
        solid=cq.Workplane('XY').box(w-2*r,d,h,centered=(True,True,False)).union(
            cq.Workplane('XY').box(w,d-2*r,h,centered=(True,True,False)))
        for x in (-w/2+r,w/2-r):
            for y in (-d/2+r,d/2-r):
                solid=solid.union(cq.Workplane('XY').center(x,y).circle(r).extrude(h))
        return solid.translate((0,0,z))
    solid=prism(326,326,60,0,56).cut(prism(320,320,57,3,50))
    step=out/'rounded.step';dxf=out/'rounded.dxf'
    cq.exporters.export(solid,str(step));shutil.copyfile(ROOT/'examples/reference_net.dxf',dxf)
    return step,dxf


def inputs():
    rounded=make_rounded(WORK/'inputs')
    return {'reference':(ROOT/'examples/reference_box.step',ROOT/'examples/reference_net.dxf'),'rounded':rounded}


def label(h):return f'{h:g}'.replace('.','p')


def study(phase,host):
    from ansys_bridge import PS1,remote_run
    cases=inputs()
    if phase in ('python','all'):
        for name,(step,dxf) in cases.items():
            for h in SIZES:
                out=WORK/'python'/f'{name}_{label(h)}'
                r=run_shell(step,dxf,out,h,10)
                print(name,h,r['official_metrics']['pizza_temperature_20min_C'],r['runtime_s'],flush=True)
        # Additional time refinement at fixed geometry/discretization.
        for name,(step,dxf) in cases.items():
            for dt in (5.,2.5):run_shell(step,dxf,WORK/'python'/f'{name}_6_dt{label(dt)}',6.,dt)
    if phase in ('ansys','all'):
        if not host:raise ValueError('--host is required for live ANSYS validation')
        for name,(step,dxf) in cases.items():
            for h in (24.,12.,6.,math.sqrt(21),3.5):
                out=WORK/'ansys'/f'{name}_{label(h)}'
                meta=prepare_shell(step,dxf,out,h,10)
                (out/'run_ansys.ps1').write_text(PS1,encoding='utf-8')
                print('ANSYS',name,h,meta['mesh'],flush=True)
                remote_run(out,host,nproc=1)
    write_report()


def estimate_gci(temperature, spacing):
    """Generalized Richardson estimate using actual mean element length."""
    from scipy.optimize import brentq
    coarse,medium,fine=temperature
    hc,hm,hf=spacing
    a,b=coarse-medium,medium-fine
    if a*b<=0 or abs(b)<1e-12:
        return {'GCI_status':'Nonmonotone/degenerate triplet; extrapolation unavailable'}
    r32,r21=hc/hm,hm/hf
    def equation(order):
        return r21**order*(r32**order-1)/(r21**order-1)-a/b
    try:
        order=brentq(equation,.01,10.)
    except ValueError:
        return {'GCI_status':'Triplet is not in the positive-order asymptotic range'}
    extra=fine+(fine-medium)/(r21**order-1)
    return {'observed_order':order,'richardson_temperature_C':extra,
            'fine_grid_GCI_relative_temperature_rise':1.25*abs((fine-medium)/(fine-22))/(r21**order-1),
            'effective_surface_spacing_m':list(spacing),'effective_refinement_ratios':[r32,r21]}


def write_report():
    REPORT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for solver in ('python','ansys'):
        for path in sorted((WORK/solver).glob('*/result.json')):
            result=json.loads(path.read_text());g=result['geometry'];n=result['numerics'];m=result['official_metrics'];e=result['heat_transferred_J'];f=result['fluxes_20min_W']
            name=path.parent.name.split('_')[0]
            meta_path=path.parent/'metadata.json'
            meta=json.loads(meta_path.read_text()) if meta_path.exists() else {}
            row=dict(design=name,solver=solver,model_version=result['model_version'],ansys_total_elements=meta.get('mesh',{}).get('total_ansys_elements'),preparation_runtime_s=meta.get('preparation_runtime_s'),runtime_scope=result.get('runtime_scope','CAD import + meshing + thermal solve' if solver=='python' else 'MAPDL solver process wall time'),mesh_mm=n['mesh_mm'],dt_s=n['dt_s'],elements=n['elements'],nodes=n['nodes'],
                     area_m2=g['area_m2'],thickness_m=g['thickness_m'],volume_m3=g['volume_m3'],thermal_mass_g=g['mass_g'],sheet_mass_g=m['mass_g'],
                     shell_vs_step_volume_error_pct=g['shell_volume_error_pct'],pizza_C=m['pizza_temperature_20min_C'],retention_pct=m['pizza_moisture_retention_20min_pct'],runtime_s=result['runtime_s'],
                     cardboard_heat_loss_J=e['cardboard_to_environment'],ambient_heat_loss_J=e['box_to_ambient'],floor_heat_loss_J=e['box_to_floor'],pizza_to_air_J=e['pizza_to_air'],pizza_to_base_J=e['pizza_to_base'],air_to_box_J=e['air_to_box'],
                     energy_balance_error_J=result['energy_balance_error_J'])
            row.update({k+'_W':v for k,v in f.items()});rows.append(row)
    if not rows:return
    with (REPORT/'results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    convergence={}
    for design in ('reference','rounded'):
        for solver in ('python','ansys'):
            rs=sorted([r for r in rows if r['design']==design and r['solver']==solver and r['dt_s']==10],key=lambda r:-r['mesh_mm'])
            if not rs:continue
            baseline=rs[0]
            if max(abs(r['volume_m3']/baseline['volume_m3']-1) for r in rs)>1e-9:raise RuntimeError('Physical volume changed with mesh')
            changes=[{'coarse_mm':a['mesh_mm'],'fine_mm':b['mesh_mm'],'delta_C':b['pizza_C']-a['pizza_C'],
                      'relative_delta_temperature_rise':abs(b['pizza_C']-a['pizza_C'])/(b['pizza_C']-22),
                      'delta_cardboard_heat_loss_J':b['cardboard_heat_loss_J']-a['cardboard_heat_loss_J'],
                      'relative_delta_cardboard_heat_loss':abs(b['cardboard_heat_loss_J']-a['cardboard_heat_loss_J'])/abs(b['cardboard_heat_loss_J'])} for a,b in zip(rs,rs[1:])]
            entry={'successive_refinements':changes}
            # GCI uses mean element spacing; the requested sizes also have a constant ratio.
            grids=(6.,3.,1.5) if solver=='python' else (6.,math.sqrt(21),3.5)
            vals={r['mesh_mm']:r['pizza_C'] for r in rs}
            if all(h in vals for h in grids):
                selected=[next(r for r in rs if r['mesh_mm']==h) for h in grids]
                spacing=[math.sqrt(r['area_m2']/r['elements']) for r in selected]
                entry.update(estimate_gci([vals[h] for h in grids],spacing))
                entry['gci_grids_mm']=grids
            convergence[design+'_'+solver]=entry
    pairs=[]
    for solver in ('python','ansys'):
        base={(r['mesh_mm'],r['dt_s']):r for r in rows if r['design']=='reference' and r['solver']==solver}
        for r in rows:
            if r['design']=='rounded' and r['solver']==solver and (r['mesh_mm'],r['dt_s']) in base:
                b=base[(r['mesh_mm'],r['dt_s'])]
                pairs.append({'solver':solver,'mesh_mm':r['mesh_mm'],'dt_s':r['dt_s'],'rounded_minus_reference_C':r['pizza_C']-b['pizza_C'],'rounded_minus_reference_retention_pp':r['retention_pct']-b['retention_pct']})
    for name in ('reference','rounded'):
        time_rows=sorted([r for r in rows if r['design']==name and r['solver']=='python' and r['mesh_mm']==6.],key=lambda r:-r['dt_s'])
        convergence[name+'_time_refinement']=[{'dt_s':r['dt_s'],'pizza_C':r['pizza_C']} for r in time_rows]
    (REPORT/'convergence.json').write_text(json.dumps({'convergence':convergence,'design_differences':pairs},indent=2))
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for name in ('reference','rounded'):
        for solver in ('python','ansys'):
            rs=sorted([r for r in rows if r['design']==name and r['solver']==solver and r['dt_s']==10],key=lambda r:r['mesh_mm'])
            axes[0].plot([r['mesh_mm'] for r in rs],[r['pizza_C'] for r in rs],'.-',label=name+' '+solver)
        ps=sorted([r for r in pairs if r['solver']=='python' and r['dt_s']==10],key=lambda r:r['mesh_mm'])
    axes[1].plot([r['mesh_mm'] for r in ps],[r['rounded_minus_reference_C'] for r in ps],'.-')
    axes[0].set(xlabel='Surface size (mm)',ylabel='Pizza at 20 min (C)');axes[0].legend(fontsize=8)
    axes[1].set(xlabel='Surface size (mm)',ylabel='Rounded - reference (C)')
    for ax in axes:ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(REPORT/'convergence.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=['all','python','ansys','report'],default='python')
    parser.add_argument('--host')
    args=parser.parse_args()
    if args.phase=='report':write_report()
    else:study(args.phase,args.host)
