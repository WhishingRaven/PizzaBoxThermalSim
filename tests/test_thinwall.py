from pathlib import Path
import math
import cadquery as cq
import numpy as np
import pytest
from scipy.sparse.linalg import spsolve
from simulator import DEFAULTS, GeometryError
from thinwall import build_shell, assemble, half_wall_film, solve, sensors, write_shell_apdl
from validate_thinwall import make_rounded

BASE=Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('thickness',[.0015,.003,.006])
def test_slab_resistance_and_storage_are_physical(thickness):
    # A single isothermal triangle between a hot contact and ambient has an exact series solution.
    p=DEFAULTS.copy();p['h_internal_W_m2K']=1e-15
    mesh={'xyz':np.array([[0.,0,0],[1.,0,0],[0,2.,0]]),'triangles':np.array([[0,1,2]]),
          'area':np.array([1.]),'floor':np.array([False]),'contact':np.array([True]),
          'thickness_m':thickness,'free_air_volume_m3':.004}
    sys=assemble(mesh,p);n=3;hot=90.;cold=p['T_environment_C']
    wall=spsolve(sys['H'][:n,:n],sys['source'][:n]-sys['H'][:n,n+1].toarray().ravel()*hot)
    q=sys['h']['exterior']*(wall.mean()-cold)
    exact=(hot-cold)/(1/p['h_contact_W_m2K']+thickness/p['box_k_W_mK']+1/p['h_external_W_m2K'])
    assert q==pytest.approx(exact,rel=1e-12)
    assert sys['capacity'][:n].sum()==pytest.approx(200*1400*thickness)
    assert half_wall_film(10,.003,.08)==pytest.approx(1/(.1+.003/.16))


@pytest.mark.parametrize('rounded',[False,True])
def test_cad_area_mass_contact_and_energy_invariant(tmp_path,rounded):
    step,dxf=make_rounded(tmp_path/'cad') if rounded else (BASE/'examples/reference_box.step',BASE/'examples/reference_net.dxf')
    meshes=[build_shell(step,dxf,h,DEFAULTS) for h in (24.,12.)]
    a,b=meshes
    for key in ('area_m2','thickness_m','volume_m3','mass_g','cavity_volume_m3','free_air_volume_m3'):
        assert a[key]==pytest.approx(b[key],rel=1e-10)
    for mesh in meshes:
        assert mesh['area'][mesh['contact']].sum()==pytest.approx(math.pi*.15**2,rel=1e-10)
        assert abs(mesh['shell_volume_error_pct'])<.01
        sys=assemble(mesh,DEFAULTS)
        assert sys['capacity'][:-2].sum()==pytest.approx(mesh['mass_g']*.001*1400,rel=1e-12)
        assert sys['capacity'][-2]==pytest.approx(1.2*1005*mesh['free_air_volume_m3'])
        assert sys['capacity'][-1]==pytest.approx(.75*3000)
        assert np.max(abs((sys['H']-sys['H'].T).data),initial=0)<1e-12
    p=DEFAULTS.copy();p['duration_s']=20
    result=solve(a,p,2,sensors(a,16))
    assert abs(result['energy_balance_error_J'])<.01


def test_openings_and_overlapping_thickness_are_rejected(tmp_path):
    outer=cq.Workplane('XY').box(326,326,56,centered=(False,False,False))
    inner=cq.Workplane('XY').box(320,320,50,centered=(False,False,False)).translate((3,3,3))
    solid=outer.cut(inner)
    vent=cq.Workplane('XY').center(163,163).circle(10).extrude(10).translate((0,0,50))
    path=tmp_path/'vent.step';cq.exporters.export(solid.cut(vent),str(path))
    with pytest.raises(GeometryError,match='openings/vents'):
        build_shell(path,BASE/'examples/reference_net.dxf',24,DEFAULTS)
    tab=cq.Workplane('XY').box(100,30,3,centered=(False,False,False)).translate((0,0,56))
    cq.exporters.export(solid.union(tab),str(path))
    with pytest.raises(GeometryError,match='Overlapping tabs'):
        build_shell(path,BASE/'examples/reference_net.dxf',24,DEFAULTS)


def test_export_is_shell_with_explicit_thickness_and_interfaces(tmp_path):
    mesh=build_shell(BASE/'examples/reference_box.step',BASE/'examples/reference_net.dxf',24,DEFAULTS)
    meta=write_shell_apdl(mesh,tmp_path,DEFAULTS,10,sensors(mesh,16))
    deck=(tmp_path/'run.inp').read_text()
    assert 'SHELL131' in deck and 'SOLID70' not in deck and 'SOLID87' not in deck
    assert 'SECDATA,0.003,' in deck
    assert meta['wall_thickness_m']==.003
    assert 'heat_flux_raw' in deck and 'CW(' in deck


def test_concave_panel_junction_has_no_duplicate_material(tmp_path):
    from thinwall import _offset
    cavity=cq.Workplane('XY').polyline([(0,0),(320,0),(320,310),(310,310),(310,320),(0,320)]).close().extrude(50).translate((3,3,3)).val()
    solid=_offset(cavity,3).cut(cavity)
    path=tmp_path/'junction.step';cq.exporters.export(solid,str(path))
    mesh=build_shell(path,BASE/'examples/reference_net.dxf',24,DEFAULTS)
    assert mesh['area'].sum()==pytest.approx(mesh['area_m2'],rel=1e-12)
    assert abs(mesh['shell_volume_error_pct'])<.1


def test_ranking_preserves_moisture_direction_and_blocks_thickness_mix(tmp_path):
    import json
    from ansys_bridge import compare
    paths=[]
    for name,retention in [('crisp',90),('wet',95)]:
        path=tmp_path/name/'result.json';path.parent.mkdir()
        obj={'model_version':'thinwall-v1','solver':'Python thin-wall FEM','geometry_mesh_mode':'shell',
             'numerics':{'dt_s':10,'mesh_mm':3},'geometry':{'thickness_m':.003},
             'official_metrics':{'mass_g':170,'pizza_temperature_20min_C':69,'pizza_moisture_retention_20min_pct':retention,'temperature_threshold_pass':True}}
        path.write_text(json.dumps(obj));paths.append(path)
    assert compare(paths,tmp_path/'rank')[0]['name']=='crisp'
    obj['geometry']['thickness_m']=.006;paths[1].write_text(json.dumps(obj))
    with pytest.raises(ValueError,match='Incompatible'):compare(paths,tmp_path/'mixed')


def test_richardson_uses_actual_unequal_refinement_ratios():
    from validate_thinwall import estimate_gci
    h=[.006,.0032,.0015]
    temperature=[70+1000*x*x for x in h]
    result=estimate_gci(temperature,h)
    assert result['observed_order']==pytest.approx(2,abs=1e-8)
    assert result['richardson_temperature_C']==pytest.approx(70,abs=1e-10)
    assert 'GCI_status' in estimate_gci([70,69,69.5],h)


def test_failed_ansys_output_cannot_be_accepted(tmp_path):
    from thinwall import finalize_shell
    (tmp_path/'ansys.out').write_text('*** FATAL ***\nOld CSV files might still exist\n')
    with pytest.raises(RuntimeError,match='log contains errors'):
        finalize_shell(tmp_path,{})
