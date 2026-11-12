from pathlib import Path
import csv
import json
import shutil

import numpy as np
import pytest

from simulator import extract_flat, inspect_step, GeometryError
from ansys_bridge import compare, finalize

BASE = Path(__file__).resolve().parents[1]


def test_reference_dxf_mass_and_geometry():
    flat=extract_flat(BASE/'examples/reference_net.dxf')
    assert flat['fold_count']==5
    assert flat['connected_sheet']
    assert flat['mass_g']==pytest.approx(171.3456,abs=.01)
    shape,bbox,step_mass=inspect_step(BASE/'examples/reference_box.step')
    assert len(shape.Shells())==2
    assert step_mass==pytest.approx(166.2912,abs=.01)
    assert bbox[2,1]-bbox[2,0]==pytest.approx(.056)


def test_reference_quick_thermal_result_exists_and_is_finite():
    path=BASE/'examples/reference_result/result.json'
    obj=json.loads(path.read_text(encoding='utf-8'))
    m=obj['official_metrics']
    assert m['pizza_temperature_20min_C']==pytest.approx(69.299410,abs=.0005)
    assert m['pizza_moisture_retention_20min_pct']==pytest.approx(93.1916,abs=.002)
    assert m['temperature_threshold_pass']
    assert obj['numerics']['water_balance_error_g']<1e-5


def test_ansys_deck_contains_physics_and_sensors():
    folder=BASE/'examples/reference_ansys'
    text=(folder/'run.inp').read_text(encoding='ascii')
    assert 'ET,1,SOLID70' in text
    assert 'ET,2,SURF152' in text
    assert 'ET,4,MASS71' in text
    assert 'ET,5,LINK33' in text
    assert 'NSUBST,120,120,120' in text
    assert 'D,' in text and 'TEMP,22.0' in text
    assert 'SN1(1)=' in text
    assert '*DO,j,1,120' in text
    assert 'SFE,' in text


def test_ansys_postprocess_with_mock_sensor_data(tmp_path):
    # This exercises postprocessing only. DO NOT mistake it for ANSYS validation.
    folder=BASE/'examples/reference_ansys'
    shutil.copy(folder/'metadata.json',tmp_path/'metadata.json')
    hist=np.genfromtxt(BASE/'examples/reference_result/temperature.csv',delimiter=',',skip_header=1)[1:]
    np.savetxt(tmp_path/'history.csv',hist,delimiter=',',
               header='time_s,pizza_C,air_C',comments='')
    count=len(json.loads((tmp_path/'metadata.json').read_text())['sensors'])
    with (tmp_path/'wall_samples.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','sensor_id','wall_C'])
        for t,p,a in hist:
            for sensor_id in range(1,count+1):
                w.writerow([t,sensor_id,22+(a-22)*.3])
    d=finalize(tmp_path)
    assert d['official_metrics']['pizza_temperature_20min_C']==pytest.approx(69.29941)
    assert 80 < d['official_metrics']['pizza_moisture_retention_20min_pct'] <100
    assert (tmp_path/'moisture.csv').exists()


def test_rank_rejects_incompatible_fidelity(tmp_path):
    q=BASE/'examples/reference_result/result.json'
    a=BASE/'examples/reference_ansys/metadata.json'
    src=json.loads(q.read_text())
    x=tmp_path/'a'/'result.json';x.parent.mkdir();x.write_text(json.dumps(src))
    alt=json.loads(q.read_text());alt['numerics']['requested_voxel_mm']=2.
    y=tmp_path/'b'/'result.json';y.parent.mkdir();y.write_text(json.dumps(alt))
    with pytest.raises(ValueError,match='Incompatible'):
        compare([x,y],tmp_path/'ranking')
    got=compare([x,x],tmp_path/'ranking')
    assert len(got)==2
    assert (tmp_path/'ranking'/'ranking.csv').exists()
