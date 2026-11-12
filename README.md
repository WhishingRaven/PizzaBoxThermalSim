# 피자 박스 설계 컨테스트 열해석

STEP 접힘 형상과 DXF 전개도를 입력받아 20분 후 온도·수분 잔존율을 계산합니다. Python과 Windows ANSYS MAPDL을 지원하며, 기본 벽 모델은 **두께가 명시된 CAD 중간면 열 셸**입니다. 조교의 비공개 채점 모델과는 별도의 독립 추정 모델입니다.

## 설치

Python 3.11~3.13을 권장합니다. Mac에서:

```bash
cd PizzaBoxTempSolver
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Gmsh/OpenCascade는 로컬에서 CAD **표면** 메시를 생성합니다. Windows 워크스테이션에는 메시 생성용 Python/Gmsh가 필요하지 않습니다. 새 터미널에서는 가상환경을 다시 활성화하고, 종료할 때는 `deactivate`를 사용합니다.

## 입력과 대회 조건

- STEP 길이는 mm입니다. 하나의 유효한 연결 solid, 하나의 밀폐된 내부 공간, 평평한 외부 바닥이 필요합니다.
- 지름 300 mm·두께 15 mm 피자가 중심에 완전히 들어가고, 전체 바닥 면적에서 평평한 받침과 접촉해야 합니다. 피자는 공간적으로 균일한 온도를 갖습니다.
- 재료는 3 mm 판 한 장입니다. DXF의 `CUT` 레이어에 닫힌 `LWPOLYLINE` 외곽/구멍을, `FOLD` 레이어에 직선 접힘 표시를 넣습니다. 단위 선언도 필요합니다. 곡선 CUT은 직선 폴리라인으로 분할합니다.
- 공식 비교에 사용할 질량은 **전개도 전체 면적 × 두께 × 밀도**입니다. 탭과 겹치는 판 면적도 포함하며 83~249 g이어야 합니다. STEP 체적 질량과 DXF 질량 차이는 12% 이내여야 합니다.
- 셸에 실제로 배치되는 열적 재료 질량은 CAD 중간면 면적 × 두께 × 밀도이며 결과에 따로 기록합니다. 질량 차이를 숨기거나 전개도 질량으로 열용량을 보정하지 않습니다.
- **여러 층이 겹친 STEP, 두께가 변하는 접합부는 단층 3 mm로 축소하지 않고 오류로 중단합니다.** 이들은 명시적인 적층 모델이 필요합니다. 필요하면 보존된 체적 STEP 메시 경로로 별도 평가할 수 있으며 셸 결과와 혼합 순위를 매기지 않습니다.
- 열린 통풍구는 밀폐 공기 가정과 다르므로 오류로 중단합니다. 셸을 이어 구멍을 막거나 접합 저항을 임의로 추가하지 않습니다.
- DXF 연결성과 FOLD 표시가 실제 접힘 가능성을 증명하지 않습니다. 치수·절단선·접는 선을 표시하고, 제조 가능한 전개도와 밀폐·간섭 여부는 별도로 입증해야 합니다. 예제 전개도 역시 제작 검증품은 아닙니다.

20분 온도 68°C 미만은 탈락입니다. 내부 순위는 질량은 작을수록, 온도는 높을수록, 수분 잔존율은 **낮을수록** 유리한 세 순위의 합입니다. 조교의 공식 가중치는 알 수 없습니다. “작게 만들면 온도 손해가 없다”, “뚜껑·옆면 두께는 효과가 없다”는 설계 가설이며 모델에 정답으로 입력하지 않습니다.

## Python 실행

```bash
python simulator.py --step examples/reference_box.step \
  --dxf examples/reference_net.dxf --out jobs/reference_python \
  --mesh-mm 3 --wall-thickness-mm 3 --dt 10
```

`--mesh-mm`는 표면 내 요소 크기입니다. `--voxel-mm`도 호환 별칭으로 받지만 기본 셸 경로에서 두께를 정하지 않습니다. `--wall-thickness-mm`는 STEP 실제 두께와 일치해야 합니다. GUI는 `python gui.py`로 실행합니다.

생성 파일은 `result.json`, `temperature.csv/png`, `heat_flux.csv/png`, `moisture.csv/png`, `thermal_field.npz`입니다. 셸의 NPZ 좌표/연결은 `xyz_m`, `triangles`로 저장됩니다. 기존 voxel NPZ와 다른 스키마입니다.

## 벽의 물리 모델

CAD 내부 공간을 실제 두께의 절반만큼 바깥으로 offset해 중간면을 얻습니다. 전체 두께 offset으로 원래 STEP을 복원할 수 있는지 검사하여 얇은 벽과 겹친 층을 구분합니다. 곡면·모서리에서도 판 면적을 한 번만 계산하고, 맞닿는 표면은 공유 절점을 사용합니다. 열린 경계·비다양체 접합은 거부합니다.

셸 온도는 중간면에 하나씩 두고, 표면 내 전도는 `k*t`, 면적당 열용량은 `rho*c*t`로 계산합니다. 두께 방향은 각 면에 다음 직렬 저항을 둡니다.

```text
내부 유체/피자 -- 1/(h_i*A) -- t/(2*k*A) -- 셸 중심
셸 중심       -- t/(2*k*A) -- 1/(h_o*A) -- 외기
```

따라서 전체 벽 저항은 격자와 무관하게 `t/(k*A)`입니다. 외부 바닥은 22°C로 고정된 floor와 접하므로 셸 중심과 floor 사이에 `t/(2*k*A)`만 둡니다. 중간면 전체를 22°C로 고정하여 반쪽 벽 저항을 없애지 않습니다.

곡면의 삼각형 면적은 CAD 면의 정확한 면적으로 정규화합니다. 이 면적 적분 가중치를 전도·열용량·대류 모두에 사용합니다. ANSYS에서는 같은 가중치를 재료 계수와 표면 필름 계수에 반영하며 **물리 두께는 바꾸지 않습니다**. 이는 곡면의 다각형 근사로 질량이 달라지는 것을 방지하는 수치 적분 보정입니다. 셸 부피 `A*t`와 원래 STEP 부피의 차이가 0.5%를 넘으면 얇은 셸 근사가 부적절하다고 판단해 거부합니다.

피자 접촉 원은 CAD 표면을 분할해 정확한 면적으로 유지합니다. 격자 중심이 원 안에 들어가는지로 접촉 면적을 선택하지 않습니다. 내부 공기는 기존과 같이 **정확한 CAD 공간 체적에서 피자 체적을 뺀 부피**와 `rho_air*c_air*V_air` 열용량을 갖는 완전혼합 체적 모델입니다. 공간적인 유동/CFD 격자는 원래부터 없으며, 공기·피자를 셸 재료로 변환하지 않습니다.

시간 적분은 후진 Euler입니다. 벽 두께 내 일시적 온도 분포는 한 개 중간면 열용량과 분석적 반벽 저항으로 축약합니다. 매우 짧은 시간의 두께 방향 온도 이력을 풀어야 하면 다층 셸이 필요합니다. 복사·누설·잠열 냉각은 기존 모델과 같이 포함하지 않습니다. 수분은 열해석 후 내부 벽 표면 온도로 계산합니다.

## ANSYS 원격 실행

```bash
python ansys_bridge.py prepare --step examples/reference_box.step \
  --dxf examples/reference_net.dxf --out jobs/reference_ansys \
  --mesh-mode shell --mesh-mm 6 --wall-thickness-mm 3 --dt 10
python ansys_bridge.py ssh --job jobs/reference_ansys --host workstation --nproc 1
```

SHELL131의 한 층·중간면 단일 온도 옵션, 명시적 3 mm section, 체적 열용량을 가진 MASS71 공기/피자, SURF152 내부 인터페이스, LINK33 피자–공기 전도를 사용합니다. Python과 같은 반벽 저항을 필름 계수에 직렬로 합칩니다. [SHELL131 공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_elem/Hlp_E_SHELL131.html), [SURF152 공식 문서](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_elem/Hlp_E_SURF152.html).

`workstation`은 Mac의 SSH 별칭이며 다른 환경에서는 실제 별칭/계정으로 바꿉니다. Windows에서 `C:\Program Files\ANSYS Inc`와 `ANSYS Student` 아래 버전 폴더를 탐색합니다. 확인된 실행 파일은 `C:\Program Files\ANSYS Inc\ANSYS Student\v261\ANSYS\bin\winx64\ANSYS261.exe`입니다. 필요하면 `--ansys-exe`로 지정합니다. Student 판은 셸 외에도 표면 요소를 합산한 전체 요소 수 제한이 있으므로 ANSYS 기본 표면 크기는 6 mm로 설정했습니다.

순서는 SSH 폴더 생성 → SCP 입력 전송 → MAPDL 실행 → 이력·벽 온도·열유량·실행 시간 회수 → 로컬 수분/결과 처리입니다. `ANSYS_RESULT_READY`와 결과 파일을 확인합니다. `--dry-run`은 전송/실행 명령만 반환합니다. 수동으로 회수했다면 `python ansys_bridge.py finalize --job jobs/reference_ansys`를 사용합니다. 셸 작업의 로컬 `shell_mesh.npz`와 `metadata.json`도 보관해야 합니다.

## 설계 비교 및 검증

```bash
python ansys_bridge.py rank jobs/team_A/result.json jobs/team_B/result.json --out jobs/ranking
python -m pip install -r requirements-dev.txt
python -m pytest -q
python validate_thinwall.py --phase python
python validate_thinwall.py --phase ansys --host workstation
python validate_thinwall.py --phase report
```

새 결과에는 `model_version=thinwall-v1`이 기록됩니다. 순위 계산은 셸/과거 모델, 두께, 표면 크기, 시간 간격, 솔버가 다른 결과의 혼합을 기본적으로 거부합니다. 수치 검증 목적으로 혼합할 때만 `--allow-mixed`를 명시합니다.

[검증 데이터](validation/thinwall/results.csv), [수렴성 계산](validation/thinwall/convergence.json), [수렴 그림](validation/thinwall/convergence.png), [검증 설명](validation/thinwall/README.md)을 함께 제공합니다. 기준 상자와 둥근 모서리 대조군에 대해 24/12/6/3/2.5/2/1.5 mm를 계산하고, 일정 비율 6→3→1.5 mm에서 Richardson 외삽/GCI를 산출합니다. 둥근 STEP은 곡률에 대한 열해석 대조군이며 그 DXF의 실제 제작 가능성을 증명하지 않습니다.

이전 voxel 결과는 재현을 위해 남겨둡니다. `simulator.py --wall-model voxel --voxel-mm ...`와 ANSYS `--mesh-mode voxel`은 과거 모델, `--mesh-mode step`은 과거 체적 사면체 모델입니다. 새 기본 셸 모델의 수렴 결과와 혼동하지 마세요. 예제에 저장된 `reference_result*`도 이전 voxel 결과입니다.
