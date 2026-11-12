# 피자 박스 설계 컨테스트 시뮬레이터

2026년 2학기 Heat Transfer Project #2의 공통 조건을 적용하는 **독립 평가 도구**입니다. 조교의 비공개 FEM 구현이나 실제 점수 산식과 동일하다고 보장하지 않습니다. 빠른 로컬 FEM과 Windows ANSYS Mechanical APDL 배치 솔버를 지원합니다.

## 1. 어떤 파일을 넣나요?

- **접힌 상태**: 밀폐된 내부 공간 하나를 가진 단일 솔리드 `my_box.step` (`.stp`도 가능). CAD 길이 단위는 **mm**. 바닥은 평평해야 하고, 내벽 안에 지름 300 mm, 두께 15 mm의 피자가 들어가야 합니다.
- **전개도**: `my_box.dxf`. 단위 **mm** (DXF 파일에 단위 설정 권장). **CUT** 레이어에 폐곡선 `LWPOLYLINE`으로 전체 외곽을 표시하세요. 뚫린 구멍도 CUT 폐곡선으로 표현할 수 있습니다. **FOLD** 레이어에 접는 선을 LINE으로 표시하세요. 탭과 겹치는 부분을 포함한 **전체 판재 면적**이 질량에 반영됩니다.
- STEP과 DXF가 동일한 판재 설계를 표현하도록 만드세요. 프로그램은 STEP의 부피를 이용한 질량과 DXF의 면적을 이용한 질량의 차이가 12%를 넘으면 거부합니다. 이 검사는 **실제 접힘 가능 여부, 겹침, 누수 여부를 완전히 증명하지는 않습니다.** CAD에서 간섭·전개 가능성을 추가 검사해야 합니다.
- 예시 `examples/reference_box.step`과 `examples/reference_net.dxf`는 수치 시험을 위한 **설명용 형상**입니다. 전개도의 접힘 운동학과 모든 접합 탭을 검증한 실제 제작 가능 설계는 아닙니다.

## 2. 설치 (Mac 또는 Python이 있는 컴퓨터)

Python 3.11~3.13 환경을 권장합니다. 터미널에서:

```bash
cd pizza_box_contest
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

ANSYS용 STEP 형상 추종 2차 사면체 메시 기능도 사용할 경우, 가상환경을 활성화한 상태에서 추가 설치:

```bash
python -m pip install -r requirements-exact-step.txt
```

이후 GUI와 CLI 명령은 `.venv`를 활성화한 터미널에서 실행하세요. 가상환경을 종료하려면 `deactivate`를 실행합니다.

Gmsh 및 CAD 종속 라이브러리를 사용할 수 있어야 합니다. Windows 워크스테이션에 Gmsh/Python을 설치할 필요는 없습니다. ANSYS 입력 파일은 로컬에서 생성합니다.

**GUI 실행:** `python gui.py` 실행 후 STEP, DXF를 선택하고 `Quick FEM run` 또는 `Prepare ANSYS` → `Remote ANSYS via SSH`를 누르면 됩니다. Tkinter가 없는 Python 빌드는 아래 CLI를 사용하세요.

## 3. 빠른 로컬 결과 (STEP → 3D 육면체 FEM)

```bash
python simulator.py \
  --step examples/reference_box.step \
  --dxf examples/reference_net.dxf \
  --out results/reference_fast \
  --voxel-mm 3 --dt 10
```

결과: `result.json`, `temperature.csv/png`, `moisture.csv/png`, `heat_flux.csv/png`, `thermal_field.npz`.

- `--voxel-mm 3`: 빠른 초기 탐색. 3 mm 판재를 경계에서 계단형으로 근사하는 오차가 있습니다.
- `--voxel-mm 2` 또는 `1.5`: 더 고해상도 계산. 모델에 따라 소요 시간이 크게 늘어납니다.
- 형상이 작게 바뀌어 근소하게 성능이 달라지는 경우 **최소 두 공간 해상도**와 `--dt 10`/`--dt 5`를 비교하세요.

예시 기준 상자의 **수치 민감도**를 실제 계산했습니다(시간 간격 모두 10 s):

| 공간 격자 | 20분 후 피자 온도 | 수분 잔존율 |
|---|---:|---:|
| 3.0 mm | 69.30°C | 93.19% |
| 2.5 mm | 68.62°C | 93.25% |

3 → 2.5 mm에서 피자 온도가 **0.67°C 변했습니다.** 아직 메시 수렴이 증명된 것이 아니며, 68°C 합격선에 가까운 설계를 이 수준에서 확정하면 위험합니다. **더 작은 격자 또는 STEP 형상 추종 ANSYS FEM으로 수렴성을 확인해야 합니다.** 이 값은 모두 조교의 공식 정답이 아닌 독립 추정치입니다.

## 4. ANSYS로 더 정교한 FEM (Windows SSH)

### 4-1. ANSYS 입력 파일 생성

**권장: 실제 STEP 경계에 맞추는 사면체 요소 (Gmsh 필요)**

```bash
python ansys_bridge.py prepare \
  --step designs/team_A.step --dxf designs/team_A.dxf \
  --out jobs/team_A --mesh-mode step --mesh-mm 2.0 --dt 10 --sensors 64
```

이 방식은 Gmsh로 STEP을 따라가는 10절점 사면체 체적 메시를 만들고, **ANSYS SOLID87**로 열전도를 계산합니다. 단, 대류 표면은 SURF152 선형 삼각면으로 근사하며 피자 접촉 영역도 메시 크기의 영향을 받습니다. **이 STEP 메시 경로는 제공된 개발 환경에 Gmsh가 없어 실제 ANSYS에서 검증되지 않았습니다.** 먼저 아래 `voxel` 경로로 ANSYS 연동을 확인한 뒤 시험하세요.

**대안: 구조가 복잡하거나 STEP 메시가 실패하면 계단형 육면체 메시**

```bash
python ansys_bridge.py prepare \
  --step designs/team_A.step --dxf designs/team_A.dxf \
  --out jobs/team_A --mesh-mode voxel --mesh-mm 2.0 --dt 10 --sensors 64
```

이 방식은 **ANSYS SOLID70** 요소를 사용합니다. 해석 솔버는 ANSYS이지만 형상의 표면은 여전히 voxel 근사이므로 "정확한 STEP 형상 FEM"이 아닙니다.

두 방식 모두 피자와 공기는 **MASS71 집중 열용량**, 피자-공기 열전달은 **LINK33**, 내벽-공기·피자-바닥·외벽 대류/접촉은 **SURF152** 요소를 이용합니다. 외측 바닥 온도는 22°C 고정입니다. 전체 1200초 비정상 열전달을 풉니다. ANSYS 도구에서 사용할 수 있는 요소 및 명령은 공식 Element Reference를 확인하세요.

생성 파일: `jobs/team_A/run.inp`, `metadata.json`, `run_ansys.ps1`, `wall_patches.csv`.

### 4-2. Windows 워크스테이션 준비

워크스테이션에서 이미 Tailscale + SSH를 사용한다면 기존 설정을 그대로 이용할 수 있습니다. 워크스테이션에는 **Mechanical APDL 실행 파일과 유효한 ANSYS 라이선스**가 있어야 합니다. `ansys261.exe`, `ansys252.exe` 등 버전에 따라 이름이 다릅니다. **Workbench만 설치되어 있고 MAPDL 실행 파일/라이선스가 없으면 이 배치 작업은 실행되지 않습니다.**

Mac에서 먼저 확인:

```bash
ssh WINDOWS_USER@TAILSCALE_IP
```

IP는 기존 Tailscale에서 확인된 예시이므로 현재 실제 IP를 `tailscale status`로 재확인하세요. `WINDOWS_USER`는 워크스테이션의 실제 SSH 로그인 계정으로 바꾸세요. Windows OpenSSH Server가 동작해야 합니다. SSH 키 인증을 권장합니다.

### 4-3. 원격 자동 실행

```bash
python ansys_bridge.py ssh \
  --job jobs/team_A \
  --host WINDOWS_USER@TAILSCALE_IP \
  --nproc 8
```

실행 순서는 **SSH로 원격 작업 폴더 생성 → SCP 업로드 → Windows PowerShell로 ANSYS 배치 실행 → 결과 CSV/로그 다운로드 → 로컬 수분 계산 및 평가 JSON 생성**입니다. 비밀번호는 코드에 저장되지 않습니다.

프로그램이 `C:\Program Files\ANSYS Inc`와 그 아래 `ANSYS Student` 폴더에서 MAPDL 실행 파일을 자동 탐색합니다. 설치 경로가 다른 경우:

```bash
python ansys_bridge.py ssh --job jobs/team_A \
  --host WINDOWS_USER@TAILSCALE_IP --nproc 8 \
  --ansys-exe 'D:\Ansys\v261\ANSYS\bin\winx64\ansys261.exe'
```

예를 들어 ANSYS Student 2026 R1은 `C:\Program Files\ANSYS Inc\ANSYS Student\v261\ANSYS\bin\winx64\ANSYS261.exe`에 설치될 수 있으며 자동 탐색 대상입니다.

전송/실행 명령만 확인하려면 `--dry-run`을 붙이세요. Windows의 SSH/SCP 경로 설정에 따라 경로 표기가 조정되어야 할 수 있습니다. 기본 원격 경로는 **Windows SSH 로그인 사용자의 홈을 기준으로** `pizza_box_contest/jobs/팀이름`입니다.

원격 계산에서 `history.csv` 및 `wall_samples.csv`를 수동으로 받아 둔 경우:

```bash
python ansys_bridge.py finalize --job jobs/team_A
```

### 4-4. 정밀도 검증

동일 설계를 `--mesh-mm 3`, `2`, `1.5` 등으로 다시 계산하여 피자 온도와 수분 잔존율이 안정되는지 확인하세요. 또한 `--dt 10`과 `--dt 5`를 비교하세요. **메시/시간 수렴을 확인하지 않은 단일 ANSYS 결과를 '정확한 정답'으로 취급하면 안 됩니다.** 3D FEM의 수치 오차와 공기/접촉/복사 생략 같은 *모델링 오차*는 별개입니다.

## 5. 팀별 순위 비교

같은 솔버·메시 모드·시간 간격·메시 해상도로 얻은 결과만 비교하세요.

```bash
python ansys_bridge.py rank \
  results/team_A/result.json \
  results/team_B/result.json \
  results/team_C/result.json \
  --out results/scoreboard
```

생성 파일: `ranking.csv`, `ranking.json`, `ranking.png`.

**과제에서 지정한 공식적인 세 지표**를 기준으로 먼저 68°C 기준 미달을 구분하고, 기준을 통과한 팀들의 질량/온도/수분 각각의 순위를 계산한 뒤 **세 순위의 합**으로 *팀 내부용* 순위를 매깁니다. 과제에는 세 지표의 정량적 배점이 명시되어 있지 않으며, **설계 근거와 최적화 방법은 사람이 별도 평가**해야 합니다. 복잡한 CFD, 잠열 연성, 복사나 벽면 재증발 등은 기본 모델에 없습니다.

다른 솔버 또는 해상도를 억지로 섞어 순위 계산하려면 `--allow-mixed`가 필요하지만, 컨테스트에는 권장하지 않습니다.

## 6. 사용한 공통 물성

| 항목 | 값 |
|---|---|
| Pizza | diameter 0.300 m, thickness 0.015 m, mass 0.750 kg, cp 3000 J/(kg K) |
| Board | k 0.08 W/(m K), rho 200 kg/m³, cp 1400 J/(kg K), nominal thickness 3 mm |
| Air | rho 1.2 kg/m³, cp 1005 J/(kg K) |
| Internal h / external h / pizza contact hc | 3 / 10 / 10 W/(m² K) |
| Initial pizza/air/box, ambient/floor | 90/22/22 °C, 22 °C |
| Moisture | initial pizza water 40%, air RH 50%, aw = 1, km 1e-3 m/s, Rv 461.526 J/(kg K) |
| Duration | 1200 s |
| Box mass allowance | 83–249 g; pizza temperature threshold 68°C |

열해석에서는 수분에 의한 잠열 피드백을 무시합니다(과제 지시사항). 피자 표면과 측면은 내부 공기와 대류하고, 피자 바닥은 상자와 접촉 열전달을 합니다. 내부 증발·상자 벽 응축은 과제 식에 따라 순차 계산하며, 수분 재증발은 기본 모델에서 무시합니다. 실제 교실 채점 프로그램의 세부 가정과 다르면 결과 역시 달라질 수 있습니다.

## 7. 팀 운영 추천

1. 모든 팀이 **동일한 STEP, DXF 규칙**에 따라 파일을 제출합니다.
2. `--voxel-mm 3`으로 빠르게 후보를 걸러냅니다.
3. 상위 후보를 **동일한 고해상도 설정**으로 ANSYS에 올립니다.
4. 최종 후보의 공간·시간 수렴 확인 후 순위를 계산합니다.
5. 수치 결과와 함께 **실제로 접을 수 있는지**, **간섭이 없는지**, **피자가 바닥에 닿는지**, **내부가 밀폐되는지**를 사람이 최종 확인합니다.

### 일반적인 오류

- `Expected one sealed air cavity`: STEP에 구멍이 있거나, 뚜껑이 완전히 닫히지 않았거나, 메시가 너무 큽니다.
- `Flat pattern must represent exactly one connected sheet`: DXF에 분리된 판재가 있습니다.
- `Pizza support base varies in height`: 피자 전체 바닥 접촉이 만족되지 않습니다.
- `FEM matrix ...`: 잘못된 형상/요소 또는 솔버가 수렴하지 않습니다.
- `No module named gmsh`: `--mesh-mode step` 대신 `voxel`을 쓰거나 gmsh를 설치하세요.
- `ANSYS_RESULT_READY`가 출력되지 않음: 워크스테이션의 `ansys.out`과 라이선스 로그를 확인하세요.

## 8. 신뢰성 범위

**검증한 것:** 로컬 FEM 기준 상자 실행, ANSYS용 voxel APDL deck 생성, SSH 명령 구성, ANSYS 결과 CSV 후처리(합성 테스트), 팀별 순위 계산. 2026-11-13에 `workstation`으로 SSH/SCP 연결해 ANSYS Student 2026 R1 자동 탐색과 기준 예제 voxel 해석(31,106개 체적 요소, 120개 시간 단계)을 실행했으며, MAPDL이 오류 0건으로 완료되고 결과 CSV와 평가 JSON을 생성했습니다.

**검증하지 못한 것:** 실제 설계별 계산, 정확한 STEP 기반 Gmsh 사면체 생성, 메시/시간 간격 수렴성, 조교 채점 모델과의 결과 일치. 기준 예제 성공이 다른 입력 설계의 성공이나 공식 점수와의 일치를 보장하지 않습니다.
