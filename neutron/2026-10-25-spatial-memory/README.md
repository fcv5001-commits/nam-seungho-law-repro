# 남승호법칙 중성자 공간 기억응답 공개 재현

저자·프로젝트: 남승호  
공개 기록일: 2026-10-25 Asia/Seoul

저장된 공간·기억응답의 수치 재현, 공개입력 연결 및 선형 안정조건 검증을 공개한다.
오픈데이터 사용을 허용하며 완전 무입력을 요구하지 않는다.
질량·수명 입력의 재현과 입력하지 않은 관측량에 대한 독립 예측은 구분한다.

## 확인된 결과

| 항목 | 결과와 범위 |
| --- | --- |
| 4개 파수에서 각각 12개 응답극 | 행렬과 모드 다항식의 최대 차이 약 9.05e-15 |
| 원래 응답행렬의 극 잔차 | 최대 상대 잔차 약 8.58e-14 |
| 저장된 기준 극 재현 | 현재 실행에서 최대 차이 0 |
| 공개 질량 에너지 연결 | 939.56542194 MeV 재현, 공개 입력 |
| 공개 고유수명 연결 | 878.4 s 재현, 공개 입력 |
| 시간발전 교차검산 | ODE와 행렬지수의 차이 약 1.12e-10 이하 |
| 선형 안정조건 | 주어진 양의 정부호 행렬과 eta,tau>0에서 모든 실수 k에 대해 충족 |

수치·재현 검산 완료를 기록한다. 중성자 전체의 원작용 직접 유도,
스핀·전하·약전이·전체 중력장 또는 독립 실험 검증이 완료됐다고 주장하지 않는다.

## 실행

저장소 루트에서 실행한다. Python 3.10 이상이 필요하다.

```bash
cd neutron/2026-10-25-spatial-memory
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows에서는 .venv\Scripts\activate
python -m pip install -r requirements.txt
python run.py --output reproduced_results.json
python run.py --lapse 1.000001 --output reproduced_lapse_results.json
python verify.py --output verification_reproduced.json
```

`reference_results.json`은 기존 저장 기준이며 수정하지 않았다.
`results.json`, `results_lapse.json`, `verification.json`은 이번 공개 전 실행 기록이다.
환경별 부동소수점 차이는 코드의 허용오차로 검사한다.
외부 계정, 비밀키 또는 비공개 데이터 없이 포함된 입력으로 실행할 수 있다.
공개 원자료를 매번 다운로드하는 구현은 아니며 출처와 채택값을 입력파일에 기록했다.

## 원자료와 입력 계보

`run.py`와 `inputs.json`은 기존 「남승호법칙 중성자 공간 기억응답 수치구현」에서
바이트를 유지해 복사했다. 입력 JSON은 이전 행렬 기록을 정리한 결합 입력이며,
원래 모든 섹터의 보편 작용함수 원본이라고 표시하지 않는다.

- 질량: NIST 2022 CODATA, 939.56542194(48) MeV.
  https://physics.nist.gov/cuu/Constants/Table/allascii.txt
- 수명: PDG 2025, 878.4(5) s. 이는 해당 판본의 UCN 중심 평균이며,
  서로 다른 수명 측정방법의 차이가 모두 해결됐다는 뜻은 아니다.
  https://pdg.lbl.gov/2025/listings/rpp2025-list-n.pdf
- 입력의 원래 수명 출처: https://pdg.lbl.gov/2025/tables/rpp2025-sum-baryons.pdf
- hbar: 기존 입력의 6.582119569e-22 MeV s.
- H, Kt, Kx, V, A, Q 및 eta,tau: 기존 프로젝트 공간·기억응답 기록.

## 응답식과 공간 투영

```text
D_R(w,k) = H + k^2 Kx - w^2 Kt - i eta w Kt/(1-i tau w)
eta = 0.08233435425000193
tau = 0.21000022206351907
z' = v
v' = -Kt^(-1)(H+k^2 Kx)z - eta r
tau r' = v-r
deltaZ7 = Q q,   Q = (V; A V)
```

전체상태의 결합은 행렬 전체로 다룬다. `Q`는 네 모드를 일곱 좌표로 옮기는 행렬이며,
물리적 중성자의 3차원 밀도·전하·스핀 유도와 동일시하지 않는다.
시간과 파수는 무차원이며 그대로 초나 역미터 단위로 읽지 않는다.

## 공개값 연결과 중력

```text
w0 = lowest positive-frequency k=0 response pole
s = hbar/[2 tau_n (-Im w0)]
E(w) = E_n + s (w-Re w0)
H_g = E_n phi,   phi = Phi/c^2
coordinate_lifetime = tau_n/N
```

마지막 두 관계는 약한 정적 전체 에너지 결합과 외부에서 지정한 일정 lapse N의 검산이다.
공간에 따라 변하는 계량, 전체 T_whole의 계량 변분 또는 중력장 해를 계산하지 않았다.
질량·수명은 보정 입력이므로 위 연결의 일치를 독립 예측으로 사용하지 않는다.

## 선형 안정조건

Kt, H, Kx가 대칭 양의 정부호이면 모든 실수 k에 대해 H+k^2 Kx도 양의 정부호다.
일반화 고유값 m>0인 각 모드의 특성식은 다음과 같다. m은 입자 질량이 아니다.

```text
tau s^3 + s^2 + (tau m + eta)s + m = 0
Routh-Hurwitz determinant: (tau m + eta) - tau m = eta > 0
```

주어진 선형식에서 성장모드가 없다는 조건부 결론이다.
비선형 전체상태의 안정성이나 물리적 중성자라는 동일성을 이 결론으로 증명하지 않는다.

## 보존과 판정

단일물통 출력식, 터빈 회전항을 출력 근거로 사용하는 식, 구성입자 질량의 단순합,
판테온 람다항과 동등 보정항을 사용하지 않는다. 중력의 전체상태 결합 원칙을 유지한다.
터빈 이론 1254 W와 실측 보고 400~700 W는 사용자 제공 프로젝트 기록이며,
본 중성자 계산의 입력이나 이번 검증 대상이 아니다.

이 공개는 대조 후 코드·입력·결과를 고정하는 실행 기록이다.
사전 예측 봉인이나 외부기관 인증, 독립 제3자 실행 재현을 뜻하지 않는다.
`SHA256SUMS.txt`는 파일 동일성 확인용이며 물리적 타당성 인증이 아니다.

© 2026 남승호. 저장소의 기존 권리 표기를 유지한다.
