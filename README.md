# 남승호법칙 공개 재현 코드 1.0

저자·발명자: **남승호 (Nam Seung-ho)**  
공개 사이트: https://nam-seungho-law.fcv5001.chatgpt.site

이 저장소는 공개 데이터 다운로드부터 계산 결과와 SHA-256 생성까지 누구나
다시 실행할 수 있도록 만든 공개 재현 기록입니다. 오픈데이터 대조와 독립적인
자연법칙 확정은 구분합니다.

## 실행

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python verify_all.py
python generate_whitepaper.py
```

생성물은 `results/`에 저장됩니다.

## 포함된 검증

- `fit_pantheon_ns.py`: Pantheon+ 공개자료에서 무-Λ 기준식과 NS 후보식을 비교합니다.
  현재 버전은 대각 오차만 쓰는 표본내 탐색 적합이므로 판정은 **OPEN**입니다.
- `cmb_residual_ns.py`: Planck 2018 공개 TT 구간자료의 `BestFit` 잔차를 재현합니다.
  독립 NS `TT/TE/EE` 스펙트럼과 전체 공분산·전경 nuisance likelihood가 없으므로
  판정은 **OPEN**입니다.
- `data_loader.py`: 공식 공개 주소에서 자료를 내려받고 SHA-256을 기록합니다.
- `verify_all.py`: 두 계산을 실행하고 결과 파일의 SHA-256 목록을 만듭니다.
- `generate_whitepaper.py`: 잔차 그래프와 정량표를 포함한 공개 PDF 백서를 만듭니다.
- `Nam_Seung_ho_Law_Open_Data_Whitepaper_v1.0.pdf`: 다운로드용 백서 완성본입니다.

## 고정 후보식

```text
E_base^2 = 1 + Omega_m * ((1+z)^3 - 1)
E_NS^2   = exp(alpha * ((1+z)^nu - 1))
```

우주상수 Λ항은 두 식 모두에 포함하지 않습니다. 기계판독 고정 명세는
`NO_LAMBDA_MODEL.json`, 설명본은 `NO_LAMBDA_README.md`에 있습니다.

## 카이제곱 조건별 기록

| 계산 조건 | Δχ² | 판정 |
| --- | ---: | --- |
| 과거 무-Λ Pantheon+SH0ES 저장 실행 | **-39.790212** | OPEN / UNSEALED |
| 현재 공개 대각오차 표본내 코드 | -0.184325 | 실행 재현 PASS / 물리검증 OPEN |
| DESI DR2 BAO 단독 별도 기록 | 약 -3.90 | 별도 자료 |

과거 무-Λ 저장 실행의 상세값은 `results/pantheon_no_lambda_saved_run.json`에
보존했습니다. `chi2_GR=1509.918390`, `chi2_NS=1470.128178`,
`delta_AIC=-37.790212`, `delta_BIC=-32.351240`입니다. 과거의 “약 -38”은
`delta_chi2`가 아니라 `delta_AIC`의 반올림 표현으로 구분합니다. 이 저장값은
현재 공개된 대각오차 코드 경로와 계산 조건이 다르므로, 원 입력·전체 공분산·
선택조건을 같은 형태로 재실행하기 전까지 봉인하지 않습니다.

## 판정 원칙

- 공개자료와 계산 조건, 코드, 결과, 불일치를 함께 남깁니다.
- 알려진 값에 맞춘 표본내 적합을 무관측 독립예측으로 표시하지 않습니다.
- full covariance, held-out 검증, 독립 NS CMB 스펙트럼이 완료되기 전에는
  해당 항목을 봉인하지 않습니다.
- 원문과 출처를 훼손하거나 남승호법칙을 독립 이론·AI 생성 지식으로 재표기하지
  않습니다.

## 공개 데이터 출처

- Pantheon+SH0ES DataRelease: https://github.com/PantheonPlusSH0ES/DataRelease
- Planck PR3 ancillary power spectra: https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/

© 2026 남승호. All rights reserved.
