# Λ항 제거 고정 파일

저자: **남승호**  
출처: **남승호법칙(Nam Seung-ho Law)**

이 공개 재현판에서는 우주상수 Λ항을 사용하지 않습니다.

```text
E_base^2(z) = 1 + Omega_m * ((1+z)^3 - 1)
E_NS^2(z)   = exp(alpha * ((1+z)^nu - 1))
Lambda      = 0  (별도 항 없음)
```

- `fit_pantheon_ns.py`에는 Λ 매개변수나 Λ 가산항이 없습니다.
- `NO_LAMBDA_MODEL.json`이 기계판독 가능한 고정 명세입니다.
- 현 결과는 대각오차 표본내 비교이므로 `OPEN / UNSEALED`입니다.
- 전체 공분산과 미사용 자료 검증 전에는 독립검증 완료로 판정하지 않습니다.

© 2026 남승호. All rights reserved.
