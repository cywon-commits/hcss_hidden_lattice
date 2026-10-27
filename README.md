# hcss_hidden_lattice

2차원 강체 코어 + 계단 숄더(HCSS) 모형에서 나타나는 정사각–삼각(12회, δ = 1.40σ)과 삼각형–80° 마름모(18회, δ = 1.27σ) 타일링을 다룹니다. 다루는 내용은 세 가지입니다.

- **숨은 삼각 격자:** 족 회전 메커니즘과 대각선 덮개 정리
- **Z₃ 구역 구조:** 소용돌이와 페이존 보존량
- **18회 상의 안정화 기구:** 축퇴 압력과 무른 결함의 자유 부피

출발점은 입자 Voronoi 그물 위의 Ising 좌절을 disclination의 Z₂ 사영으로 보는 연구 제안서입니다(`docs/proposal_notes.md`).

관련 저장소:

- [cywon-commits/hardcore_square_shoulder](https://github.com/cywon-commits/hardcore_square_shoulder): 모형 정의, λ* 논문, 예비 실행 1–10
- [cywon-commits/hcss_lambda_star](https://github.com/cywon-commits/hcss_lambda_star): λ* = 2cos15° 타일링 후속. 이 저장소의 전달 행렬 결과(TASK4, `results06.md`, `reference80/`, `tm80/`)는 그쪽 커밋 `2930f75`에 있고, 여기에는 요약만 둡니다.

## 주 결과 (자세한 내용은 `results01.md`)

| # | 결과 | 종류 | 근거 |
|---|---|---|---|
| 1 | 정사각–삼각 타일링의 홀수 배위 비율은 조성만의 함수 n_odd/N = 2/(1 + x/2)이고 페이존 뒤집기에 불변이다(12회: 4√3 − 6 ≈ 0.928). | 해석 | `geometry/check_cover_patch.py` |
| 2 | **족 정리:** 모서리 방향 30°·m의 족 m mod 2 가운데 한 족을 θ만큼 돌리는 변형은 모든 정사각–삼각 타일링(4⁴ 포함)의 정확한 평평한 메커니즘이다. 따름정리로, 정사각형마다 대각선 하나를 골라 모든 입자를 6배위로 만드는 덮개가 항상 존재하고 Q의 성분마다 2개다(평면과 토러스 모두). | 정리 + 수치 | `geometry/check_cover_patch.py`, `check_holonomy.py`, `develop.py` |
| 3 | 평평한 변형 공간의 차원은 (모서리 류 수 − 1)이다. 12회 패치는 2류(메커니즘 1개), 3³.4²는 줄마다 하나씩 더해진다. | 정리 + 수치 | `geometry/shear_correlation.py` |
| 4 | HCSS에서 메커니즘의 유효 퍼텐셜: 압력은 늘 기울임을 선호하고, δ < √2이면 짧은 대각선이 숄더 벽 θ_s = 90° − 2 arcsin(δ/2)를 만든다. 문헌의 12회 준결정(δ = 1.40, U/N = 2.536, ρ* ≈ 1.07)과 정확히 일치한다. | 해석 | `geometry/diagonal_distances.py` |
| 5 | **꼭짓점 산술:** 예각 90° − θ 마름모에서 o ≠ a 꼭짓점은 θ = 10°, 15°, 18°에서만 가능하다. 18회의 80° 마름모는 정사각형 셋이 만나는 꼭짓점이 평평해지는 가장 작은 기울기다. | 해석 | `results01.md` 2절 |
| 6 | **Z₃ 구역:** 80° 타일링에서 모서리 족 (방향/20°) mod 3, 마름모 표지(빠진 족), 소용돌이 전하 q = (a − o)/3, 벽은 a·a·o·o 꼭짓점. MC에서 예측대로 분류되고(322개 중 320개), 벽형 : 내부형 = 1.93(조합 수 2). | 정리 + MC | `mc/zthree.py` |
| 7 | **상승된 넓이 정리:** 12개 타일 류의 ∧²Z⁶ 2-형식은 선형 독립이다. 경계가 류별 타일 수를 고정하므로 Z₃ 구역 조성은 페이존 보존량이고, 타일 5개 이하에는 국소 재배치가 없다. 삼각형 + 80° 마름모만으로는 두 수직 공간의 넓이를 동시에 0으로 만들 수 없다(σ₅: 1.4845, σ₇: 0.7899). | 정리 + 열거 | `geometry/lifted_rank.py`, `moves.py`, `perp_area.py` |
| 8 | **전달 행렬(hcss_lambda_star TASK4):** 확장 엔트로피는 조성 2.2–2.4의 두 족 영역에만 있고, 이는 기울인 정사각–삼각 무작위 타일링과 엔트로피가 5e-14 안에서 같다. 세 족·소용돌이 영역(≲ 1.56)은 1차원 사슬 엔트로피뿐이다. | 전달 행렬 | hcss_lambda_star `results06.md` |
| 9 | **축퇴 압력:** 삼각형–벽 마름모 타일링은 P_c(0) = ε/[σ²(sin α − √3/2)](δ = 1.27에서 8.7)에서 모든 조성이 축퇴한다. | 해석 | `results01.md` 5절 |
| 10 | **18회 상의 안정화:** 18회 상은 P_c 근처에서 유체를 고정 압력으로 식혀도 저절로 생긴다(Z₃ 무질서, 조성 ≈ 1.4). 이를 받치는 것은 무른 결함(형태 자유도 n − 3인 다각형)의 자유 부피 엔트로피다. 교차 압력은 T* = 0.20–0.30에서 10.8 ± 0.2로 평평하며, 자유 부피 모형 Δs ≈ ln(βf̄ℓ_f)과 맞고 상수 엔트로피 모형은 기각된다. | MC + 해석 | `mc/run_npt_cool.py`, `cage.py`, `fl_any.py` |

## 구성

```
engine/hcss_mc.py        NPT/NVT MC 엔진 (출처는 SOURCE.md)
geometry/                타일링 기하: sqtri.py(Stampfli 팽창·면·덮개), develop.py(메커니즘 전개),
                         check_cover_patch.py, check_holonomy.py, check_matching.py, shear_correlation.py,
                         diagonal_distances.py, moves.py(정수 좌표 패치 열거), lifted_rank.py, perp_area.py,
                         make_sd.py(σ상을 10° 기울인 단일 구역 주기 타일링)
mc/                      입자 MC 드라이버와 분석: run_cool/run_slow/run_anneal(NVT), run_npt(압력 계열),
                         run_npt_cool(고정 압력 냉각), run_pc_T.sh(P_c(T) 스캔), fl_any.py(Frenkel–Ladd),
                         cage.py(케이지 엔트로피), analyze.py(타일 분해), zthree.py(Z₃ 분류), famtri.py, xcomp.py
data/configs/            결과표의 배치 스냅숏 (s, box, lam, T[, P]) — results01.md 부록 표 참조
data/measurements/       실행 로그(*.log)와 자유에너지·케이지 측정(*.json)
figures/                 발표 그림 f1–f8과 make_figs.py
docs/proposal_notes.md   연구 제안서(수정판, 1–14절) 사본
results01.md             이번 단계의 결과 보고
```

## 빠른 재현

```bash
pip install -r requirements.txt
cd geometry && python check_cover_patch.py      # 12회 패치: 덮개 좌절 0, 덮개 후 모두 6배위
python lifted_rank.py && python perp_area.py    # 계수 12, 두 수직 공간 조성 모순
python moves.py 4                               # 타일 4개 이하 재배치 0
cd ../mc
python xcomp.py ../data/configs/coolP105_T0.240.npz          # 조성, 결함, Z₃ 자화
python zthree.py ../data/configs/coolPc_T0.240.npz           # 소용돌이, 벽, 표지
python run_npt_cool.py ../data/configs/fluid_start.npz 10.5 out 150000 1 0.40,0.32,0.28,0.24   # 약 10분
```

N ≈ 1150인 실행은 압력·온도마다 하나씩입니다. 수치 오차와 한계는 `results01.md`에 적었습니다.
