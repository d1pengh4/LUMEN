# 보충 정보 (Supplementary Information)

*다음 논문의 부록: "어려운 표적에서 de novo 설계가 할 수 있는 것과 없는 것: 폐동맥고혈압을
겨냥한 비면역억제성 FKBP12 리간드를 위한 정직하고 재현 가능한 파이프라인." 아래 모든
표는 각 캡션에 명시된 결과 파일에서 직접 옮긴 것이며, 해당 파일은 프로젝트 저장소의
`results/`에 있다(gitignore됨; 출처와 생성 스크립트는 `docs/PROJECT_LOG.md`와
`docs/research-log/`에 기록).*

---

## 표 S1. 도킹 벤치마크 검증 세트 (재도킹 대조군 + drug-like 앵커)

출처: `results/phase1/docking/docking_scores.csv`; §2.2, §4.4. FK506은 천연 재도킹
대조군이고, 6개 앵커는 정제된 ChEMBL 세트 중 가장 강력한 drug-like(MW ≤ 500) FKBP12
리간드다.

| 리간드 | pChEMBL | 도킹 점수 (kcal/mol) |
|---|---|---|
| FK506 (천연, 재도킹) | 높음 | **−10.30** |
| CHEMBL4090599 | 9.35 | −7.52 |
| CHEMBL4072643 | 8.15 | **−7.57** (앵커 6개 중 최고) |
| CHEMBL4075704 | 8.15 | −7.55 |
| CHEMBL4101268 | 8.15 | −6.93 |
| CHEMBL33094 | 9.00 | −6.85 |
| CHEMBL4102121 | 8.15 | −6.73 |

§2.2에서 논의한 순위 역전에 주목: CHEMBL4090599(pChEMBL 9.35, 앵커 중 가장 강력)는
최고점을 얻지 못했고, CHEMBL4072643(pChEMBL 8.15)이 얻었다. 이 6분자 세트가 정식
60분자 enrichment 시험(표 1)의 동기가 됐다.

---

## 표 S2. GNINA 순위화 전체 37개 후보 풀 (경로 B, 배치 1)

출처: `results/phase2b/FINAL_candidates.csv`; §2.3, §2.5, 그림 1, 표 4. 생성 배치 1에서
나온 developable 33개 + TIER-1 4개 DiffSBDD 후보를 GNINA CNNaffinity(exhaustiveness 4
스크리닝 패스)로 순위화. `ad_tanimoto`는 `developable` 출처 행에만 채워져 있고(ML 효능
모델에 쓰인 적용범위 게이트, §2.2), `dock_score`는 `tier1` 행에만 채워져 있다(GNINA
채점 전에 도킹+developability로 독립 플래그됨, §2.3). 빈 칸은 데이터 누락이 아니라 해당
분자가 어떤 상류 필터를 거쳤는지를 반영한다.

| CNNaffinity | 출처 | SMILES | binder_like | QED | SA | logP | MW | AD Tanimoto | 도킹(kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| 6.560 | tier1 | `CCC[C@@H]1CC(=O)N[C@@]2(C)[C@H](C1)[C@@H](Cc1ccc(F)cc1N)CC[C@@H]2O` | True | 0.719 | 4.39 | 3.42 | 362.2 | — | −8.54 |
| 6.336 | developable | `O=C(O)CCCCN(CCO)c1cccc(F)c1F` | True | 0.711 | 2.28 | 2.02 | 273.1 | 0.224 | — |
| 6.336 | tier1 | `CC1([C@H](Cc2cccc(Cl)c2)OCCNC(=O)O)CCC1` | True | 0.754 | 3.10 | 3.73 | 311.1 | — | −7.58 |
| 6.273 | developable | `CCSC1=NCC[C@@H](CN(CF)C(=O)NCCc2ccccc2)N1` | False | 0.742 | 3.45 | 2.64 | 352.2 | 0.338 | — |
| 6.164 | tier1 | `OCCCC[C@H]1CCCC[C@H](CC2c3ccccc32)N1` | False | 0.776 | 3.43 | 3.59 | 273.2 | — | −7.57 |
| 6.155 | developable | `C=C(C)NNS(=O)(=O)c1ccccc1Cc1ncccn1` | False | 0.788 | 2.63 | 1.38 | 304.1 | 0.197 | — |
| 6.067 | developable | `O=C1[C@@H](Cc2ccccc2)C[C@H]2C[C@H]2N1CC1CCOCC1` | False | 0.855 | 3.46 | 2.89 | 299.2 | 0.291 | — |
| 5.941 | developable | `Cc1ccc(Cn2ccnn2)cc1NC=C1CC1` | False | 0.893 | 2.75 | 2.72 | 240.1 | 0.213 | — |
| 5.820 | developable | `CCCCCN[C@@H]1[C@@H](c2ccccc2)CCCC[C@@H]1O` | False | 0.610 | 3.07 | 3.85 | 275.2 | 0.333 | — |
| 5.803 | developable | `NCCCCS[C@H]1[C@H](N2CCCCC2)CCCS1(=O)=O` | False | 0.757 | 3.78 | 1.85 | 320.2 | 0.169 | — |
| 5.761 | developable | `Cc1ccc(CN2CCC[C@@H](CCCC(=O)O)C2)cc1` | False | 0.864 | 2.27 | 3.46 | 275.2 | 0.328 | — |
| 5.647 | tier1 | `C[C@H]1CCCN(c2ccc(F)cc2)[C@@H]1CCCCC(=O)O` | False | 0.804 | 2.92 | 4.08 | 293.2 | — | −7.34 |
| 5.625 | developable | `C[C@]1(c2cccc(O)c2)CCOC[C@H](CC#N)N1` | False | 0.836 | 3.92 | 1.90 | 246.1 | 0.191 | — |
| 5.620 | developable | `CNC(=O)[C@H]1CCO[C@H](C[C@H]2CCCO2)[C@@H]1c1ccc(C#N)cc1` | False | 0.921 | 3.86 | 2.36 | 328.2 | 0.222 | — |
| 5.583 | developable | `C[C@H]1CC[C@@H]([C@@]2(O)CCC[C@H](N)C2)CC1` | False | 0.700 | 3.47 | 2.45 | 211.2 | 0.159 | — |
| 5.580 | developable | `CCCOC[C@@H]1CC[C@H](O)C[C@@H]1NCC(=O)C1CCCC1` | False | 0.675 | 3.71 | 2.29 | 297.2 | 0.222 | — |
| 5.510 | developable | `Nc1cc(F)ccc1[C@@H]1CCCC[C@@H](C=CC(=O)O)CO1` | False | 0.662 | 3.47 | 3.30 | 293.1 | 0.208 | — |
| 5.476 | developable | `C[C@H]1CCC[C@@](Oc2ccccc2)(C(=O)NCC(=O)O)[C@@H]1O` | False | 0.759 | 3.55 | 1.19 | 307.1 | 0.286 | — |
| 5.464 | developable | `O=C(O)CCCN1Cc2cc(F)ccc2[C@H](F)C1` | False | 0.898 | 2.90 | 2.52 | 255.1 | 0.228 | — |
| 5.183 | developable | `C1CCC(C[C@H]2CCCS2)CC1` | False | 0.628 | 3.20 | 3.85 | 184.1 | 0.173 | — |
| 5.166 | developable | `NC(=O)CCO[C@@H]1CC[C@H]2C[C@@H]2CC(=O)C1` | False | 0.776 | 3.59 | 1.03 | 225.1 | 0.206 | — |
| 5.137 | developable | `NC(=O)[C@@H]1CCCCCCc2c(OOF)cccc2O1` | False | 0.683 | 3.69 | 2.62 | 283.1 | 0.219 | — |
| 5.116 | developable | `CC[C@@H]1C[C@@]1(C)c1ccc(F)cc1N` | False | 0.717 | 3.43 | 3.10 | 193.1 | 0.176 | — |
| 5.060 | developable | `Nc1ccc(Cl)c(CC2NOCS2)c1` | False | 0.764 | 3.97 | 2.02 | 230.0 | 0.150 | — |
| 5.043 | developable | `CCCC[C@H]1CC=C[C@H]2CC[C@@H](C(=O)O)C[C@@H]12` | False | 0.752 | 3.86 | 3.87 | 236.2 | 0.212 | — |
| 4.878 | developable | `NCCc1cccc2c1C(=O)SC2=O` | False | 0.792 | 2.82 | 1.22 | 207.0 | 0.137 | — |
| 4.753 | developable | `COC(=O)N1CCC[C@@H]1C1=CCC=CC1=O` | False | 0.676 | 3.42 | 1.67 | 221.1 | 0.274 | — |
| 4.657 | developable | `CCCCCN1OC(=O)CCC12CC2` | False | 0.647 | 3.67 | 2.26 | 197.1 | 0.172 | — |
| 4.597 | developable | `CC[C@@H](CCC(=O)O)N(C(=O)O)[C@@H]1C(=O)[C@H]1C` | False | 0.730 | 3.77 | 1.20 | 243.1 | 0.188 | — |
| 4.459 | developable | `C[C@@H]1SCC[C@@H]1C(=O)O` | False | 0.604 | 3.88 | 1.21 | 146.0 | 0.170 | — |
| 4.307 | developable | `CCC(=O)CC[C@H](C)C(=O)O` | False | 0.659 | 2.73 | 1.47 | 158.1 | 0.191 | — |
| 4.289 | developable | `C[C@@H]1C(=O)COC(F)(F)[C@H]1c1ccccc1` | False | 0.735 | 3.56 | 2.60 | 226.1 | 0.226 | — |
| 4.202 | developable | `FC(F)(F)O[C@@H](c1ccccc1)C1CC1` | False | 0.748 | 2.79 | 3.67 | 216.1 | 0.230 | — |
| 4.027 | developable | `O=C(O)CNc1ccccc1` | False | 0.682 | 1.45 | 1.18 | 151.1 | 0.277 | — |
| 3.690 | developable | `C=C(C)N(CCC)C(=O)O` | False | 0.655 | 2.99 | 1.91 | 143.1 | 0.155 | — |
| 3.434 | developable | `O=C(O)Cc1ccccc1` | False | 0.665 | 1.21 | 1.31 | 136.1 | 0.343 | — |
| 3.094 | developable | `CN1COc2ccc(O)cc21` | False | 0.604 | 2.56 | 1.18 | 151.1 | 0.216 | — |

상위 3개 행만 active급 임계값(≥6.3)을 넘는다; 나머지 34개는 미끼 수준(5.80 ± 0.53, 표
1) 이하로 채점된다 — §2.3에서 언급한, 포켓 조건부 생성의 정직한 다수결 결과다.

---

## 표 S3. 분자동역학 원자료 궤적 요약 (리드 + 실패한 합성가능 후보)

출처: `results/md/g3_md_rmsd.csv`; §2.4, §2.5, 그림 2. 모든 실행: 1 ns, OpenCL 플랫폼
(§4.6, 본 문서 하단 §S3 참조), 저장 프레임 100개.

| 리드 | CNNaffinity | 리간드 중원자 수 | RMSD 평균 (Å) | RMSD 최대 (Å) | RMSD 최종 (Å) | 판정 |
|---|---|---|---|---|---|---|
| L1 | 7.21 | 32 | 1.15 | 1.70 | 1.67 | PASS |
| L2 | 7.04 | 28 | 1.70 | 2.91 | 1.52 | PASS |
| L3 | 6.95 | 32 | 1.05 | 1.60 | 0.59 | PASS |
| S1 (합성가능 백업, §2.5) | 6.34 | 19 | 2.93 | 6.04 | 5.45 | **FAIL** |

L3의 최종 프레임 RMSD(0.59 Å)는 자체 평균(1.05 Å)보다 눈에 띄게 낮다. 즉 1 ns 구간
끝에서 오히려 초기 포즈로 다시 접근하고 있다는 뜻으로, 멀어지기만 하다가 아직 이탈하지
않은 상태가 아니라 진정으로 안정된 포즈임을 뒷받침한다. S1의 최종 RMSD(5.45 Å)는 자체
평균(2.93 Å)을 초과한다 — 즉 시뮬레이션이 끝나는 시점에도 여전히 포켓에서 멀어지고
있다는 뜻으로, 일시적 이탈이 아니라 진행 중인 이탈과 일치한다.

---

## 표 S4. 선택성 핸들 접촉 기하 (재도킹 원자료)

출처: `results/selectivity/leads_handle_contacts.csv`; §2.6, 표 3, 그림 3. 재도킹은
exhaustiveness 8의 smina를 사용했다(§4.8); 접촉 임계값 ≤ 4.5 Å.

| 리드 | 도킹 점수 (kcal/mol) | Trp59 최소 거리 (Å) | Ile90 최소 거리 (Å) | Trp59 접촉 | Ile90 접촉 |
|---|---|---|---|---|---|
| L1 | −7.70 | 3.44 | 5.16 | yes | no |
| L2 | −6.90 | 3.87 | 3.83 | yes | yes |
| L3 | **−9.00** | 3.67 | 3.59 | yes | yes |

L3는 세 리드 중 재도킹 점수가 가장 우수하고(−9.00 kcal/mol) 표 S3의 최고 MD 안정성과
두 선택성 핸들을 모두 갖췄다 — 실험 로드맵(§3)에서 L3를 최우선으로 둔 종합 근거다.

---

## 표 S5. 리드 3종의 구조적/합성가능성 세부 심사

출처: `results/synthesis/leads_synth_feasibility.csv`; §2.5, §4.7. `synthesizable_easy`는
다른 기준과 함께 SA ≤ 3.5를 요구한다(§4.4); 리드 3종 모두 역합성 시도 전 이 엄격한
심사부터 통과하지 못하며, 이는 AiZynthFinder 결과(표 3)를 예고한다.

| 리드 | SA 점수 | 입체중심 | 총 고리 | 방향족 고리 | 최대 고리 크기 | 변형 고리 | 스피로 중심 | "쉬운" 합성가능성 |
|---|---|---|---|---|---|---|---|---|
| L1 | 4.34 | 3 | 4 | 2 | 7 | No | 0 | No |
| L2 | 4.05 | 3 | 3 | 0 | 6 | No | 0 | No |
| L3 | 4.39 | 3 | 4 | 2 | 7 | No | 0 | No |

---

## 표 S6. 세 생성 풀에 걸친 developability 패널

출처: `docs/research-log/2026-07-01-developability-and-verification.md`
(`results/developability/*.csv`); §2.3, §2.8. RDKit 전용 규칙 기반 developability
대리지표(QED, SA, Lipinski/Veber, PAINS — §4.4)를 도킹이나 GNINA 필터링 이전에 적용,
표 4의 전체 게이트 배치와는 다른 초기 파일럿 규모에서.

| 풀 | n | developable(SA≤4.5, QED≥0.5, Lipinski+Veber, PAINS 없음) | SA(중앙값) | QED(중앙값) | PAINS 경보 |
|---|---|---|---|---|---|
| 경로 A, 상위 50 (BRICS) | 50 | 3 | 4.64 | 0.60 | 0 |
| 경로 B, 전체 도킹 (DiffSBDD 파일럿) | 49 | 16 | 4.75 | 0.57 | 0 |
| 경로 B, 도킹-유망 부분집합 | 8 | 4 | 4.88 | 0.72 | 0 |

도킹-유망 DiffSBDD 부분집합조차 절반만 동시에 쉽게 합성가능하다 — §2.5·그림 1에서
전체 규모로 정식화한 친화도↔합성가능성 긴장을 더 이른, 더 작은 규모에서 미리 확인해준다.

---

## 표 S7. 확장 재현성 노트 (본문 Supplementary Methods §S1–S4도 참조)

| 이슈 | 증상 | 해결 |
|---|---|---|
| gnina 바이너리 크기 | 설정이 멈춘 것처럼 보임 | 1.28 GB 다운로드였을 뿐, 멈춘 게 아니었음; 진행 로깅 추가 |
| `obabel --gen3d` | 특정 SMILES에서 행(hang) | RDKit ETKDG 임베딩으로 전환 |
| 긴 동기식 Colab 실행 | WebSocket 연결 끊김으로 중단 | detached 실행(`nohup`), 폴링, 완료 후 다운로드 |
| 출력 파일명 재사용 | 실행 간 로그가 조용히 손상됨 | 실행/세션마다 고유 파일명 |
| arm64 Colab CLI 래퍼 | x86_64 Python.framework를 강제해 `cryptography` wheel이 dlopen 실패 | 네이티브 uv 설치 `colab` 바이너리로 전환 |
| 실행 중 VM 회수 | ~40분 동기식 `colab exec`가 VM 전체를 잃음(L2 첫 MD 시도) | 클라이언트 연결을 유지하는 harness 추적 백그라운드 실행 선호 |

---

*보충 정보 끝. 위에서 참조한 모든 원본 CSV/로그는 프로젝트 저장소의 `results/`와
`docs/research-log/` 디렉토리에 있다.*
