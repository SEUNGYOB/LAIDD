# Proteogenomic Pipeline (MATLAB) — Steps (0–8)

재현 가능한 유전체–전사체–단백체–인산화단백체 통합 분석을 위한 MATLAB 파이프라인입니다.  
본 저장소는 0–8단계 스크립트와 보조 함수, 그리고 예시 데이터 스키마를 포함합니다.

원본 코드는 https://github.com/doyoungh/Hyeon_et_al_PDAC_Nature_Cancer 에서 발췌했습니다. 

---

## 1. 개요 (원리)

이 파이프라인은 mRNA–단백질–인산화단백질 데이터를 공통 스케일로 정렬·정규화한 뒤,
변이–표현형 연관성과 아형(Subtype)을 단계적으로 추정하고, 최종적으로 외부 코호트에 mRNA 시그니처를 적용하는 흐름으로 설계되어 있습니다.

### 핵심 아이디어
- 동일 환자 집단에서 얻은 여러 오믹스(mRNA / Protein / Phospho)를 동일 분포로 맞춰(quantile/median 기반) 비교 가능하게 만듭니다.
- TMT(proteomics/phospho) 실험 디자인(세트별 레퍼런스 채널)에 맞춰 within-set FC(접기값)를 계산하고, 세트 간 배치를 보정해 세트 통합(merge)합니다.
- 검출률 기준(예: 100%/50%/20%)으로 분석 목적(클러스터링/시그니처/연관성)에 맞는 피처만 선별하여 통계적 검정·클러스터링의 안정성을 높입니다.
- 변이 ↔ (phospho/protein) 차이는 순위합 검정(비모수) + FDR 보정으로, mRNA–Protein 상관은 피어슨/스피어만 + FDR 보정으로 평가합니다.
- 아형은 aoNMF 기반 2단계 접근(고세포밀도 샘플 선별 → 전체 재학습)으로 보다 안정적이고 해석 가능한 시그니처를 도출합니다.
- 마지막으로 mRNA 시그니처를 외부 코호트에 적용할 때, 시그니처–샘플 간 상관(ρ)과 경험적 p-value로 소속 아형을 판정합니다.

### 단계별 원리 요약
1) 정규화(0단계)
   - mRNA: log2(x+1) → quantile normalization → 중앙값(center) 기준 FC
   - Protein/Phospho(TMT):
     - 스캔 단위 PIP ≥ 70 등 품질 필터 → log2 변환
     - 동일 세트 내 레퍼런스 채널 대비 FC 산출 → 세트 내 quantile 정규화
     - Stack/세트 간 정규화로 배치 효과 완화 → 모든 세트를 채널 마스크(ind_used)로 병합
   - 피처 선별: 분석 목적별로 검출률 기준을 다르게 사용
     - 100%: 환자 클러스터링(결측 최소화)
     - 50%: 시그니처 도출(표현 안정성)
     - 20%: 변이–표현 연관(커버리지 우선)

2) 변이–인산화 연관(1단계)
   - 변이 유무(0/1)에 따른 phosphopeptide FC 분포 차이를 rank-sum test로 평가
   - 효과크기(Δmedian) + p + FDR(BH) 함께 저장해 해석 가능성 확보

3) 변이–단백질/인산화 연관(2단계)
   - 검출률 20% 기준의 protein/phospho 집합에서 변이와의 표현 차이 탐색
   - 단백체/인산화단백체 모두에 대해 양/음 연관 리스트를 도출

4) mRNA–Protein 상관(3단계)
   - 동일 유전자의 mRNA FC와 protein FC 간 상관계수(피어슨/스피어만) 계산
   - BH-FDR로 다중검정 보정 → 유의 유전자 리스트 산출
   - 샘플 매칭은 교집합 기준으로 페어링

5–7) 아형(Subtype) 규명(4–6단계)
   - aoNMF(결측 허용 변형 포함) 으로 mRNA → protein → phospho 순서로 아형 추정
   - 2단계 전략: (i) 고세포밀도 샘플에서 초기 구조 학습 → (ii) 전체 샘플로 고정 재학습(fix)
   - 각 오믹스에서 클러스터 멤버십 + 시그니처 유전자/단백/피처 도출
   - 이후 통합(7단계)에서 세 오믹스의 멤버십을 결합(표결/합의/유사도 융합 중 설계안)하여 최종 아형 정의

8) 외부 코호트 적용(8단계)
   - 내부에서 도출한 mRNA 시그니처를 외부 발현 행렬에 투사
   - 시그니처–샘플 간 상관 ρ 계산 + 경험적 p-value 로 통계적 유의성 평가
   - 가장 일치하는 시그니처 기준으로 외부 샘플의 클러스터 소속을 판정

### 데이터 흐름(요약)
원시데이터 → (정규화/FC) → 피처 선별(검출률) → 
[분기A] 변이–표현 연관(1–2) → 결과표
[분기B] mRNA–Protein 상관(3) → 유의 유전자
[분기C] 아형 규명(4–6) → 통합(7) → mRNA 시그니처 외부 적용(8)

### 설계 선택의 이유
- Quantile/FC 설계: 다중 배치(TMT set/채널) 간 분포를 맞추고, 레퍼런스 대비 상대량으로 기술적 변동을 억제
- 검출률 별도 기준: 목적별 bias–variance 절충(클러스터링은 결측 최소화, 연관 분석은 커버리지 확대)
- 비모수 검정 + FDR: 분포 가정 위배 가능성을 줄이고 다중 비교를 통제
- aoNMF 2-단계: 군집 안정도와 해석 가능성(양수 부하) 확보, 세포성·배치 영향 완화

### 전제/의존성(짧게)
- MATLAB(Statistics Toolbox 필요), quantilenorm2(커스텀) 경로 추가
- 입력 변수(FPKM/IDs/raw/pids/eids/ind_used/PATH_scans/QUANT_PATH) 사전 로드
- 재현성: rng(2025) 고정, 단계별 .mat/.csv 출력 유지

---

## 2. 환경 및 의존성

- MATLAB R2020a 이상 권장  
- 필요 Toolboxes:
  - Statistics and Machine Learning Toolbox  
  - Bioinformatics Toolbox (선택)  
- OS: Windows / macOS / Linux  
- 권장: random seed 고정

---

## 3. 데이터 구조

### 입력 데이터 (예시)
- **FPKM** (#genes × #patients): 유전자 발현량  
- **raw / pids / eids**: 각 TMT set의 peptide–protein–gene 정보  
- **ind_used**: 사용된 TMT 채널 인덱스  
- **mutation_index / ind_mut**: 변이 정보 (patients × genes)  
- **cellularity_rank**: 환자 세포밀도 순위  

### 출력 데이터 (예시)
- **fc_qnorm**, **fc_merged_qnorm**, **fc_merged_qnorm2**: 정규화된 fold-change 행렬  
- **IDs_exp**, **proteins_union_***, **peptides_union_***: 검출된 feature 리스트  
- **clus_***, **sig_***: 클러스터 결과 및 시그니처 세트  
- **int_clus**: 통합 클러스터 결과  

---

## 4. 실행 순서

각 단계는 독립적으로 실행 가능하며, 0단계부터 순차적으로 수행하는 것을 권장합니다.  
아래 예시는 모든 `.m` 파일이 `code/` 폴더 안에 있는 경우를 기준으로 합니다.

```matlab
% 0. 데이터 정규화
run('code/00_Data_normalization.m')

% 1. 변이–인산화 연관 분석
run('code/01_Mutation-phosphorylation_correlation.m')

% 2. 변이–단백질 및 인산화단백질 연관 분석
run('code/02_Mutation-protein_and_phosphopeptide_association.m')

% 3. mRNA–단백질 상관 분석
run('code/03_mRNA-protein_correlation.m')

% 4–6. 아형(Subtype) 규명 단계
run('code/04_Subtype_identification_mRNA.m')
run('code/05_Subtype_identification_protein.m')
run('code/06_Subtype_identification_phospho.m')

% 7. 통합 클러스터링
run('code/07_Integrated_clustering.m')

% 8. 외부 코호트에 mRNA 시그니처 적용
run('code/08_Application_of_mRNA_signatures.m')
