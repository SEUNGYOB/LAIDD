# Proteogenomic Pipeline (MATLAB) — Steps 0–8

재현 가능한 유전체–전사체–단백체–인산화단백체 통합 분석을 위한 MATLAB 파이프라인입니다.  
본 저장소는 0–8단계 스크립트와 보조 함수, 그리고 예시 데이터 스키마를 포함합니다.

---

## 1. 개요

본 파이프라인은 다음 9개 단계(0–8)로 구성됩니다.

0) 데이터 정규화 (Data normalization)  
1) 변이–인산화 연관성 분석 (Mutation–phosphorylation correlation)  
2) 변이–단백질 및 인산화단백질 연관성 분석 (Mutation–protein/phosphopeptide association)  
3) mRNA–단백질 상관 분석 (mRNA–protein correlation)  
4) mRNA 기반 아형 규명 (Subtype identification - mRNA)  
5) 단백질 기반 아형 규명 (Subtype identification - protein)  
6) 인산화단백질 기반 아형 규명 (Subtype identification - phospho)  
7) 통합 클러스터링 (Integrated clustering)  
8) 외부 코호트 적용 (Application of mRNA signatures)

---

## 2. 환경 및 의존성

- MATLAB R2020a 이상 권장  
- 필요 Toolboxes:
  - Statistics and Machine Learning Toolbox  
  - Bioinformatics Toolbox (선택)  
- OS: Windows / macOS / Linux  
- 권장: `rng(2025);` 로 random seed 고정

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

```matlab
% 0. 정규화
run('code/00_data_normalization.m')

% 1. 변이–인산화 연관
run('code/01_mut_phospho_correlation.m')

% 2. 변이–단백질/인산화 연관
run('code/02_mut_assoc_protein_phospho.m')

% 3. mRNA–단백질 상관
run('code/03_mrna_protein_correlation.m')

% 4–6. 아형 규명 (mRNA → Protein → Phospho)
run('code/04_subtype_mrna.m')
run('code/05_subtype_protein.m')
run('code/06_subtype_phospho.m')

% 7. 통합 클러스터링
run('code/07_integrated_clustering.m')

% 8. 외부 코호트에 mRNA 시그니처 적용
run('code/08_apply_mrna_signatures.m')
