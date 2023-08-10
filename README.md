# Integrative analysis of metabolomics and flow cytometry data for understanding myelosuppression mechanisms in a bone marrow organ-on-a-chip model
##  Background
The bone marrow niche serves as a specialized microenvironment crucial for regulating adult stem cell renewal and
differentiation, playing a vital role in maintaining physiological homeostasis. However, the proliferative nature of the
bone marrow niche makes it susceptible to the adverse effects of chemotherapy drugs. These drugs target rapidly
dividing cancer cells and inadvertently affect normal, healthy cells in the bone marrow, leading to myelosuppression.
Myelosuppression is a serious and potentially life-threatening side effect characterized by the destruction of rapidly
proliferating progenitor cells responsible for the production of mature blood cells and platelets in the peripheral cir-
culation. Data from two experiments previously conducted on a bone marrow organ-on-chip (BM-OOC) system were
utilised to understand the impact of chemotherapeutic agents on immune cell count dynamics and biogenic amine
concentration. We specifically investigated the effects of three chemotherapeutic agents known to cause myelosup-
pression: the PARP inhibitor niraparib and two tyrosine kinase inhibitors, imatinib and sunitinib.

## Computational Methods
1. Smoothing Splines: The dynamics of cell counts over the experimental period are modeled using a smoothing spline approach. By determining optimal smoothing parameters through general cross-validation, reliable trends are extracted without overfitting.
2. Exploratory Data Analysis (EDA): An exploratory approach is taken to analyze drug-treated samples, confirming dose-dependent changes in biogenic amine concentrations and cell counts.
3. Sparse Canonical Correlation Analysis (SCCA): Integrating datasets, SCCA identifies highly correlated subsets of metabolites and cell types. This methodology ensures robustness through penalty values determined using a permutational approach, backed by statistical significance confirmed via a permutational test.

## Results
1. A splines based time-series modelling approach revealed patterns of BM differentiation and maturation.
2. Exploratory data analysis (EDA) reveals differences in the dynamics of metabolites and cell counts between control and drug conditions
3. Sparse canonical correlation analysis revealed associations between cell types and metabolites in the BM-OOC

