KIF1A DATA ANALYSIS - CODE SUMMARY
=====================================
Source: Code.R
Generated: 2026-04-07
=====================================


1. OVERALL PURPOSE
------------------
A comprehensive multi-omics analysis of KIF1A (Kinesin Family Member 1A) gene variants and their
clinical/functional phenotypes. The pipeline integrates variant-level molecular characterization
(motor protein dynamics), computational pathogenicity predictions, clinical outcomes, and machine
learning models to predict phenotypic outcomes (VABS scores, seizure susceptibility, trajectory
clusters) in patients with KIF1A-associated neurological disorder (KAND).


2. INPUT DATA SOURCES
---------------------
- "Deidentified KIF1A Data November 2025.xlsx"
    - Sheet: Cohort Overview     -> patient demographics, mutation type, inheritance, ACMG classification
    - Sheet: Vinelands           -> VABS adaptive behavior scores (longitudinal)
    - Sheet: Yufeng variables    -> 55+ clinical phenotypes (seizure, imaging, developmental regression,
                                   motor milestones, neurological features)
- "KIF1A mutants parameters (2026-04-02) v2.xlsx"
    -> Motor protein dynamics: Movement, Diffusion, Velocity, Run_length, Dwell_time, biochemical properties
- Pathogenicity scores embedded in variant data: VEP, MisFit (v1.5), ESM, CADD, REVEL, gMVP, AlphaMissense
- KIF1A_anno.xlsx  -> variable metadata/annotations


3. MAIN ANALYSIS STEPS
-----------------------
Step 1 - Data Loading & Preprocessing
  - Load variant, clinical, and molecular data from Excel files
  - Standardize values (Yes/No/Unknown/Abnormal/Normal)
  - Handle missing values (floor values for motor scores)
  - Correct clinical records; derive seizure status, EEG results, movement abnormalities

Step 2 - Variable Definition
  - Create lifetime variables (aggregate across visits)
  - Calculate VABS endpoints: final scores, change from baseline (Diff), rate of change (Rate)
  - Create derived variables: Age_VABS_1st, EEG_or_seizure, Dwell_time_group
  - Build endpoint dataset with first and last measurements

Step 3 - Trajectory Analysis
  - VABS trend lines stratified by demographic, genetic, and molecular variables
  - Patient clustering using k-means on mean VABS and slope features
  - Group-based trajectory modeling (GBTM) with 1-4 latent classes (lcmm package)
  - Log-transformed age to account for earlier disease stages

Step 4 - Clustering Analysis
  - Hierarchical clustering of scaled molecular features (ward.D2)
  - PCA of molecular data alone and combined with clinical variables
  - Molecular cluster annotation with VABS phenotypes via boxplots
  - Pairwise correlation heatmaps (pathogenicity vs. molecular features)

Step 5 - Statistical Association Testing
  - Scatter plots: molecular/pathogenicity features vs. VABS outcomes
  - Violin plots comparing features across seizure/EEG status groups
  - Histograms of pathogenicity scores and VABS by age groups
  - Pearson correlation analysis across all variable types

Step 6 - Predictive Modeling (Machine Learning)
  - Leave-One-Variant-Out (LOVO) cross-validation for robust performance
  - Model architectures: Random Forest, Elastic Net (glmnet), XGBoost
  - Ablation studies: Full | Remove Age | Remove Molecular | Remove Both
  - Four prediction targets:
      a) VABS final score (regression)
      b) GBTM trajectory class (classification)
      c) VABS change from baseline (regression)
      d) VABS rate of change (regression)

Step 7 - Clinical Association Testing
  - Logistic regression: molecular phenotype vs. EEG/seizure lifetime occurrence
  - Descriptive statistics by GBTM clusters
  - Linear regression of multiple predictors vs. VABS (age-stratified)
  - Seizure resistance analysis


4. KEY STATISTICAL METHODS & R PACKAGES
-----------------------------------------
Data I/O & Wrangling:
  openxlsx, dplyr, data.table

Bioinformatics:
  Biostrings, seqinr, rentrez (NCBI GenBank codon analysis)

Clustering & Dimensionality Reduction:
  prcomp() + factoextra   - PCA
  cluster (silhouette, clusGap) - clustering validation
  kmeans()                - k-means (nstart=25, iter.max=100)
  hclust()                - hierarchical clustering (ward.D2)
  lcmm                    - group-based trajectory modeling (GBTM)

Machine Learning:
  caret                   - model training framework
  glmnet                  - Elastic Net
  randomForest            - Random Forest (ntree=500)
  xgboost                 - XGBoost gradient boosting
  pROC                    - ROC curves and AUC
  doParallel              - parallel computing

Statistics:
  lm(), glm(), cor.test() (Pearson)
  lme4::lmer()            - linear mixed effects
  splines::ns()           - natural splines

Visualization:
  ggplot2, ggrepel, ggpmisc, ggpubr
  pheatmap, ComplexHeatmap
  ggcorrplot, ggdendro
  cowplot, patchwork, plotly


5. KEY PARAMETERS & THRESHOLDS
--------------------------------
VABS:
  - Scales: ABC (Adaptive Behavior Composite), Communication, Daily Living Skills, Motor, Socialization
  - Exclusion: VABS scores > 100 removed from trajectory modeling

Age cutoffs: 7 and 10 years (for stratified analyses)

PCA Clustering:
  - k_max = 8; optimal k selected by Gap statistic and Silhouette

GBTM:
  - ng_max = 4 (test 1-4 latent trajectory groups)
  - Model selection by entropy

Cross-Validation:
  - LOVO-CV: all samples sharing a variant are withheld together
  - Global tuning: 5-fold repeated CV, 3 repeats

Hyperparameter Grids:
  - Elastic Net: alpha 0-1 (step 0.2), lambda 10^-4 to 10^1
  - Random Forest: mtry varies by feature count
  - XGBoost: nrounds 100-300, max_depth 3-7, eta 0.01-0.1

Other:
  - Correlation significance: p < 0.05
  - Residual annotation cutoffs: |residual| <= 10 or 20
  - Recurrent variant threshold: frequency >= 5


6. MAIN OUTPUTS
---------------
Trajectory Plots:
  - VABS ABC trends colored by Sex, Death, Variants, Movement, Diffusion
  - Trends by pathogenicity scores (REVEL, ESM, AlphaMissense, MisFit)
  - Trends by molecular features (Velocity, Run_length, Dwell_time, conservation, AA properties)
  - Separate trajectories by amino acid position
  - Highly recurrent variants (frequency >= 5)

Clustering Results:
  - Patient clustering (k=2, k=3) by VABS trajectory
  - GBTM trajectory groups (ng=2-4)
  - Molecular feature clustering dendrogram
  - PCA cluster assignments (molecular and molecular+clinical)
  - Heatmaps with VABS boxplot annotations

Association Analysis:
  - Pairwise correlation heatmaps (pathogenicity vs. molecular)
  - 30+ scatter plot combinations (features vs. VABS/pathogenicity)
  - Scatter plots grouped by molecular cluster, movement, and diffusion status
  - Violin plots: seizure/EEG status vs. features
  - Run_length vs. Dwell_time plots with marginal density

Predictive Models:
  - LOVO predictions: actual vs. predicted plots for 4 outcomes
  - Feature importance rankings (top 20 features)
  - Confusion matrices (classification tasks)
  - ROC curves (seizure prediction)
  - Residual plots and per-variant prediction summaries
  - Model ablation comparison tables

Clinical Analysis:
  - Descriptive statistics by GBTM clusters
  - Linear and logistic regression summary tables
  - Statistical report PDFs


7. BIOLOGICAL ENDPOINTS ANALYZED
----------------------------------
Primary Clinical Endpoints:
  - VABS ABC Score (final)             - developmental level
  - VABS Change (Delta)                - longitudinal decline from baseline
  - VABS Rate of Change                - slope of decline per year
  - GBTM Trajectory Classes (2-4)      - developmental subgroups

Secondary Clinical Phenotypes:
  - Seizure/EEG status (lifetime and current)
  - Movement abnormalities (ataxia, dystonia, neuropathy)
  - Motor milestones (sitting, walking, speech onset)
  - Developmental regression (hand skills, language, motor, daily living)
  - Neurological features (optic nerve, hypertonia, hypotonia, microcephaly)

Molecular/Functional Parameters:
  - Movement (binary): motor protein capable of movement
  - Diffusion (binary): diffusion competence
  - Velocity: maximum velocity (um/s)
  - Run_length: distance per run (um)
  - Dwell_time: pause duration (seconds)
  - Amino acid properties: hydrophobicity, pI, molecular weight differences

Pathogenicity Scores:
  - CADD_phred, REVEL, gMVP, ESM, AlphaMissense, MisFit_S, MisFit_D

Clustering Features:
  - Molecular cluster (4 groups, hierarchical)
  - Patient trajectory cluster (GBTM: 2-4 groups)
  - PCA cluster (k-means: 2-4 groups, molecular or molecular+clinical)
