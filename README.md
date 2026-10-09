# uncertainty_time_series

A repository for analysing uncertainty in long-term ecological time series and checking how common statistical modelling assumptions affect trend estimation.

## Project Overview

This project evaluates whether population time-series models can reliably recover long-term trends under realistic levels of temporal autocorrelation, inter-annual variability, and observation duration.

The analysis combines:
- Ecological abundance time series from the BIOTIME database
- Simulated Poisson log-normal time series
- Autoregressive (AR1) mixed models using `glmmTMB`
- Reliability summaries: estimation error, bias, power, false positives, Type M error, and Type S error

## Repository Structure

```
├── data/
│   ├── biotime_table.csv          # Model results from BIOTIME dataset
│   ├── conversor.csv              # Mapping between model IDs and spatial coordinates
│   └── data_filtered.csv          # Filtered and processed BIOTIME data
├── scripts/
│   ├── 0_biotime_clean_data.R     # Data cleaning and filtering
│   ├── 1_analyse_biotime.R        # Fit models to empirical data
│   ├── 2_create-simulation.R      # Generate synthetic time series
│   ├── 3_model-simulation.R       # Fit models to simulated data
│   └── 4_analysis-ms.R            # Generate figures and summaries
├── LICENSE
└── README.md
```

## Script Documentation

### 0_biotime_clean_data.R

**Purpose:** Cleans, filters, and prepares the BIOTIME ecological database for analysis.

**Key Functions:**
- `add_plot_id()`: Clusters spatial sampling locations using hierarchical clustering to create unique plot identifiers per study
- `filter_data()`: Applies multiple quality filters to time series

**Main Steps:**
1. Loads BIOTIME raw data and metadata
2. Filters by taxa (Birds, Mammals, Invertebrates, Reptiles, Amphibians) and terrestrial realm
3. Taxonomically classifies species using GBIF backbone
4. Groups observations into spatial plots within each study
5. Removes time series that:
   - Don't have ≥3 consecutive years of data
   - Have fewer than 2 detections of a species
   - Have <50% detection rate
   - Are all zeros across all years in a study

**Outputs:**
- `data_filtered.csv`: Cleaned dataset ready for modelling
- `conversor.csv`: Lookup table for model IDs to spatial coordinates

**Dependencies:** `reshape2`, `dplyr`, `tidyr`, `ggplot2`, `data.table`, `terra`

---

### 1_analyse_biotime.R

**Purpose:** Fits autoregressive Poisson mixed models to empirical BIOTIME time series and extracts model summaries.

**Model Specification:**
```
ABUNDANCE ~ YEAR2 + ar1(YEAR3 + 0 | group)
family = poisson
```
Includes:
- Linear trend component (YEAR2)
- Autoregressive correlation structure of order 1 (AR1)
- Log-normal random effects for inter-annual variability

**Key Functions:**
- `theta2phi()`: Converts `glmmTMB` theta parameterization to correlation coefficient
- `get_AR()`: Extracts AR1 correlation coefficient with 95% confidence intervals

**Extracted Metrics:**
- Fixed effect estimates and standard errors (intercept, slope)
- Growth rate (r) with confidence intervals
- Standard deviation of random effects (inter-annual variability)
- AR1 autocorrelation coefficient (phi)
- Model comparison (AIC, BIC, log-likelihood)
- Reliability measures (conditional/marginal R²)
- Convergence diagnostics (Hessian positive-definite, convergence flag)

**Note:** Designed for parallel execution using POD_ID environment variable to split models into 100 batches.

**Outputs:** `Model_group_[pod_id].RData` (one per replica)

**Dependencies:** `dplyr`, `data.table`, `glmmTMB`

---

### 2_create-simulation.R

**Purpose:** Generates synthetic Poisson log-normal time series across a grid of ecological parameters.

**Simulation Method:**
- 100 years of observations per replicate
- 100 replicates per parameter combination
- Parameters varied:
  - **Temporal trend (r):** -0.25 to +0.25 by 0.01 (51 levels)
  - **Autocorrelation (phi):** -0.9 to +0.9 by 0.15 (13 levels)
  - **Inter-annual variability (sdev):** 0 to 4 (17 levels)
- Mean abundance (mu_zero): random uniform 0.2–90

**Simulation Steps:**
1. Generates correlated random effects using multivariate normal distribution
2. Creates expected abundance with exponential trend and random effects
3. Draws observation counts from Poisson distribution
4. Stores simulated data and true parameter values

**Mathematical Model:**
```
μ_t = μ₀ × (1 + r)^(t-1) × exp(ε_t)
ε_t ~ MVN(0, σ² × R(φ))
y_t ~ Poisson(μ_t)
```
where R(φ) is the AR1 correlation matrix.

**Outputs:** `simulated_data_[repi].RData` (one per replica)

**Dependencies:** `dplyr`, `mvtnorm`, `glmmTMB`, `MASS`, `data.table`

---

### 3_model-simulation.R

**Purpose:** Fits the same AR1 Poisson model to simulated data and evaluates model performance.

**Key Analysis:**
1. Loads simulated data from script 2
2. For each combination of true parameters, truncates to varying durations (3–57 years)
3. Fits the AR1 Poisson model to the truncated data
4. Extracts estimated parameters and model fit statistics

**Estimated Metrics:**
- Same as script 1 (trend, variance, autocorrelation with CIs)
- Additional: comparison to true simulated values
- Residual variance (full vs. conditional)
- Convergence diagnostics

**Outputs:** `[repi_i]_model.RData` containing results for all parameter combinations and durations

**Dependencies:** `glmmTMB`, `dplyr`, `data.table`

---

### 4_analysis-ms.R

**Purpose:** Synthesizes results and generates publication-ready figures and summary statistics.

**Main Analyses:**

1. **Reliability Metrics on Simulated Data**
   - Filters converged models
   - Computes:
     - Estimation **error** (RMSE of slope)
     - **Bias** in slope estimation
     - **Power** to detect a 2%/year decline
     - **False positive rate** (Type I error)
     - **Type M error** (magnitude exaggeration)
     - **Type S error** (sign reversal)
   - Summarizes by duration, inter-annual variability, and autocorrelation

2. **Empirical BIOTIME Comparison**
   - Loads and filters BIOTIME model results
   - Rounds estimated parameters to nearest simulated scenario
   - Joins with simulated reliability metrics
   - Visualizes expected reliability for each BIOTIME time series

3. **Figure Generation**
   - **Fig 2:** Reliability heatmaps (error, bias, power, false positives, Type M/S errors)
   - **Fig 3:** Distribution of BIOTIME characteristics (duration, variability, autocorrelation)
   - **Fig 4:** Expected vs. empirical reliability in BIOTIME data
   - **Supplementary figures:** Sensitivity analyses, correlation matrices, power by parameter, partial pooling effects

**Key Outputs:**
- `errors_simulations.csv`: Reliability summary across all simulated scenarios
- `error_biotime_filtered.csv`: Reliability metrics for each BIOTIME time series
- `Fig_2.png` through `Fig_S7.png`: Publication figures
- `tableS1.csv`: Supplementary reliability table

**Color Scheme:**
- Insecta: #F2C400, Mammalia: #8B4513, Aves: #ACA0DC, Squamata: #98FB98, Arachnida: #2B0078, Amphibia: #41B6C4

**Dependencies:** `glmmTMB`, `ggplot2`, `dplyr`, `cowplot`, `data.table`, `patchwork`, `scales`, `corrplot`, `tidyverse`

---

## Data Workflow

```
BIOTIME Raw Data
        ↓
0_biotime_clean_data.R  ← Filters, clusters, quality checks
        ↓
data_filtered.csv       ← Clean empirical data
        ↓
1_analyse_biotime.R     ← Fit models to empirical time series
        ↓
biotime_table.csv       ← Model results on real data

        ↓↓↓ SIMULATIONS ↓↓↓

2_create-simulation.R   ← Generate synthetic data across parameter grid
        ↓
simulated_data_*.RData  ← 100 replicates per parameter combo
        ↓
3_model-simulation.R    ← Fit models to simulated data
        ↓
*_model.RData           ← Model results on simulated data
        ↓
4_analysis-ms.R         ← Combine, compare, visualize, produce figures
        ↓
Figures, tables, summaries
```

## Key Analytical Concepts

### Reliability Metrics

- **Estimation Error:** Root mean squared error in estimated growth rate (%)
- **Bias:** Mean difference between estimated and true growth rate
- **Power:** Probability of detecting a significant 2%/year decline when truly present
- **False Positives:** Probability of detecting a significant trend when none exists (r = 0)
- **Type M Error:** Ratio of |estimated effect| to |true effect| among significant results (exaggeration factor)
- **Type S Error:** Probability that sign of estimated effect differs from true effect among significant results

### AR1 Correlation Structure

The autoregressive order-1 structure models temporal dependence:
- **φ > 0:** Positive autocorrelation (similar values in consecutive years)
- **φ = 0:** No autocorrelation
- **φ < 0:** Negative autocorrelation (alternating patterns)

### Inter-Annual Variability

Modeled as log-normal random effects. The standard deviation (sigma) on log-scale is converted to coefficient of variation (CV) for interpretation.

## Software Requirements

**R Packages:**
- Statistical modelling: `glmmTMB`, `MASS`
- Data manipulation: `dplyr`, `data.table`, `tidyr`, `reshape2`
- Visualization: `ggplot2`, `cowplot`, `patchwork`, `scales`, `colorspace`, `corrplot`
- Spatial: `terra`
- Utilities: `here`, `stringr`

**R Version:** 4.0+

## Running the Analysis

### Sequential execution (local machine):

```bash
Rscript scripts/0_biotime_clean_data.R
Rscript scripts/1_analyse_biotime.R
Rscript scripts/2_create-simulation.R
Rscript scripts/3_model-simulation.R
Rscript scripts/4_analysis-ms.R
```

### Parallel execution (HPC cluster):

Scripts 1 and 3 are designed for parallel execution. Set the `POD_ID` environment variable (0–99) to process each replica:

```bash
for i in {0..99}; do
  POD_ID=$i Rscript scripts/1_analyse_biotime.R &
done
wait

for i in {0..99}; do
  POD_ID=$i Rscript scripts/3_model-simulation.R &
done
wait

Rscript scripts/4_analysis-ms.R
```

## Notes

- **Hard-coded paths:** Many scripts contain hard-coded directory paths (e.g., `setwd("/workdir/")`). Adjust these to match your local environment.
- **Data availability:** BIOTIME data requires download from the BIOTIME project website (https://www.biotime.org/).
- **Computational time:** The full simulation and modelling pipeline (especially scripts 2–3 with 100 replicates × ~850 parameter combinations) is computationally intensive and benefits from parallelization.
- **Convergence:** Not all models converge successfully. Scripts filter for convergence using Hessian positive-definiteness and convergence flags.

## License

MIT License. See LICENSE file for details.

## References

For more information on the BIOTIME database, visit: https://www.biotime.org/

For `glmmTMB` documentation, see: https://cran.r-project.org/web/packages/glmmTMB/
