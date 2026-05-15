# Repro Samples — Short Course Demo

This repository accompanies a short course on **simulation-based statistical inference**, covering three methods built on the Repro Samples framework and one method from the Simulation-Based Inference (SBI) literature.

The core idea across all methods: instead of deriving analytical confidence sets, we generate artificial samples from the data-generating process and invert a test statistic to obtain valid frequentist confidence sets.

---

## Repository Structure

```
.
├── demo/
│   ├── repro_demo_presentation.Rmd   # Main R Markdown demo (GMM + HD Linear + Classification)
│   ├── lf2i_demo.ipynb               # LF2I / Waldo demo (Python)
│   └── sbi_nre_demo.ipynb            # Classifier-based SBI / NRE demo (Python)
│
├── gmm/
│   └── functions.R                   # Repro Samples for GMM with unknown K
│
├── hd_linear/
│   └── functions.R                   # Repro Samples for high-dimensional linear regression
│
├── classification/
│   └── functions.R                   # Repro Samples for high-dimensional binary classification
│
├── lf2i-main/                        # LF2I package source (Waldo method)
│
└── paper/                            # Reference papers
```

---

## Methods Covered

### 1. GMM with Unknown Number of Components (`gmm/`)

**Setting:** Mixture of Gaussians $Y \sim \sum_k \pi_k \mathcal{N}(\mu_k, \sigma^2_k)$ with unknown $K$.

**Repro Samples approach:**
- Generate artificial samples $Y^* = G(\theta, U^*)$ using the same noise structure
- Estimate $K$ via modified BIC on a flexmix regression of $[Y \mid U^*]$
- Nuclear mapping: tail probability of the $\hat{K}$ distribution
- Confidence set: all $K$ values not rejected by the nuclear test

### 2. High-Dimensional Linear Regression (`hd_linear/`)

**Setting:** $Y = X\beta + \varepsilon$, $p \gg n$, unknown active set $\tau$ and coefficients $\beta_\tau$.

**Repro Samples approach:**
- Fisher inversion via adaptive Lasso on $[X \mid U^*]$
- Model confidence set via Fisher–Dempster $p$-values
- Coefficient confidence intervals for each $\beta_j$

### 3. High-Dimensional Binary Classification (`classification/`)

**Setting:** Logistic regression with $p \gg n$, unknown active set $\tau$.

**Repro Samples approach:**
- Fisher inversion via logistic adaptive Lasso
- Waldo-like nuclear mapping over the discrete model space
- Coefficient confidence intervals via betaj_cs_wald

### 4. LF2I / Waldo (`demo/lf2i_demo.ipynb`)

**Setting:** Simulation-based inference for $X \mid \theta \sim \mathcal{N}(\theta, \sigma^2)$.

**Waldo approach (regression-based SBI):**
- Train ML regressors to estimate $\mathbb{E}[\theta \mid X]$ and $\mathrm{Var}[\theta \mid X]$
- Test statistic: $T(\theta, X) = (\hat{\mu}(X) - \theta)^2 / \hat{v}(X)$
- Calibrate critical values via quantile regression
- Confidence set via Neyman inversion

### 5. Neural Ratio Estimation (`demo/sbi_nre_demo.ipynb`)

**Setting:** Same Gaussian model with analytical solution for verification.

**NRE approach (classifier-based SBI):**
- Train a binary classifier to distinguish joint samples $p(\theta, X)$ from the marginal product $p(\theta)p(X)$
- Classifier output estimates the likelihood ratio: $\hat{r}(\theta, x) = d / (1 - d)$
- **Bayesian route:** posterior $\propto \hat{r}(\theta, x) \cdot p(\theta)$ → 90% HPD interval
- **Frequentist route:** $T(\theta, x) = -\log \hat{r}(\theta, x)$ → Neyman inversion → 90% CS
- Results verified against analytical solutions

---

## Running the Demos

### R Markdown demo (Methods 1–3)

Requires R with packages: `ClusterR`, `flexmix`, `dplyr`, `glmnet`, `MASS`, `tidyverse`, `pbapply`, `intervals`, `mixAK`, `ggplot2`.

```r
# In RStudio or VS Code with R extension:
rmarkdown::render("demo/repro_demo_presentation.Rmd")
```

### Python demos (Methods 4–5)

Requires Python 3.10+ with: `numpy`, `scipy`, `matplotlib`, `scikit-learn`, `torch`, `xgboost`.

```bash
# Install lf2i (for lf2i_demo.ipynb)
pip install -e lf2i-main/

# Open in Jupyter
jupyter notebook demo/lf2i_demo.ipynb
jupyter notebook demo/sbi_nre_demo.ipynb
```

The NRE demo (`sbi_nre_demo.ipynb`) has no lf2i dependency and runs with standard packages only.

---

## Reference Papers

See `paper/` for the corresponding methodology papers.
