# Inorganic arsenic in beef cattle — population PBPK + risk assessment (concentration-based)

This repository contains the R code used to derive a health-based **Maximum
Limit (ML)** for inorganic arsenic (iAs) in beef cattle feed. This is the
**concentration-based** variant of the analysis: the risk assessment is driven
directly by the **population distribution of steady-state tissue iAs
concentrations** predicted by the PBPK model.

1. A **physiologically-based pharmacokinetic (PBPK) model** that simulates
   AsIII / AsV / MMA / DMA distribution in cattle tissues (Hung 2021 control
   diet) and yields, per virtual animal, the steady-state iAs (AsIII + AsV)
   concentration in muscle, liver, kidney, and other offal.
2. A **Monte Carlo risk assessment** that resamples these per-animal
   steady-state tissue concentrations and back-calculates a cancer-based feed
   ML for adult consumers, protective at the 1-in-1,000,000 lifetime risk
   level.

---

## Repository layout

```
.
|-- scripts/
|   |-- PBPK_ML_cattle_Adult_concentration.R   # Integrated pipeline (concentration-based)
|   `-- check_steady_state.R                   # Time-to-steady-state diagnostic
|-- data/
|   |-- food_intake_iAs.xlsx                    # Inputs: food intake, bodyweight, CSF
|   `-- parameter_conversion.csv                # Human -> cattle kinetic parameter scaling
|-- output_concentration/                       # Generated artefacts (git-ignored)
|-- README.md
|-- LICENSE
`-- .gitignore
```

The script auto-detects its own location and reads `data/` / writes
`output_concentration/` relative to the repository root.

### Cross-species parameter scaling

Kinetic parameters are extrapolated from a validated human arsenic PBPK model
(El-Masri & Kenyon 2008; Yu 1999; Mann et al. 1996) to cattle (human 70 kg ->
cattle 621 kg). Following standard allometric theory, **first-order rate
constants** (absorption `Ka`, redox `K_red`/`K_ox`, urinary/biliary/faecal
excretion) scale as **BW^-0.25**, **maximum metabolic velocities** (`Vmax`)
scale as **BW^0.75**, and **Michaelis-Menten constants** (`Km`) are invariant
(`BW^0`). The full per-parameter mapping (human value, scaling rule, factor,
cattle value) is tabulated in `data/parameter_conversion.csv`.

---

## Requirements

- R >= 4.2
- R packages: `deSolve`, `tidyverse`, `parallel`, `gridExtra`, `readxl`,
  `MASS`, `truncnorm`, `scales`, `conflicted`; plus `sensitivity` for the
  sensitivity-analysis script.

Install:

```r
install.packages(c("deSolve", "tidyverse", "parallel", "gridExtra",
                   "readxl", "MASS", "truncnorm", "scales", "conflicted",
                   "sensitivity"))
```

---

## How to reproduce

From a terminal at the repository root:

```bash
Rscript scripts/PBPK_ML_cattle_Adult_concentration.R
```

Or from within RStudio: open the script and `source()` it. Runtime is a few
minutes (PBPK population on multiple cores, then the risk Monte Carlo). The
PBPK section runs in parallel (`parallel::makeCluster` on Windows,
`mclapply` on Unix).

### Steady-state diagnostic

After the main pipeline has produced the deterministic time course, check when
tissue iAs reaches steady state:

```bash
Rscript scripts/check_steady_state.R
```

It reports, per tissue, the time to reach 90/95/99/99.9% of the plateau and a
rate-based time (when the relative rate of change falls below 1.0 / 0.5 /
0.1 %/day). A custom path to `A_pbpk_deterministic.csv` can be passed as an
argument.

### Supplementary analyses

Two further scripts reuse the model definitions from the main pipeline (they
re-evaluate only the definition prefix, so the full population run is **not**
repeated) and write their outputs to `output_concentration/`:

```bash
# Sensitivity analysis: local one-at-a-time (normalized SC) + global Morris
# screening of the 48 kinetic/partition parameters. Requires the `sensitivity`
# package. Outputs: SA_sensitivity_coefficients.csv, SA_sensitivity.{tiff,png}
Rscript scripts/sensitivity_analysis.R

# Linear-system / eigenvalue analysis (Appendix A.1): finite-difference Jacobian,
# eigenvalue spectrum (relaxation rate constants), and the eigenvector species
# composition of the slow modes.
# Outputs: S1_eigen_summary.csv, S1_linear_system.{tiff,png}
Rscript scripts/linear_system_analysis.R
```

---

## Workflow inside the main script

| Part | Steps | Outputs |
|------|-------|---------|
| **A. PBPK** | Deterministic ODE + mass balance + 10,000-run population Monte Carlo | `A_*` CSVs and figures in `output_concentration/` |
| **B. Concentrations** | Build per-animal steady-state tissue iAs concentration table | in-memory inputs for Part C |
| **C. Risk / ML** | 10,000-trial Monte Carlo cancer-risk model using steady-state concentrations | `C_*` CSVs and figures |

---

## License

Released under the MIT License (see `LICENSE`). The repository will be archived
on Zenodo with a permanent DOI upon acceptance.
