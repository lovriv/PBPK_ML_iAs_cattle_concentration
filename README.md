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
|   `-- food_intake_iAs.xlsx                    # Inputs: food intake, bodyweight, CSF
|-- output_concentration/                       # Generated artefacts (git-ignored)
|-- README.md
|-- LICENSE
`-- .gitignore
```

The script auto-detects its own location and reads `data/` / writes
`output_concentration/` relative to the repository root.

---

## Requirements

- R >= 4.2
- R packages: `deSolve`, `tidyverse`, `parallel`, `gridExtra`, `readxl`,
  `MASS`, `truncnorm`, `scales`

Install:

```r
install.packages(c("deSolve", "tidyverse", "parallel", "gridExtra",
                   "readxl", "MASS", "truncnorm", "scales"))
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
