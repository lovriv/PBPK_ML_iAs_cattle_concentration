# Inorganic arsenic in beef cattle — integrated PBPK + MRL risk assessment

This repository contains the R code used to derive a Maximum Residue Limit
(MRL) for inorganic arsenic (iAs) in beef cattle feed via:

1. A **physiologically-based pharmacokinetic (PBPK) model** that simulates
   AsIII / AsV / MMA / DMA distribution in cattle tissues (Hung 2021 control
   diet).
2. A **Monte Carlo risk assessment** that converts the PBPK-derived tissue
   transfer factors (TFs) into a cancer-based MRL for Adult consumers,
   protective at the 1-in-1,000,000 lifetime risk level.

> **Status:** under peer review — private repository.

---

## Repository layout

```
.
|-- scripts/
|   `-- PBPK_MRL_cattle_Adult.R    # Single integrated pipeline
|-- data/
|   `-- food_intake_iAs.xlsx       # Inputs: food intake, bodyweight, CSF
|-- output/                         # Generated artefacts (git-ignored)
|-- README.md
|-- LICENSE
`-- .gitignore
```

The script is fully self-contained: it auto-detects its own location and
reads `data/` / writes `output/` relative to the repository root.

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
Rscript scripts/PBPK_MRL_cattle_Adult.R
```

Or from within RStudio: open the script and `source()` it.

Runtime: ~6 minutes on 7 cores (PBPK ~5 min, MRL ~1 min). The PBPK section
uses parallel execution; on Windows it spawns a `parallel::makeCluster`,
on Unix it uses `mclapply`.

Optional: set `N_RUNS` to change the PBPK Monte Carlo size (default 1000):

```bash
N_RUNS=5000 Rscript scripts/PBPK_MRL_cattle_Adult.R
```

---

## Workflow inside the script

| Part | Steps | Outputs |
|------|-------|---------|
| **A. PBPK** | Deterministic ODE + mass balance + 1000-run Monte Carlo | `A_*` CSVs and figures in `output/` |
| **B. TF bridge** | Compute TF mean + 95% CI from PBPK population | `B_TF_summary.csv` |
| **C. MRL** | 10000-trial Monte Carlo cancer risk model | `C_*` CSVs and figures |

Console summary reports:
- PBPK mass balance (two independent checks, both should pass within
  `stopifnot` tolerances)
- TF point estimates and 95% CI per tissue
- Final MRL recommendation in ug iAs / kg feed dry matter

---

## Key outputs

Generated under `output/` after running the script:

- `A4_Cattle_Tissue_TimeCourse.tiff` — publication-grade time course (Fig. 1)
- `A5_Cattle_SteadyState_Distribution.tiff` — population distribution (Fig. 5)
- `C1_MRL_distribution.png` — Monte Carlo MRL histogram (Fig. 4)
- `C_Table4_MRL_by_endpoint.csv` — Table 4 raw values
- `C_regulatory_summary.csv` — final regulatory recommendation

---

## Data source

`data/food_intake_iAs.xlsx` contains three sheets:

- **food intake** — population-level food intake (g/day) per age group,
  per tissue (meat / liver / kidney / others), per population (general
  public / consumer only). Source: Taiwan National Food Consumption Survey.
- **bodyweight** — bodyweight (kg) distribution per age group.
- **csf** — cancer slope factors for skin, lung, bladder
  (units: (ug/kg-day)^-1).

PBPK control dose corresponds to the Ching-Chi Hung (2021) iAs control diet
fed to cattle (9420 g DM/day, 103 ug iAs/kg DM = AsIII 14.3 + AsV 89).

---

## Citation

If you use this code, please cite the accompanying manuscript (forthcoming).
A DOI for the archived release will be issued via Zenodo upon acceptance.

---

## License

See [LICENSE](LICENSE) (MIT).
