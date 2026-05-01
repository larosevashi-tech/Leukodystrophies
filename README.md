BSA-curve-analysis/
│
├── data/
│   ├── raw/                    # original data, never touch this
│   └── processed/              # cleaned data ready for analysis
│
├── scripts/
│   ├── 01_clean_data.R         # clean & prepare data
│   ├── 02_linear_model.R       # fit and plot linear curve
│   ├── 03_quadratic_model.R    # fit and plot quadratic curve
│   ├── 04_model_comparison.R   # compare models, pick the best
│   └── 05_extract_unknowns.R   # calculate unknown concentrations
│
├── outputs/
│   ├── plots/                  # saved graphs
│   └── results/                # tables, unknown values
│
├── README.md                   # project description
└── BSA-curve-analysis.Rproj    # RStudio project file
