# Reconstruction of Truncated Vertical Profiles using FDA 

This repository provides the R implementation used in Couteyen Carpaye et al. (in prep) to reconstruct truncated vertical profiles (e.g. temperature, salinity...) with functional data analysis (FDA). This method is based on the principal components analysis through conditional expectation (PACE) described by Yao et al. (2005).

The method is designed for substantial datasets of vertical profiles, that may be irregularly sampled, with some profiles being truncated, or partially missing.

Parts of the implementation are adapted from the R package **fda** (Ramsay, 2024), in compliance with its GPL license.

---

## Repository structure

```

.
├── main.R                              # Main script to run the reconstruction
├── data/
│   └── example_data.csv                # Example sparse vertical profiles
├── src/
│   ├── my_fda.R                        # Adapted FDA-related functions
│   ├── profile_reconstruction.R        # Source functions for reconstruction
│   ├── run_pace_on_sparse_profiles.R   # Core FPCA reconstruction workflow
│   ├── visualisation.R                 # Plotting and visualisation functions
│   └── my_app.R                        # Optional Shiny app for parameter tuning
└── README.md

````

---

## Requirements

The following R packages are required:

- `fda`
- `ggplot2`
- `dplyr`
- `shiny` (optional, for interactive visualisation)
- `bslib`, `shinycssloaders` (optional)

---

## Input data format

Input data must be provided as a data frame with at least the following columns:

| Column name | Description |
|------------|-------------|
| `INDEX`    | Profile identifier |
| `DEPTH`    | Vertical coordinate |
| `<VAR>`    | Variable to reconstruct (e.g. `TEMPERATURE`, `SALINITY`) |

Each `INDEX` corresponds to one vertical profile, which may be incomplete or irregularly sampled.

---

## Running the reconstruction

The reconstruction pipeline is executed via the script `main.R`.

Before running the script, please ensure that:

- The main branch of this repository is set as your working directory
- Your input data file is correctly specified in the object `sparse_data`
- The name of the variable to reconstruct is defined in `my_var`
- Reconstruction parameters (e.g. smoothing parameter, number of basis functions, number of PCs) are set according to your application


The main steps of the reconstruction are:

1. Identification of complete and incomplete profiles
2. Construction of a functional data object from complete profiles using B-splines
3. FPCA on functional data object
4. Projection of incomplete profiles in the FPCA space with PACE
5. Reconstruction of full vertical profiles
6. Computation of confidence intervals

---

## Output objects

### `reconstructed_profiles`

`reconstructed_profiles` is a **data frame containing reconstructed vertical profiles for initially incomplete observations**.

Each row corresponds to a reconstructed value at a given depth for a given profile.

#### Structure

| Column name | Description                                                 |
| ----------- | ----------------------------------------------------------- |
| `INDEX`     | Profile identifier                                          |
| `DEPTH`     | Depth at which the reconstruction is evaluated              |
| `<VAR>`     | Reconstructed value of the target variable                  |
| `CI_low`    | Lower bound of the pointwise confidence interval            |
| `CI_up`     | Upper bound of the pointwise confidence interval            |
| `CB_low`    | Lower bound of the confidence band                          |
| `CB_up`     | Upper bound of the confidence band                          |

---

### Other outputs

* `complete_data`: observed complete profiles used to fit the FPCA model
* `fdobj`: Functional data object from complete profiles (`fobj`)
* `cov_fda`: Functional covariance estimate from fdobj 
* `fpca`: FPCA object from complete profiles (`fda::pca.fd`)
* `scr`: Projected scores of incomplete profiles
* `new_fdobj `: Functional data object from complete and reconstructed profiles (`fobj`)
* `all_data`: combined observed complete and reconstructed incomplete profiles

---

## Visualisation

The repository provides helper functions to:

* Visualise FPCA factorial plane (`plot_scores_projected`)
* Visualise observed and reconstructed profiles (`plot_one_profile`)
* Inspect covariance (`plot_cov_fun`) and correlation functions (`plot_cor_fun`)

---

## License

This code is released under the **GNU General Public License v3.0 (GPL-3.0)**.

Parts of this code are adapted from the R package **fda** (Ramsay, 2024).

Any redistribution or modification must comply with the terms of the GPL.

---

## Citation

If you use this code in your research, please cite the `fda` package and both **the associated scientific article** and **this software repository**.

Please cite the software as follows:
TO-DO

Please cite the associated article:
TO-DO


---

## References

Ramsay J. _fda: Functional Data Analysis_. R package version 6.2.0, https://CRAN.R-project.org/package=fda, 2024.

Yao, F., Müller, H. G., and Wang, J. L.: Functional data analysis for sparse longitudinal data, JASA, 100, 577–590, https://doi.org/10.1198/016214504000001745, 2005.


