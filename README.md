# Validation of conditional average treatment effects models <a href="https://riccardo-df.github.io/valiCATE/"><img src="man/figures/logo.svg" align="right" height="30"/></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![CRAN](https://www.r-pkg.org/badges/version/valiCATE)](https://CRAN.R-project.org/package=valiCATE) [![Downloads](https://cranlogs.r-pkg.org/badges/valiCATE)](https://CRAN.R-project.org/package=valiCATE) 

The `valiCATE` package provides a suite of tools for validating machine learning models estimating Conditional Average Treatment Effects (CATEs). Models are validated by estimating the best linear predictor of the actual CATEs using the estimated CATEs, the sorted group average treatment effects, and the rank-weighted average treatment effects induced by the estimated CATEs.

------------------------------------------------------------------------

## Why use `valiCATE`?

| Feature                             | Benefit                                                                                                                                        |
|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| **Heterogeneity detection** | Tests whether estimated CATE variation reflects genuine treatment heterogeneity or estimation noise.                              |
| **Multiple estimation strategies** | Implements weighted residuals, Horvitz-Thompson, AIPW, and nonparametric procedures.                                 |
| **Graphical analysis**              | Provides intuitive visualizations of CATE estimates, GATES results, and TOC curves.                                                         |

------------------------------------------------------------------------
## 🚀 Installation

To install the latest stable release from CRAN, run:

```         
# install.packages("valiCATE") # Not yet available!
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```         
devtools::install_github("riccardo-df/valiCATE")
```

------------------------------------------------------------------------

## Contributing

We welcome contributions! If you encounter issues, have feature requests, or want to contribute to the package, please follow the guidelines below.

📌 **Report an issue:** If you encounter a bug or have a suggestion, please open an issue on GitHub:
[Submit an issue](https://github.com/riccardo-df/valiCATE/issues)

📌 **Contribute code:** We encourage contributions via pull requests. Before submitting, please:
1. Fork the repository and create a new branch.
2. Ensure that your code follows the existing style and documentation conventions.
3. Run tests and check for package integrity.
4. Submit a pull request with a clear description of your changes.

📌 **Feature requests:** If you have ideas for new features or extensions, feel free to discuss them by opening an issue.

------------------------------------------------------------------------

## Citation

If you use `valiCATE` in your research, please cite the corresponding paper:

> **Author(s).** *Title of Paper.* arXiv, 2025
