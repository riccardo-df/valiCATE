# Validation of conditional average treatment effects models <a href="https://riccardo-df.github.io/valiCATE/"><img src="man/figures/logo.svg" align="right" style="height:130px;"/></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) [![CRAN](https://www.r-pkg.org/badges/version/valiCATE)](https://CRAN.R-project.org/package=valiCATE) [![Downloads](https://cranlogs.r-pkg.org/badges/valiCATE)](https://CRAN.R-project.org/package=valiCATE) 

The `valiCATE` package validates machine learning predictions of Conditional Average Treatment Effects (CATEs) using the Centered-Weighted Average Treatment Effect (CWATE) and its normalized counterpart (NCWATE). These estimands unify existing validation tools — BLP, AUTOC, QINI, and AUC-HVL — as special cases of a single framework. AIPW-based estimation with cross-fitted nuisance functions provides valid inference under mild regularity conditions.

------------------------------------------------------------------------

## Why use `valiCATE`?

| Feature                             | Benefit                                                                                                                                        |
|-------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------|
| **Heterogeneity detection**         | CWATE one-sided test: does estimated CATE variation reflect genuine treatment heterogeneity or estimation noise?                                |
| **CATE recovery**                   | NCWATE two-sided test: do predicted CATEs recover the true treatment effect function?                                                          |
| **Unified framework**               | Four built-in weight functions (AUTOC, AUC-HVL, BLP, QINI) nest existing tools as special cases of a single estimand.                         |

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

> **Di Francesco, R., & Knaus, M. C. (2025).** *Validating Machine Learning Predictions of Heterogeneous Treatment Effects via Centered-Weighted Average Treatment Effects.* arXiv, 2025.
