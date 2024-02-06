# Evaluating Machine Learning Estimated Conditional Average Treatment Effects in RCTs <a href="https://riccardo-df.github.io/evaluCATE/"><img src="man/figures/logo.svg" align="right" height="200" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/riccardo-df/evaluCATE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/riccardo-df/evaluCATE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Quality evaluation of machine learning estimated conditional average treatment effects (CATEs). The quality is assessed by estimating the best linear predictor of the actual CATEs using the estimated CATEs, the sorted group average treatment effects, and the rank-weighted average treatment effects induced by the estimated CATEs. 

To get started, please check the online [short tutorial](https://riccardo-df.github.io/evaluCATE/articles/evaluCATE-short-tutorial.html).

## Installation  
The current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/evaluCATE") # run install.packages("devtools") if needed.
```

## References

- Imbens, G. W., & Rubin, D. B. (2015).
<b>Causal inference for statistics, social, and biomedical sciences: An introduction.</b>
<i>Cambridge University Press</i>.
[<a href="https://www.cambridge.org/core/books/causal-inference-for-statistics-social-and-biomedical-sciences/71126BE90C58F1A431FE9B2DD07938AB">book</a>]

- Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2017).
<b>Generic machine learning inference on heterogeneous treatment effects in randomized experiments.</b>
<i>National Bureau of Economic Research</i>.
[<a href="https://arxiv.org/abs/1712.04802">paper</a>]

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized random forests.</b>
<i>Annals of Statistics</i>.
[<a href="https://projecteuclid.org/journals/annals-of-statistics/volume-47/issue-2/Generalized-random-forests/10.1214/18-AOS1709.full">paper</a>]

- Künzel, S. R., Sekhon, J. S., Bickel, P. J., & Yu, B. (2019).
<b>Metalearners for estimating heterogeneous treatment effects using machine learning.</b>
<i>Proceedings of the National Academy of Sciences</i>.
[<a href="https://www.pnas.org/doi/abs/10.1073/pnas.1804597116">paper</a>]

- Imai, K., & Li, M. (2021).
<b>Statistical inference for heterogeneous treatment effects discovered by generic machine learning in randomized experiments.</b>
<i>arXiv preprint</i>.
[<a href="https://arxiv.org/abs/2203.14511">paper</a>]

- Yadlowsky, S., Fleming, S., Shah, N., Brunskill, E., & Wager, S. (2021).
<b>Evaluating treatment prioritization rules via rank-weighted average treatment effects.</b>
<i>arXiv preprint</i>.
[<a href="https://arxiv.org/abs/2111.07966">paper</a>]
