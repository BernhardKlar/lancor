## lancor: Statistical Inference via Lancaster Correlation

[![CRAN version](https://www.r-pkg.org/badges/version/lancor)](https://cran.r-project.org/package=lancor)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/lancor)](https://cran.r-project.org/package=lancor)


### Synopsis 

The [lancor package](https://cran.r-project.org/web/packages/lancor/index.html) implements the methods described in [Holzmann, Klar (2024)](https://onlinelibrary.wiley.com/doi/full/10.1111/sjos.12733).
Lancaster correlation is a correlation coefficient which equals the absolute value of the Pearson correlation for the bivariate normal distribution,
and is equal to or slightly less than the maximum correlation coefficient for a variety of bivariate distributions. Rank and moment-based estimators and corresponding confidence intervals are implemented, as well as independence 
tests based on these statistics.

### Examples 
 The Lancaster correlation coefficient and the linear Lancaster correlation coefficient are estimated via the `lcor` function:

 ```R
 n <- 1000 
 x <- matrix(rnorm(n*2), n)
 lcor(x, type = "rank")
 lcor(x, type = "linear")
 ```

 Confidence intervals are given by the `lcor.ci` function:

```R
 n <- 1000
 x <- matrix(rnorm(n*2), n)
 nu <- 2
 y <- x / sqrt(rchisq(n, nu)/nu)
 lcor(y, type = "rank")
 lcor.ci(y, type = "rank")
```

Finally the Lancaster correlation test of bivariate independence `lcor.test`:

```R
 n <- 200
 x <- matrix(rnorm(n*2), n)
 nu <- 2
 y <- x / sqrt(rchisq(n, nu)/nu)
 cor.test(y[,1], y[,2], method = "spearman")
 lcor.test(y, type = "rank") 
```


### Installation via CRAN 

The `lancor` package can be installed from within R via

```R
install.packages("lancor")
```

### Authors 

Hajo Holzmann, Bernhard Klar

### License 

GPL-2
