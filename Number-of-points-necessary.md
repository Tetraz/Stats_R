Number of point necessary
================
Lemariey
2023-07-11

\#test du nombre de points necessaire pour fit un model

``` r
library("EnvStats")
```

    ## Warning: le package 'EnvStats' a été compilé avec la version R 4.2.3

    ## 
    ## Attachement du package : 'EnvStats'

    ## Les objets suivants sont masqués depuis 'package:stats':
    ## 
    ##     predict, predict.lm

create data

``` r
y <- c(1,2,3)
x1 <- c(4,5,7)
x2 <- c(4,4,7)
```

\#on peut fitter mais on a pas l’erreur

``` r
m0<- lm(y[-1]~x2[-1])
summary(m0)
```

    ## 
    ## Call:
    ## lm(formula = y[-1] ~ x2[-1])
    ## 
    ## Residuals:
    ## ALL 2 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   0.6667        NaN     NaN      NaN
    ## x2[-1]        0.3333        NaN     NaN      NaN
    ## 
    ## Residual standard error: NaN on 0 degrees of freedom
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 1 and 0 DF,  p-value: NA

``` r
m1 <- lm(y~x1)
summary(m1)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x1)
    ## 
    ## Residuals:
    ##        1        2        3 
    ## -0.14286  0.21429 -0.07143 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  -1.4286     0.6776  -2.108    0.282
    ## x1            0.6429     0.1237   5.196    0.121
    ## 
    ## Residual standard error: 0.2673 on 1 degrees of freedom
    ## Multiple R-squared:  0.9643, Adjusted R-squared:  0.9286 
    ## F-statistic:    27 on 1 and 1 DF,  p-value: 0.121

``` r
anovaPE(m1)
```

    ## Error in anovaPE(m1): There must be at least two replicate values for at least one value of the predictor variable.

We can fit the model and have the error but not the pure error cause
there is no replicate

``` r
m2 <- lm(y~x2)
summary(m2)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x2)
    ## 
    ## Residuals:
    ##          1          2          3 
    ## -5.000e-01  5.000e-01  1.943e-16 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)  -0.5000     1.5000  -0.333    0.795
    ## x2            0.5000     0.2887   1.732    0.333
    ## 
    ## Residual standard error: 0.7071 on 1 degrees of freedom
    ## Multiple R-squared:   0.75,  Adjusted R-squared:    0.5 
    ## F-statistic:     3 on 1 and 1 DF,  p-value: 0.3333

``` r
anovaPE(m2)
```

    ## Error in anovaPE(m2): Not enough replicate values for predictors relative to degrees of freedom for residual sums of squares.

We can fit the model and have the error. But we can’t have the
decomposition of the error (the can’t evaluate the fit error)

``` r
y <- c(1,2,3,5)
x3 <- c(4,4,6,5)
m3 <- lm(y~x3)
summary(m3)
```

    ## 
    ## Call:
    ## lm(formula = y ~ x3)
    ## 
    ## Residuals:
    ##         1         2         3         4 
    ## -1.00e+00  1.11e-16 -1.00e+00  2.00e+00 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)   -2.000      5.036  -0.397    0.730
    ## x3             1.000      1.044   0.957    0.439
    ## 
    ## Residual standard error: 1.732 on 2 degrees of freedom
    ## Multiple R-squared:  0.3143, Adjusted R-squared:  -0.02857 
    ## F-statistic: 0.9167 on 1 and 2 DF,  p-value: 0.4394

``` r
anovaPE(m3)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## x3             1   2.75    2.75     5.5 0.2566
    ## Lack of Fit    1   5.50    5.50    11.0 0.1864
    ## Pure Error     1   0.50    0.50

In summary  
-pour fitter sans l’erreur, il faut au moins p points qui forment une
matrice de plein rang  
-pour fitter sans pure erreur: p points qui forment une matrice de plein
rang + 1 point pour erreur  
-pour fitter p parameter correctement, il faut au moins p+1 point
different (pour error of fit) + 1 points suppl?mentaire r?pliqu? pure
erreur (avec uniquement 1 repliqu? suppl?mentaire, il n’y a pas assez de
degre de libert? pour estimer la pure erreur)
