variance d’estimateurs de moyennes
================
Lemariey
2023-07-11

In this example the means of differents batch follow a law N(M,s1). M is
the real mean mi the mean of each batch

Then the observation in a batch follow the law N(mi,s2)

``` r
m <- rnorm(10000,0,2)

n <- 10
mhat <- numeric(length(m))
for ( i in 1:length(m)){
  mhat[i] <- mean(rnorm(n,m[i],1))
}
```

Y= mean estimator

The total variance formula gives:  
Var(Y)=Ex(Var(Y\|X))+Var(E(Y\|X))  
Var(E(mhat sachant m))=s1  
E(Var(mhat sacahnt m))=s2/n  
(s2/n)+ s1

``` r
sd(m)**2
```

    ## [1] 4.031466

``` r
sd(mhat)**2
```

    ## [1] 4.145089
