
# GpGp4RG

GpGp4RG is an R package for fast approximate Gaussian process computation, 
adapted to be compatible with RobustGaSP2 for use in parallel partial emulation (PPE).
The changes from the original package GpGp are summarized below.
The package includes implementations of the Vecchia's (1988) original 
approximation, as well as several updates to it, including the reordered 
and grouped versions of the approximation outlined in Guinness (2018).

## Changes from GpGp

The likelihood function has been changed from the profiled likelihood to
the integrated likelihood as detailed by the accompanying paper.

The Matérn functions have been changed to be consistent with those in RobustGaSP. 
Additionally, the correlation matrix $\mathbf{R}$ has been changed to be computed by 
multiplying across input dimensions. Below is an example of the Matérn 3/2 function.
Similar changes have been implemented for the Matérn 5/2 and power exponential functions.

Previous:
$M(x,y) = \sigma^2 (1 + || D^{-1}(x - y) || ) exp( - || D^{-1}(x - y) || )$
Current:
$M(x,y) = \sigma^2 \prod_{i=1}^p (1 + \sqrt{3}|| (x - y)/\lambda_i || ) exp( - \sqrt{3}|| (x - y)/\lambda_i || )$

The gradient functions used in VPPE have been changed to compute the partial derivatives
with respect to the inverse range parameters to be consistent with RobustGaSP,
rather than the derivatives of the range parameters. By the chain rule, we simply multiply the original
derivatives (with respect to $\lambda$) by $-1/\lambda^2$.

## Installing

The package can be installed from GitHub for the latest version

```{r}
devtools::install_github("sideofman/GpGp4RG")
```

## Basic Use

The main function for fitting models is called 'fit_model', and the
main function for doing predictions is called 'predictions'.

See this youtube video for a tutorial:
https://www.youtube.com/watch?v=phyB4n0CDWg&t=4s
