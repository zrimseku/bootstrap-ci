# bootstrap-ci

***Toolbox for bootstrap sampling and estimation of confidence intervals.***

You can choose between hierarchical and non-parametric sampling and combine them 
with multiple bootstrap methods for estimation of confidence intervals. 

# Bootstrap sampling
Bootstrap can be divided into two separate steps. The first one is **bootstrap sampling**, that enables us to get the 
bootstrap distribution, which approximates the distribution of the observed parameter. 
There are different approaches to bootstrap sampling, differing primarily in their underlying data assumptions and 
parameter estimation. In this package you can choose between non-parametric and hierarchical sampling. 

### Non-parametric sampling
Non-parametric sampling is assumption-free and estimates the underlying data distribution $F$ directly with 
the original sample $X$. 
This means that for each bootstrap sample, we are sampling with replacement directly from our original sample. 
There are $n^n$ different possible samples that can arise with such procedure, but because of computational 
intensiveness, we limit ourselves on $B$ independent samples. 
On each of them we calculate the value of the observed parameter, which gives us the bootstrap distribution.

### Hierarchical sampling

# Bootstrap methods
After bootstrap sampling, we use one of the **bootstrap methods** to construct a confidence interval from the acquired 
bootstrap distribution. 

### Percentile
The percentile method is the original bootstrap method. 
Even though multiple improvements were made, it is probably still the most used one.
The percentile estimation of confidence level $\alpha$ is obtained by taking the $\alpha$ quantile of the bootstrap 
distribution, which we annotate by 

$$\hat{\theta}_{perc}\[\alpha\] = \hat{\theta}^*_\alpha.$$

In all of our implementations of methods that use quantiles, we used the "median-unbiased" version of quantile calculation.

### Standard
The standard method, sometimes also called the normal method, assumes that the bootstrap distribution is normal and 
estimates standard deviation based on that. We get the estimations of confidence levels with

$$\hat{\theta}_{std}\[\alpha\] = \hat{\theta} + \hat{\sigma} z_\alpha,$$
where $\hat{\theta}$ is the parameter value on the original sample, $\hat{\sigma}$ is the standard deviation estimate 
from the bootstrap distribution and $z_\alpha$ is the z-score of standard normal distribution.

### Basic
In the basic method, also sometimes called the reverse percentile method, the observed bootstrap distribution, 
$\theta^\*$, is replaced with $W^\* = \theta^* - \hat{\theta}$. This results in 
$$\hat{\theta}\_{bsc}\[\alpha\] = 2\hat{\theta} - \hat{\theta}^*\_{1 - \alpha}.$$


### BC
$BC$ does an important correction to the percentile interval. It removes the bias that arises from $\hat{\theta}$ not being the median of the bootstrap distribution, and is thus better in non-symetric problems, where the percentile method can fail.
The confidence level is estimated by:

```math
\hat{\theta}_{BC}[\alpha] = \hat{\theta}^*_{\alpha_{BC}}, \\

\alpha_{BC} = \Phi\big(2\Phi^{-1}(\hat{b}) + z_\alpha \big),
```
where $\Phi$ is the CDF of standard normal distribution and $\hat{b}$ is the bias, calculated as the percentage of values from bootstrap distribution that are lower than the parameter's value on the original sample, $\hat{\theta}$.

### BC<sub>a</sub>
$BC_a$ does another correction to the $BC$ interval, by computing the acceleration constant $a$, which can account for the skewness of the bootstrap distribution.

This further adjusts the $\alpha_{BCa}$, which is then calculated by:
```math
\hat{\theta}_{BCa}\[\alpha\] = \hat{\theta}^*_{\alpha_{BCa}} \\

\alpha_{BCa} = \Phi\Big(\Phi^{-1}(b) + \frac{\Phi^{-1}(\hat{b}) + z_\alpha}{1 + \hat{a} (\Phi^{-1}(\hat{b}) + z_\alpha)} \Big),
```
where $\hat{a}$ is the approximation of the acceleration constant, that can be calculated using leave-one-out jackknife:

```math
\hat{a} = \frac{1}{6}\frac{\sum_{i=1}^n U_i^3}{(\sum_{i=1}^n U_i^2)^\frac{3}{2}} \\

U_i = (n-1)(\hat{\theta}_. - \hat{\theta}_{(i)}),
```
where $\hat{\theta}\_{(i)}$ is the estimation of $\theta$ without the $i$-th datapoint and $\hat{\theta}\_.$ is the mean 
of all $\hat{\theta}_{(i)}$.

### Smoothed
The smoothed method replaces bootstrap distribution with a smoothed version of it ($\Theta^*$), by adding random noise, 
with a normal kernel centered on 0. 
We determined the kernel's size by a rule of thumb width selection: 
$h = 0.9 \min \big( \sigma^\*, \frac{iqr}{1.34} \big),$
where $iqr$ is the inter quartile range of bootstrap distribution, the difference between its first and third quartile.

The estimation of the confidence level is then obtained by taking the $\alpha$ quantile of the smoothed distribution:
$$ \hat{\theta}\_{smooth}\[\alpha\] = \hat{\Theta}^\*\_\alpha. $$

### Studentized
The studentized or bootstrap-t method, generalizes the Student's t method, using the distribution of $T = \dfrac{\hat{\theta} - \theta}{\hat{\sigma}}$ to estimate the confidence level $\alpha$.
It is computed by
$$\hat{\theta}\_{t}\[\alpha\] = \hat{\theta} - \hat{\sigma} T\_{1-\alpha},$$
where $\hat{\theta}$ and $\hat{\sigma}$ are calculated as described above.
But since the distribution of T is not known, we need to approximate its percentiles from the bootstrap distribution.
We do that by defining $T^* = \dfrac{\hat{\theta}^* - \hat\theta}{\hat{\sigma}^*}$, where $\hat{\theta}^*$ is the parameter's value on each bootstrap sample, and $\hat{\sigma}^*$ is obtained by doing another inner bootstrap sampling on each of the outer samples. There are other possible ways to acquire $\hat{\sigma}^*$, but we chose this way as it is very general and fully automatic.

### Double
The double bootstrap is made to adjust bias from a single bootstrap iteration with another layer of bootstraps.
We repeat the bootstrap procedure on each of the bootstrap samples to calculate the bias -- the percentage of times that the parameter on its inner bootstrap sample is smaller from the original parameter's value. We want to take such a limit that $P \{\hat{\theta} \in (-\infty, \hat{\theta}_{double}[\alpha])\} = \alpha$, which is why we need to select the $\alpha$-th quantile of biases $\hat{b}^*$ for the adjusted level $\alpha_{double}$. This leads to:
\begin{align*}
\hat{\theta}\_{double}\[\alpha\] &= \hat{\theta}^*\_{\alpha\_{double}} \\ 
\alpha_{double} &= \hat{b}^*_\alpha.
\end{align*}


# Suggestions on which method and parameters to use
General double
Percentiles standard

Lower B if you need faster, we propose 1000 for the lowest B for double
Use BCa if you need even faster method for mean...

Repository [Bootstrap-CI-analysis](https://github.com/zrimseku/Bootstrap-CI-analysis) for more detailed information.

# Example of use

### Package download


### Producing confidence intervals

```
sampple iz normalne
Bootstrap
rezultati

```