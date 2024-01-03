# bootstrap-ci

***Toolbox for bootstrap sampling and estimation of confidence intervals.***

You can choose between hierarchical and non-parametric sampling and combine them 
with multiple bootstrap methods for estimation of confidence intervals.

This initial version is adapted for a research paper on bootstrap sampling. Next version will be adjusted for ease of use.

## Table of Contents
- [Getting Started](#getting-started)
- [Bootstrap sampling](#bootstrap-sampling)
- [Boostrap methods](#bootstrap-methods)
- [Parameters](#parameters)
- [Suggestions](#suggestions-on-which-method-and-parameters-to-use)

# Getting started
Installation and a simple use case example.

## Installation
To use the `bootstrap-ci` package you will need to download it from [pip](https://pypi.org/project/pip/): 
```
pip install bootstrap-ci
```
## Simple example
Once you installed the package, you can use `ci` method to obtain confidence intervals for your chosen
statistic on a given sample:
```
import bootstrap-ci as boot
import numpy as np

np.random.seed(0)
sample = np.random.normal(0, 1, size=1000)

bootstrap = boot.Bootstrap(sample, statistic=np.mean)

onesided_95 = bootstrap.ci(coverages=[0.95], nr_bootstrap_samples=1000)
print(f'One-sided 95% confidence interval for mean is equal to (-inf, {round(onesided_95[0], 3)}).')

>>> One-sided 95% confidence interval for mean is equal to (-inf, 0.004).

twosided_95 = bootstrap.ci(coverages=0.95, side='two', nr_bootstrap_samples=1000)
print(f'Two-sided 95% confidence interval for mean is equal to ({round(twosided_95[0], 3)}, {round(twosided_95[1], 3)}).')

>>> Two-sided 95% confidence interval for mean is equal to (-0.108, 0.014).
```

To see more examples for different sampling possibilities go to [Parameters](#parameters).

# Bootstrap sampling
Bootstrap can be divided into two separate steps. The first one is **bootstrap sampling**, that produces the 
bootstrap distribution, which approximates the distribution of the observed parameter. 
There are different approaches to bootstrap sampling, differing primarily in their underlying data assumptions and 
parameter estimation. In this package you can choose between non-parametric and hierarchical sampling. 

### Non-parametric sampling
Non-parametric sampling is assumption-free and estimates the underlying data distribution $F$ directly with 
the original sample $X$. 
This means that for each bootstrap sample, it samples with replacement directly from the original sample. 
There are $n^n$ different possible samples that can arise with such procedure, but because of computational 
intensiveness, you can choose the number of independent samples, $B$, that you want to obtain. 
To obtain the bootstrap distribution, the value of the observed statistic is calculated on each of them.

### Hierarchical sampling
Hierarchical bootstrap sampling takes into account the group dependencies of the underlying data generating process.
We implemented the completely non-parametric *cases sampling*, where you can choose between all possible strategies,
and the parametric *random-effect sampling*.

#### Cases sampling
Bootstrap samples are obtained by resampling the groups on each level. 
They can be resampled with or without replacement, the latter meaning that we just take all the groups (or data points) 
on that level. Sampling strategy is selected with a vector of zeros and ones, $s = (s_1, \dots, s_{n_{lvl}})$, 
of the same length as the number of levels in the sample. The value 1 in the vector denotes sampling with replacement 
from that particular level. The value of 0 denotes sampling without replacement for that level.

#### Random-effect sampling
Random-effect sampling is a parametric sampling method that assumes that the data come from a random-effect model.
It first estimates the random effects of each group on each level of the sample, then draws those random effects with
replacement, to produce new bootstrap samples.

# Bootstrap methods
After bootstrap sampling, you can use one of the **bootstrap methods** to construct a confidence interval from the 
acquired bootstrap distribution. 

### Percentile
The percentile method is the original bootstrap method. 
Even though multiple improvements were made, it is probably still the most used one.
The percentile estimation of confidence level $\alpha$ is obtained by taking the $\alpha$ quantile of the bootstrap 
distribution,

$$\hat{\theta}\_{perc}\[\alpha\] = \hat{\theta}^*_\alpha.$$

In all the implementations of methods that use quantiles, the "median-unbiased" version of quantile calculation is used.

### Standard
The standard method, sometimes also called the normal method, assumes that the bootstrap distribution is normal and 
estimates standard deviation based on that. The estimations of confidence levels are obtained with

$$\hat{\theta}\_{std}\[\alpha\] = \hat{\theta} + \hat{\sigma} z_\alpha,$$
where $\hat{\theta}$ is the parameter value on the original sample, $\hat{\sigma}$ is the standard deviation estimate 
from the bootstrap distribution and $z_\alpha$ is the z-score of standard normal distribution.

### Basic
In the basic method, also sometimes called the reverse percentile method, the observed bootstrap distribution, 
$\theta^\*$, is replaced with $W^\* = \theta^* - \hat{\theta}$. This results in 
$$\hat{\theta}\_{bsc}\[\alpha\] = 2\hat{\theta} - \hat{\theta}^*\_{1 - \alpha}.$$

### BC
$BC$ does an important correction to the percentile interval. It removes the bias that arises from $\hat{\theta}$ 
not being the median of the bootstrap distribution, and is thus better in non-symetric problems, 
where the percentile method can fail.
The confidence level is estimated by:

$$\hat{\theta}\_{BC}\[\alpha\] = \hat{\theta}^*\_{\alpha_{BC}}, $$

$$\alpha_{BC} = \Phi\big(2\Phi^{-1}(\hat{b}) + z_\alpha \big),$$

where $\Phi$ is the CDF of standard normal distribution and $\hat{b}$ is the bias, calculated as the percentage of 
values from bootstrap distribution that are lower than the parameter's value on the original sample, $\hat{\theta}$.

### BC<sub>a</sub>
$BC_a$ does another correction to the $BC$ interval, by computing the acceleration constant $a$, which can account 
for the skewness of the bootstrap distribution.

This further adjusts the $\alpha_{BCa}$, which is then calculated by:

$$ \hat{\theta}\_{BCa}\[\alpha\] = \hat{\theta}^*\_{\alpha_{BCa}}$$

$$\alpha_{BCa} = \Phi\Big(\Phi^{-1}(b) + \frac{\Phi^{-1}(\hat{b}) + z_\alpha}{1 + \hat{a} (\Phi^{-1}(\hat{b}) + 
z_\alpha)} \Big),$$
where $\hat{a}$ is the approximation of the acceleration constant, that can be calculated using leave-one-out jackknife:

$$\hat{a} = \frac{1}{6}\frac{\sum U_i^3}{(\sum U_i^2)^\frac{3}{2}} $$
$$U_i = (n-1)(\hat{\theta}\_. - \hat{\theta}\_{(i)}),$$
where $\hat{\theta}\_{(i)}$ is the estimation of $\theta$ without the $i$-th datapoint and $\hat{\theta}\_.$ is the mean 
of all $\hat{\theta}_{(i)}$.

### Smoothed
The smoothed method replaces bootstrap distribution with a smoothed version of it ($\Theta^*$), by adding random noise, 
with a normal kernel centered on 0. 
The kernel's size is determined by a rule of thumb width selection: 
$h = 0.9 \min \big( \sigma^\*, \frac{iqr}{1.34} \big),$
where $iqr$ is the inter-quartile range of bootstrap distribution, the difference between its first and third quartile.

The estimation of the confidence level is then obtained by taking the $\alpha$ quantile of the smoothed distribution:

$$\hat{\theta}\_{smooth}\[\alpha\] = \hat{\Theta}^\*\_\alpha.$$

### Studentized
The studentized or bootstrap-t method, generalizes the Student's t method, using the distribution of 
$T = \dfrac{\hat{\theta} - \theta}{\hat{\sigma}}$ to estimate the confidence level $\alpha$.
It is computed by
$$\hat{\theta}\_{t}\[\alpha\] = \hat{\theta} - \hat{\sigma} T\_{1-\alpha},$$
where $\hat{\theta}$ is the parameter value on the original sample, $\hat{\sigma}$ is the standard deviation estimate 
from the bootstrap distribution.
Since the distribution of T is not known, its percentiles are approximated from the bootstrap distribution.
That is done by defining $T^\* = \dfrac{\hat{\theta}^\* - \hat\theta}{\hat{\sigma}^\*}$, where $\hat{\theta}^\*$ is the 
parameter's value on each bootstrap sample, and $\hat{\sigma}^\*$ is obtained by doing another inner bootstrap sampling 
on each of the outer samples. There are other possible ways to acquire $\hat{\sigma}^\*$, but we chose this way as it is 
very general and fully automatic.

### Double
The double bootstrap is made to adjust bias from a single bootstrap iteration with another layer of bootstraps.
The bootstrap procedure is repeated on each of the bootstrap samples to calculate the bias - the percentage of times 
that the parameter on its inner bootstrap sample is smaller from the original parameter's value. 
We want to take such a limit that $P \{\hat{\theta} \in (-\infty, \hat{\theta}\_{double}\[\alpha\])\} = \alpha$, 
which is why we need to select the $\alpha$-th quantile of biases $\hat{b}^*$ for the adjusted level $\alpha_{double}$. 
This leads to:

$$\hat{\theta}\_{double}\[\alpha\] = \hat{\theta}^\*\_{\alpha\_{double}}$$ 
$$\alpha_{double} = \hat{b}^\*_\alpha.$$

# Parameters
Here we describe the possible parameter values on different steps and present some additional examples.

## Initialization
First a `Bootstrap` instance needs to be initialized. Following parameters can be set:
- `data`: a `numpy` array containing values of the sample of interest.
- `statistic`: a callable function that accepts arrays of the same structure as parameter `data` and return a single 
  value.
- `use_jit`: bool that selects whether to use the `numba` library to speed up the sampling. Default value is set to 
  `False`. Change to `True` if you use a big number of bootstrap samples and want to speed up the calculations 
- `group_indices`: a parameter given only for hierarchical data. A list of lists that tells us how the data points in
`data` parameter group together. For example indices \[\[\[0, 1], \[2]], \[\[3]]] together with array \[0, 1, 2, 3] tell 
  us we have one group containing a group with points 0 and 1, and a group with point 2, and another group containing 
  a group with point 3.

You initialize an instance that will estimate the distribution of mean statistic on a given sample from normal 
distribution with the following code:
```
import bootstrap-ci as boot
import numpy as np

np.random.seed(0)
sample = np.random.normal(0, 1, size=1000)

bootstrap = boot.Bootstrap(sample, statistic=np.mean)
```

## Sampling

The method `sample` draws bootstrap samples from the original dataset. Following parameters can be used:
- `nr_bootstrap_samples`: how many bootstrap samples to draw, the size of the bootstrap distribution. Default value is 
  set to 1000, but we propose to take the largest feasible number to get the best results.
- `seed`: random seed. Default value `None` skips setting the seed value.
- `sampling`: select the type of sampling - possible to choose between *nonparametric* or *hierarchical* sampling, default
  value is *nonparametric*.
- `sampling_args`: sampling arguments, used only when doing hierarchical sampling. They should be saved in a dictionary,
  that should include key *method*. Implemented methods available to choose from are *cases* and *random-effect*. 
  For *cases* sampling a *strategy* also needs to be defined with an array of equal length as is the number of levels 
  in the dataset, containing zeroes and ones, telling us on which level we sample with replacement and where without.

For example, you can use the non-parametric sampling to get bootstrap distribution of size 1000 on the `bootstrap`
instance from above.
```
bootstrap.sample(nr_bootstrap_samples=1000, seed=0)
```

The values of the bootstrap distribution are now saved in the `bootstrap.bootstrap_values` parameter.

### Hierarchical sampling
If you are working with hierarchical data, you need to specify the group structure together with the given sample.
There are two different hierarchical sampling methods available to choose from, *random-effect* and *cases* sampling.
Here is an example of *cases* sampling where we sample with replacement on all but the last level:
```
# sample that is grouped like this: [[[0.1, -0.2], [1, -0.5]], [[10, 11]]]
sample = np.array([0.1, -0.2, 1, -0.5, 10, 11])
indices = [[[0, 1], [2, 3]], [[4, 5]]]

hierarchical_bootstrap = boot.Bootstrap(sample, statistic=np.mean, group_indices=indices)

samp_args = {'method': 'cases', 'strategy': [1, 1, 0]}
hierarchical_bootstrap.sample(nr_bootstrap_samples=1000, sampling='hierarchical', sampling_args=samp_args)
```

## Confidence intervals
After the bootstrap distribution is obtained, you can produce the confidence intervals by calling the method `ci`.
Following parameters can be set:
- `coverages`: array of coverage levels for which the values need to be computed. In the case of two-sided
                            intervals (`side`=*two*) it is a float number.
- `side`: it is possible to choose between *one* and *two* sided confidence intervals. One-sided returns
                        the left-sided confidence interval threshold x, representing CI in the shape of (-inf, x).
- `method`: which method to use for construction of confidence intervals. It is possible to select from
                       *percentile*, *basic*, *bca*, *bc*, *standard*, *smoothed*, *double* and *studentized*.
- `nr_bootstrap_samples`: number of bootstrap samples. Default value `None` should be used if the sampling was done 
  before as a separate step and you don't want to repeat it. If the sampling was not done you should specify the 
  number of samples.
- `seed`: random seed. Default value `None` skips setting the seed value.
- `sampling`: type of sampling, possible to choose between *nonparametric* and *hierarchical*. 
  Passed to the method `sample`.
- `sampling_args`: additional arguments used with hierarchical sampling, passed to the method `sample`.
- `quantile_type`: type of quantiles, possible to select from methods used in `numpy`'s quantile function.

It returns an array of threshold values for confidence intervals of corresponding coverage levels.
An example to get the one-sided 95% and 97.5% confidence intervals from the `bootstrap` instance from above, 
where sampling was already done:
```
bootstrap.ci(coverages=[0.95, 0.975], method='bca')

>>> array([0.00272853, 0.0119834 ])
```

## Jackknife after bootstrap
After bootstrap sampling you can diagnose the sampling process with the use of jackknife after bootstrap method, 
that draws a plot showing the influence each data point has on the statistic value.

# Suggestions on which method and parameters to use
For the general use case we propose to use the **double** bootstrap method. In the case of confidence interval of
extreme percentiles, we propose to use the **standard** bootstrap method.

We suggest to always use the largest number of bootstrap samples that is feasible for your sample size and statistic.
If you need to speed up the calculations, lower the number of bootstrap samples from the default value of 1000.

Go to repository [Bootstrap-CI-analysis](https://github.com/zrimseku/Bootstrap-CI-analysis) for more detailed 
information. It includes a detailed study of where bootstrap methods can be used and which one is suggested in 
certain use case.
