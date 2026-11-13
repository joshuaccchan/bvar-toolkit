# How to Estimate a VAR by Gibbs Sampling

*Code: [`guide.m`](guide.m), with
[`bvar.samplers.var_coef_joint`](../../core/+bvar/+samplers/var_coef_joint.m),
[`bvar.diag.inefficiency_factor`](../../core/+bvar/+diag/inefficiency_factor.m),
[`bvar.diag.mcse`](../../core/+bvar/+diag/mcse.m) and
[`bvar.diag.geweke`](../../core/+bvar/+diag/geweke.m). Method:
[Chan (2027)](../../CITING.md#chan-2027) and
[Geweke (1992)](../../CITING.md#geweke-1992).*

A Gibbs sampler draws the unknowns of a model in blocks, each from its distribution given the data
and the current values of the other blocks, and repeats. In
[How to Draw the VAR Coefficients under the Minnesota Prior](../coefficient_draw/) the error
covariance matrix is replaced by an estimate, and the coefficients are drawn given it. Here
$`\boldsymbol{\Sigma}`$ is unknown as well, with an inverse-Wishart prior, and the sampler
alternates between two blocks: the coefficients given $`\boldsymbol{\Sigma}`$, drawn as in that
guide, and $`\boldsymbol{\Sigma}`$ given the coefficients. We run it on the VAR(2) of the earlier
guides, check the chain with the diagnostics of `bvar.diag`, and verify the sampler under the
natural conjugate prior, whose posterior is known in closed form.

To print the numbers and draw the figure first, run from the root of the repository:

```matlab
run guides/gibbs_sampler/guide.m
```

## The Model and Its Priors

The VAR is $`\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}`$ with
$`\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_T \otimes \boldsymbol{\Sigma})`$,
in the stacked form of [How to Stack the Equations of a VAR](../var_stacking/), or
$`\mathbf{Y} = \mathbf{Z}\mathbf{A} + \mathbf{U}`$ in the matrix form, with
$`\boldsymbol{\beta} = \mathrm{vec}(\mathbf{A})`$. The independent normal and inverse-Wishart prior
is

$$\boldsymbol{\beta} \sim \mathcal{N}(\boldsymbol{\beta}_0, \mathbf{V}), \qquad \boldsymbol{\Sigma} \sim \mathcal{IW}(\nu_0, \mathbf{S}_0),$$

with $`\boldsymbol{\beta}`$ and $`\boldsymbol{\Sigma}`$ independent a priori. The inverse-Wishart
density is proportional to
$`|\boldsymbol{\Sigma}|^{-(\nu_0 + n + 1)/2}\exp\left(-\frac{1}{2}\mathrm{tr}(\mathbf{S}_0\boldsymbol{\Sigma}^{-1})\right)`$,
and its mean is $`\mathbf{S}_0/(\nu_0 - n - 1)`$. We take $`\boldsymbol{\beta}_0 = \mathbf{0}`$ and
the Minnesota variances $`\mathbf{V}`$ of the second guide, with $`c_1 = 0.2^2`$, $`c_2 = 0.1^2`$
and $`c_3 = 10^2`$, and we set $`\nu_0 = n + 2`$ and
$`\mathbf{S}_0 = \mathrm{diag}(s_1^2, \ldots, s_n^2)`$, so that the prior mean of
$`\boldsymbol{\Sigma}`$ is the estimate that the second guide uses in its place. The data are
those of the earlier guides: three variables, two lags and $`T = 200`$.

## Step 1: Split the Unknowns into Blocks

Given $`\boldsymbol{\Sigma}`$, the VAR is a linear regression with a known error covariance
matrix, and given $`\boldsymbol{\beta}`$, the errors are known. The sampler therefore cycles
through two blocks:

1. $`(\boldsymbol{\beta} \mid \mathbf{y}, \boldsymbol{\Sigma})`$, a Gaussian distribution, from
   the second guide;
2. $`(\boldsymbol{\Sigma} \mid \mathbf{y}, \boldsymbol{\beta})`$, an inverse-Wishart
   distribution.

## Step 2: Derive the Conditional Distributions

**The coefficients.** Given $`\boldsymbol{\Sigma}`$, the posterior of $`\boldsymbol{\beta}`$ is
that of the second guide, with precision
$`\mathbf{K} = \mathbf{V}^{-1} + \boldsymbol{\Sigma}^{-1} \otimes \mathbf{Z}'\mathbf{Z}`$ evaluated
at the current $`\boldsymbol{\Sigma}`$, and
[`bvar.samplers.var_coef_joint`](../../core/+bvar/+samplers/var_coef_joint.m) takes one draw.

**The covariance matrix.** Given $`\boldsymbol{\beta}`$, the residuals
$`\mathbf{U} = \mathbf{Y} - \mathbf{Z}\mathbf{A}`$ are known. The likelihood in the trace form of
the first guide, times the prior density, gives

$$p(\boldsymbol{\Sigma} \mid \mathbf{y}, \boldsymbol{\beta}) \propto |\boldsymbol{\Sigma}|^{-\frac{T}{2}}\exp\left(-\frac{1}{2}\mathrm{tr}(\boldsymbol{\Sigma}^{-1}\mathbf{U}'\mathbf{U})\right) \times |\boldsymbol{\Sigma}|^{-\frac{\nu_0 + n + 1}{2}}\exp\left(-\frac{1}{2}\mathrm{tr}(\mathbf{S}_0\boldsymbol{\Sigma}^{-1})\right),$$

the kernel of an inverse-Wishart density:

$$(\boldsymbol{\Sigma} \mid \mathbf{y}, \boldsymbol{\beta}) \sim \mathcal{IW}(\nu_0 + T, \mathbf{S}_0 + \mathbf{U}'\mathbf{U}).$$

MATLAB's `iwishrnd(S, nu)` draws from $`\mathcal{IW}(\nu, \mathbf{S})`$ in this parameterization.

## Step 3: Write the Loop

```matlab
nu0 = n + 2; S0 = diag(s2);        % Sig ~ IW(nu0, S0), with prior mean S0
iV = spdiags(1./V, 0, k*n, k*n);   % beta ~ N(0, V), the Minnesota variances
nsim = 20000; burnin = 1000;
Sig = S0;                          % the starting value
store_A = zeros(k, n, nsim); store_Sig = zeros(n, n, nsim);
for isim = 1:burnin + nsim
    A = bvar.samplers.var_coef_joint(Y, Z, Sig, ones(T, 1), iV);   % beta given Sig
    U = Y - Z*A;
    Sig = iwishrnd(S0 + U'*U, nu0 + T);                            % Sig given beta
    if isim > burnin
        store_A(:, :, isim - burnin) = A;
        store_Sig(:, :, isim - burnin) = Sig;
    end
end
```

The chain starts from $`\boldsymbol{\Sigma} = \mathbf{S}_0`$, and the first 1,000 iterations are
discarded as burn-in. `var_coef_joint` forms $`\mathbf{Z}'\mathbf{Z}`$ at every call, since in
models with stochastic volatility its weights change from one call to the next. With constant
weights, as here, a sampler for a large VAR can form $`\mathbf{Z}'\mathbf{Z}`$ and
$`\mathbf{Z}'\mathbf{Y}`$ once, before the loop.

## Step 4: Check the Chain

In Figure 1, the draws fluctuate around a stable level throughout, with no drift. The inefficiency
factor of a parameter is approximately the number of draws that carry the information of one
independent draw. Computed by
[`bvar.diag.inefficiency_factor`](../../core/+bvar/+diag/inefficiency_factor.m) with Bartlett
weights up to lag 100, the inefficiency factors lie between 0.94 and 1.14 for the 21 coefficients
and between 1.00 and 1.11 for the six distinct elements of $`\boldsymbol{\Sigma}`$, so each
posterior mean is estimated almost as precisely as from independent draws. The Monte Carlo
standard error of a posterior mean, which [`bvar.diag.mcse`](../../core/+bvar/+diag/mcse.m)
computes, is then close to the posterior standard deviation over $`\sqrt{20{,}000}`$.

The diagnostic of Geweke (1992), [`bvar.diag.geweke`](../../core/+bvar/+diag/geweke.m), compares
the mean of the first 10% of the chain with the mean of its second half; if the chain is
stationary, the statistic is approximately standard normal. It exceeds 1.96 in absolute value for
3 of the 27 parameters, against 1.35 expected by chance. Its largest value, 3.29, for the
coefficient on the first lag of variable 3 in equation 1, is large for 27 tests, and ten runs of
the sampler from other seeds give largest values between 1.66 and 2.36. With inefficiency factors
near one, the value of 3.29 is consistent with a chance event of this run.

| | True value | Posterior mean | 90% interval | Inefficiency factor | MCSE |
|---|---|---|---|---|---|
| Own first lag, equation 1 | 0.50 | 0.530 | (0.425, 0.635) | 0.99 | 0.0004 |
| First lag of variable 1, equation 2 | 0.20 | 0.134 | (0.052, 0.215) | 1.10 | 0.0004 |
| $`\Sigma_{11}`$ | 1.00 | 1.008 | (0.851, 1.193) | 1.11 | 0.0008 |
| $`\Sigma_{21}`$ | 0.30 | 0.336 | (0.226, 0.457) | 1.06 | 0.0005 |

The 90% intervals contain the true values for 25 of the 27 parameters.

![The draws of four parameters and their histograms](fig_chain.png)

*Figure 1: Top, the 20,000 draws of two coefficients and two elements of
$`\boldsymbol{\Sigma}`$ kept after the burn-in, in the order the sampler makes them. Bottom, their
histograms and the true values (dashed lines).*

## How We Check It

<details>
<summary>The same loop under the natural conjugate prior, whose posterior is known in closed form</summary>

Under the natural conjugate prior,
$`(\boldsymbol{\beta} \mid \boldsymbol{\Sigma}) \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma} \otimes \mathbf{V}_{\mathbf{A}})`$
and $`\boldsymbol{\Sigma} \sim \mathcal{IW}(\nu_0, \mathbf{S}_0)`$, where the diagonal matrix
$`\mathbf{V}_{\mathbf{A}}`$ holds 100 for the intercept and $`0.2^2/(\ell^2 s_r^2)`$ for lag
$`\ell`$ of variable $`r`$. Two lines of the loop change: the prior precision of
$`\boldsymbol{\beta}`$ becomes $`\boldsymbol{\Sigma}^{-1} \otimes \mathbf{V}_{\mathbf{A}}^{-1}`$,
which depends on $`\boldsymbol{\Sigma}`$, and the prior of the coefficients enters the draw of
$`\boldsymbol{\Sigma}`$,

<div align="center">

$`(\boldsymbol{\Sigma} \mid \mathbf{y}, \boldsymbol{\beta}) \sim \mathcal{IW}(\nu_0 + T + k, \mathbf{S}_0 + \mathbf{U}'\mathbf{U} + \mathbf{A}'\mathbf{V}_{\mathbf{A}}^{-1}\mathbf{A}).`$

</div>

The posterior is normal-inverse-Wishart in closed form, with
$`E(\mathbf{A} \mid \mathbf{y}) = \hat{\mathbf{A}} = \mathbf{K}_{\mathbf{A}}^{-1}\mathbf{Z}'\mathbf{Y}`$,
$`\mathbf{K}_{\mathbf{A}} = \mathbf{V}_{\mathbf{A}}^{-1} + \mathbf{Z}'\mathbf{Z}`$, and
$`E(\boldsymbol{\Sigma} \mid \mathbf{y}) = \hat{\mathbf{S}}/(\nu_0 + T - n - 1)`$,
$`\hat{\mathbf{S}} = \mathbf{S}_0 + \mathbf{Y}'\mathbf{Y} - \hat{\mathbf{A}}'\mathbf{K}_{\mathbf{A}}\hat{\mathbf{A}}`$.
Over 20,000 draws, the Gibbs means of the 27 parameters differ from these by at most 2.68 Monte
Carlo standard errors, and the standard deviations of the coefficient draws differ from the exact
posterior standard deviations by at most 1.3%. The check verifies the two draws of the loop,
`var_coef_joint` and the inverse-Wishart step, under a prior whose posterior is known; the
posterior under the independent prior of Step 3 has no closed form to compare with.

</details>

## Adapt It Yourself

The draws give the posterior of the parameters. How do you obtain the predictive distribution of
the next observation, $`\mathbf{y}_{T+1}`$?

<details>
<summary>Answer</summary>

For each kept draw of $`(\mathbf{A}, \boldsymbol{\Sigma})`$, draw
$`\mathbf{y}_{T+1}' = \mathbf{x}_{T+1}'\mathbf{A} + \boldsymbol{\varepsilon}_{T+1}'`$ with
$`\boldsymbol{\varepsilon}_{T+1} \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Sigma})`$, where
$`\mathbf{x}_{T+1}' = (1, \mathbf{y}_T', \mathbf{y}_{T-1}')`$. The draws of
$`\mathbf{y}_{T+1}`$ are a sample from its predictive distribution, which accounts for the
uncertainty about both $`\mathbf{A}`$ and $`\boldsymbol{\Sigma}`$.

```matlab
xT1 = [1, Y(end, :), Y(end-1, :)];   % x_{T+1}' = (1, y_T', y_{T-1}')
yT1 = zeros(nsim, n);
for r = 1:nsim
    yT1(r, :) = xT1*store_A(:, :, r) + (chol(store_Sig(:, :, r), 'lower')*randn(n, 1))';
end
```

The predictive means are 4.03, 2.67 and 3.23, with 90% intervals (2.34, 5.70), (1.15, 4.17) and
(1.93, 4.56). `guide.m` checks the mean against the average of
$`\mathbf{x}_{T+1}'\mathbf{A}`$ over the draws, to within 0.011, and the variance against the
average error variance plus the variance of $`\mathbf{x}_{T+1}'\mathbf{A}`$, to within 1.6%.

</details>

## Next

In [How to Choose the Shrinkage of a VAR by the Marginal Likelihood](../natural_conjugate_ml/) the
shrinkage of the natural conjugate prior of How We Check It is chosen from the data. Among the
examples, [ex13](../../examples/ex13_mcmc_diagnostics.m) applies the same diagnostics to two
samplers with stochastic volatility, and [ex05](../../examples/ex05_var_csv.m) lets one volatility
path scale the whole covariance matrix of the VAR.

## Using This Method in Your Research

The independent normal and inverse-Wishart prior and its Gibbs sampler are standard, and
Chan (2027) gives a textbook treatment. The diagnostics of Step 4 are due to Geweke (1992).
[`CITING.md`](../../CITING.md) lists the paper behind each function of the toolkit and how to cite
the toolkit itself.

## References

Chan, J. C. C. (2027). *Bayesian Macroeconometrics: Methods and Applications*. Chapman &
Hall/CRC. [Code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics)

Geweke, J. (1992). Evaluating the Accuracy of Sampling-Based Approaches to the Calculation of
Posterior Moments. In J. M. Bernardo, J. O. Berger, A. P. Dawid and A. F. M. Smith (Eds.),
*Bayesian Statistics 4*, 169-193. Oxford University Press.
[doi:10.1093/oso/9780198522669.003.0010](https://doi.org/10.1093/oso/9780198522669.003.0010)
