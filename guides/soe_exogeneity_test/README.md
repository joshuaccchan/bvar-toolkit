# How to Test Block Exogeneity by the Marginal Likelihood

*Code: [`guide.m`](guide.m), with [`bvar.models.var_soe`](../../core/+bvar/+models/var_soe.m),
[`bvar.ml.ncp`](../../core/+bvar/+ml/ncp.m) and
[`bvar.ml.lniwpdf`](../../core/+bvar/+ml/lniwpdf.m). Method:
[Cushman and Zha (1997)](../../CITING.md#cushman-and-zha-1997),
[Zha (1999)](../../CITING.md#zha-1999) and [Chan (2027)](../../CITING.md#chan-2027). Data:
[`soe_Q.csv`](../data/).*

[How to Estimate a VAR for a Small Open Economy](../soe_var/) treats Australia as the small open
economy and the US as the foreign economy. Its VAR has a foreign block of US variables and a
domestic block of Australian variables. That guide imposes block exogeneity: the lags of the
Australian variables are excluded from the US equations. To assess the restriction, we compare
that model with a VAR in which the Australian lags enter the US equations. The Bayes factor in
favor of the restriction is the ratio of the restricted model's marginal likelihood to the
unrestricted model's.

We compute this ratio in closed form. The prior hyperparameters are the same in both models. The
Australian block is then identical in the two models, so its marginal likelihood cancels in the
ratio. Only the two marginal likelihoods of the US block are needed. For the simulated sample of
the estimation guide, the log Bayes factor is 17.66 in favor of block exogeneity.

The Bayes factor depends on the prior variances of the Australian lags. In that sample, the log
Bayes factor is 1.41 when these variances are a hundredth of those in the estimation guide. It is
160.02 when they are a thousand times as large. With large enough prior variances, the test
favors block exogeneity in any sample. This is the Bartlett-Lindley paradox (Lindley, 1957;
Bartlett, 1957). The prior variances should therefore be defensible and not too extreme.

To reproduce the numbers, run the following command from the root of the repository:

```matlab
run guides/soe_exogeneity_test/guide.m
```

## The Two Models

The VAR with block exogeneity has two blocks. The US block is a VAR of the US variables
$`\mathbf{y}^f_t`$. The Australian block regresses the Australian variables $`\mathbf{y}^d_t`$ on
all lags and on the current US variables:

```math
\mathbf{y}^f_t = \mathbf{a}^f + \sum_{j=1}^{p}\mathbf{A}_j\mathbf{y}^f_{t-j} + \mathbf{u}^f_t, \qquad \mathbf{y}^d_t = \mathbf{c} + \mathbf{G}_0\mathbf{y}^f_t + \sum_{j=1}^{p}\left(\mathbf{G}_j\mathbf{y}^f_{t-j} + \mathbf{D}_j\mathbf{y}^d_{t-j}\right) + \mathbf{e}_t.
```

The unrestricted VAR adds the Australian lags to the US block:

```math
\mathbf{y}^f_t = \mathbf{a}^f + \sum_{j=1}^{p}\left(\mathbf{A}_j\mathbf{y}^f_{t-j} + \mathbf{C}_j\mathbf{y}^d_{t-j}\right) + \mathbf{u}^f_t.
```

Its Australian block is the same regression. Block exogeneity is the restriction
$`\mathbf{C}_1 = \cdots = \mathbf{C}_p = \mathbf{0}`$. With three lags it sets 12 coefficients to
zero in each of the four US equations, 48 in all.

Both models get the prior of the estimation guide. Each block has a natural conjugate prior. For
the US block with coefficient matrix $`\mathbf{B}`$,

```math
(\mathrm{vec}(\mathbf{B}) \mid \boldsymbol{\Sigma}^f) \sim \mathcal{N}(\mathrm{vec}(\mathbf{B}_0), \boldsymbol{\Sigma}^f \otimes \mathbf{V}), \qquad \boldsymbol{\Sigma}^f \sim \mathcal{IW}(6, \mathbf{S}_0).
```

The prior mean $`\mathbf{B}_0`$ is one for the own first lag of the federal funds rate and of the
excess bond premium, and zero elsewhere. The diagonal matrix $`\mathbf{V}`$ contains
$`c_1/(\ell^2 s_r^2)`$ for lag $`\ell`$ of variable $`r`$. As in the estimation guide, $`c_1`$
maximizes the marginal likelihood of the VAR with block exogeneity, which gives 0.070 on the
sample. In the unrestricted VAR, the lags of the Australian variables get a prior mean of zero and
variances of the same form.

## Step 1: Write the Bayes Factor with the US Block

Let $`\mathbf{x}_t`$ contain the intercept and the $`p`$ lags of all eight variables. Let
$`\boldsymbol{\theta}^f = (\mathbf{B}, \boldsymbol{\Sigma}^f)`$ denote the parameters of the US
block and $`\boldsymbol{\theta}^d`$ those of the Australian block. In each model $`M`$, the
density of $`\mathbf{y}_t`$ given the past is the product of the densities of the two blocks:

```math
p(\mathbf{y}_t \mid \mathbf{x}_t, \boldsymbol{\theta}^f, \boldsymbol{\theta}^d, M) = p(\mathbf{y}^f_t \mid \mathbf{x}_t, \boldsymbol{\theta}^f, M)\, p(\mathbf{y}^d_t \mid \mathbf{y}^f_t, \mathbf{x}_t, \boldsymbol{\theta}^d).
```

The likelihood of the sample $`\mathbf{Y}`$, given the quarters before it, is the product of these
densities over $`t = 1, \ldots, T`$. The two blocks have independent priors. The marginal
likelihood is therefore the product of one integral for each block:

```math
p(\mathbf{Y} \mid M) = m^f(M)\, m^d, \qquad m^f(M) = \int \prod_{t=1}^{T} p(\mathbf{y}^f_t \mid \mathbf{x}_t, \boldsymbol{\theta}^f, M)\, p(\boldsymbol{\theta}^f \mid M)\, \mathrm{d}\boldsymbol{\theta}^f.
```

We call $`m^f(M)`$ the marginal likelihood of the US block. In the unrestricted VAR, the
Australian lags in $`\mathbf{x}_t`$ enter the US block. Each density of that block then conditions
on the past of the Australian variables. The factor $`m^d`$ is the same integral for the
Australian block. This block and its prior are the same in both models, so $`m^d`$ cancels in the
ratio. With $`M_r`$ the VAR with block exogeneity and $`M_u`$ the unrestricted VAR, the log Bayes
factor is

```math
\log \mathrm{BF} = \log m^f(M_r) - \log m^f(M_u).
```

A positive value favors block exogeneity. With equal prior probabilities of the two models, the
posterior probability of block exogeneity is $`1/(1 + \mathrm{e}^{-\log \mathrm{BF}})`$.

The Bayes factor therefore measures one thing: whether the Australian lags help to predict the US
variables. The marginal likelihood of a model is the product of its one-step predictive densities
over the sample. The unrestricted VAR predicts better only if the Australian lags carry
information about the US variables. Otherwise its 48 additional coefficients add noise to its
predictions.

## Step 2: Compute the Two Marginal Likelihoods

Each marginal likelihood has the closed form of
[How to Choose the Shrinkage of a VAR by the Marginal Likelihood](../natural_conjugate_ml/), which
the function [`bvar.ml.ncp`](../../core/+bvar/+ml/ncp.m) computes. The two calls differ only in
the regressors and in their prior means and variances:

```matlab
VA = res.VA; B0 = res.B0;                 % the prior variances and means on x_t
nu0 = nus + 2; S0 = diag(s2(us));         % Sig_us ~ IW(nu0, S0) in both models
lml_r = bvar.ml.ncp(Ys(:, us), Zs(:, ius), VA(ius), nu0, S0, B0(ius, us));   % block exogeneity
lml_u = bvar.ml.ncp(Ys(:, us), Zs, VA, nu0, S0, B0(:, us));     % with the Australian lags
lbf = lml_r - lml_u;                      % the log Bayes factor for block exogeneity
```

Here `res` is the output of [`bvar.models.var_soe`](../../core/+bvar/+models/var_soe.m) with
`'nsim', 0`, which sets the prior and chooses $`c_1`$ as in the estimation guide. As in that
guide, `Zs` contains the intercept and all lags, and `ius` indexes the intercept and the US lags.
On the simulated sample of that guide, the log marginal likelihood of the US block is -688.77 with
block exogeneity and -706.43 without it. The log Bayes factor is 17.66.

The data also favor block exogeneity. The same computation on the 144 quarters of
[`soe_Q.csv`](../data/) gives a log Bayes factor of 9.62.

Both numbers are empirical Bayes results. The prior uses two quantities estimated from the sample
under test. The first is the set of AR(4) variances $`s_r^2`$. The second is the $`c_1`$ that
maximizes the marginal likelihood of the VAR with block exogeneity. Both are treated as fixed in
the marginal likelihoods.

## Step 3: Vary the Prior Variances of the Australian Lags

The Bayes factor depends on the prior of the 48 coefficients that the restriction sets to zero. We
multiply the prior variances of the Australian lags in the US block by factors from 0.01 to 1,000.
The table gives the log Bayes factor of the sample and of the data at each factor.

| Multiple of the prior variances | 0.01 | 0.1 | 1 | 10 | 100 | 1,000 |
|---|---|---|---|---|---|---|
| The sample: log Bayes factor | 1.41 | 4.91 | 17.66 | 53.51 | 105.16 | 160.02 |
| The data: log Bayes factor | 0.10 | 1.55 | 9.62 | 40.25 | 89.69 | 144.24 |

We note two patterns. First, the log Bayes factor is modest at small multiples. At a multiple of
0.01, the prior of the 48 coefficients concentrates near their prior mean of zero. The two models
are then nearly the same. The log Bayes factor is 1.41 in the sample, odds of about 4 to 1, and
0.10 on the data.

Second, the log Bayes factor rises with the multiple. With large prior variances, the prior places
most of its mass on coefficient values far from those consistent with the data. Each tenfold
increase then lowers the prior density of each of the 48 coefficients near those values by a
factor of about $`\sqrt{10}`$. The log Bayes factor rises by up to $`24 \log 10 = 55.26`$ with
each tenfold increase, whatever the data. From 100 to 1,000 it rises by 54.86 in the sample and
by 54.55 on the data. With large enough prior variances, the test thus favors block exogeneity
in any sample. This is the Bartlett-Lindley paradox.

A Bayes factor for a zero restriction therefore depends on the prior variances of the restricted
coefficients. These variances should be defensible and not too extreme. The prior of Step 2, at a
multiple of 1, shrinks the Australian lags exactly as it shrinks the US lags. A report of the test
should state the prior and show how the result changes with it.

## How We Check It

<details>
<summary>Bayes' theorem at the true values, and the unrestricted VAR of the estimation guide</summary>

By Bayes' theorem, at every point $`\boldsymbol{\theta}^f = (\mathbf{B}, \boldsymbol{\Sigma}^f)`$,

<div align="center">

$`\log m^f(M) = \sum_{t=1}^{T} \log p(\mathbf{y}^f_t \mid \mathbf{x}_t, \boldsymbol{\theta}^f, M) + \log p(\boldsymbol{\theta}^f \mid M) - \log p(\boldsymbol{\theta}^f \mid \mathbf{Y}, M),`$

</div>

where the prior and the posterior densities are normal-inverse-Wishart. `guide.m` evaluates the
right side at the true values with [`bvar.ml.lniwpdf`](../../core/+bvar/+ml/lniwpdf.m). It agrees
with `bvar.ml.ncp` to within 5e-12 for both models.

[How to Estimate a VAR for a Small Open Economy](../soe_var/) compares its estimates with an
unrestricted VAR. That VAR has one natural conjugate prior for all eight equations, with 10
degrees of freedom for the inverse-Wishart prior of its error covariance matrix. The posterior of
its US equations is the posterior of the US block with the Australian lags of this guide. The two
posterior means agree exactly, and the two scale matrices agree to within 2e-13.

</details>

## Adapt It Yourself

The restriction can also be tested for one Australian variable at a time. Which lines of `guide.m`
change to test whether the lags of the Australian dollar enter the US equations?

<details>
<summary>Answer</summary>

Only the regressors of the restricted US block. It keeps every column of `Zs` except the lags of
the variable under test:

```matlab
lbf_j = zeros(1, nau);
for j = 1:nau
    keep = setdiff(1:k, 1 + au(j) + (0:p-1)*n);         % x_t without the lags of variable j
    lbf_j(j) = bvar.ml.ncp(Ys(:, us), Zs(:, keep), VA(keep), nu0, S0, B0(keep, us)) - lml_u;
end
```

Each test sets 12 coefficients to zero. On the sample, the log Bayes factors are 0.92 for
Australian GDP growth and 2.80 for CPI inflation. They are 10.00 for the interest rate and 0.69
for the depreciation of the Australian dollar. All 48 true coefficients are zero. Each of these
tests is less decisive than the joint test, whose log Bayes factor is 17.66. The results for GDP
growth and the Australian dollar are below one, so the sample hardly discriminates between the two
models for those variables.

</details>

## Next

[How to Trace a Foreign Shock through a Small Open Economy](../soe_foreign_shock/) traces a shock
to the US excess bond premium through the Australian economy.

## Using This Method in Your Research

Cushman and Zha (1997) impose block exogeneity on a VAR for Canada and report a chi-square test
of the restrictions of their model, block exogeneity among them. The block form of the likelihood
follows Zha (1999). Chan (2027) gives a textbook treatment of the natural conjugate prior and its
marginal likelihood. [`CITING.md`](../../CITING.md) lists the paper behind each function of the
toolkit and how to cite the toolkit itself.

## References

Bartlett, M. S. (1957). A Comment on D. V. Lindley's Statistical Paradox. *Biometrika*, 44(3-4):
533-534. [doi:10.1093/biomet/44.3-4.533](https://doi.org/10.1093/biomet/44.3-4.533)

Chan, J. C. C. (2027). *Bayesian Macroeconometrics: Methods and Applications*. Chapman &
Hall/CRC. [Code repository](https://github.com/joshuaccchan/bayesian-macroeconometrics)

Cushman, D. O. and Zha, T. (1997). Identifying Monetary Policy in a Small Open Economy under
Flexible Exchange Rates. *Journal of Monetary Economics*, 39(3): 433-448.
[doi:10.1016/S0304-3932(97)00029-9](https://doi.org/10.1016/S0304-3932(97)00029-9)

Lindley, D. V. (1957). A Statistical Paradox. *Biometrika*, 44(1-2): 187-192.
[doi:10.1093/biomet/44.1-2.187](https://doi.org/10.1093/biomet/44.1-2.187)

Zha, T. (1999). Block Recursion and Structural Vector Autoregressions. *Journal of Econometrics*,
90(2): 291-316. [doi:10.1016/S0304-4076(98)00045-1](https://doi.org/10.1016/S0304-4076(98)00045-1)
