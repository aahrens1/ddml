{smcl}
{* *! version 17aug2023}{...}
{smcl}
{pstd}{ul:Partially-linear model with {help pystacked} and stacking}:{p_end}

{pstd}Preparation: load the data, define global macros, set the seed and initialize the model.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}

{pstd}Add supervised machine learners for estimating conditional expectations.
For simplicity, we use {help pystacked}'s default learners:
OLS, cross-validated lasso, and gradient boosting.{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}

{pstd} Cross-fitting and estimation:
The learners are iteratively fitted on the training data to obtain the estimated conditional expectations,
and then the causal coefficient of interest is estimated along with heteroskedastic-consistent SEs.
Note that the initial stacking is specified at the {help ddml crossfit:cross-fitting} stage.
In addition to the standard stacking done by {helpb pystacked},
also request short-stacking and pooled-stacking to be done by {opt ddml}.{p_end}

{phang2}. {stata "ddml crossfit, shortstack poolstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Examine the standard ({cmd:pystacked}) stacking weights as well as
the {opt ddml} short-stacking and pooled-stacking weights.{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}
{phang2}. {stata "ddml extract, show(psweights)"}{p_end}

{pstd} Re-stack without cross-fitting, using the single-best learner
instead of the default constrained nonlinear least squares.
We do this using the {help ddml estimate} command.
Since no stacking method is specified,
restacking will be done for all three methods.{p_end}

{phang2}. {stata "ddml estimate, robust finalest(singlebest)"}{p_end}

{pstd} As above, but request short-stacking only at the cross-fitting stage.
Note the speed improvement.{p_end}

{phang2}. {stata "ddml crossfit, shortstack nostdstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd} Re-stack the above without cross-fitting, using OLS as the final estimator.
Use the option {opt shortstack} since only these results are re-stacked.{p_end}

{phang2}. {stata "ddml estimate, robust shortstack finalest(ols)"}{p_end}
{phang2}. {stata "ddml estimate, robust shortstack finalest(ols)"}{p_end}

{pstd}{ul:Extended example with specified {help pystacked} learners and settings}:{p_end}

{pstd}Same example as above, but specify the base learners explicitly.
We again make use of {help pystacked} integration,
so there is a single call to {help pystacked} for each conditional expectation.
The first learner in the stacked ensemble is OLS.
We also use cross-validated lasso, ridge and two random forests with different settings.
The settings are stored in macros for readability.{p_end}

{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}
{phang2}. {stata "global rflow max_features(5) min_samples_leaf(1) max_samples(.7)"}{p_end}
{phang2}. {stata "global rfhigh max_features(5) min_samples_leaf(10) max_samples(.7)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X || method(ols) || method(lassocv) || method(ridgecv) || method(rf) opt($rflow) || method(rf) opt($rfhigh), type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X || method(ols) || method(lassocv) || method(ridgecv) || method(rf) opt($rflow) || method(rf) opt($rfhigh), type(reg)"}{p_end}

{pstd}Note: Options before ":" and after the first comma refer to {cmd:ddml}. 
Options that come after the final comma refer to the estimation command. 
Make sure to not confuse the two types of options.{p_end}

{pstd}The learners are iteratively fitted on the training data.
In addition to the standard stacking done by {helpb pystacked},
also request short-stacking to be done by {opt ddml}.
Finally, estimate the coefficients of interest.{p_end}

{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Examine the standard ({cmd:pystacked}) stacking weights as well as the {opt ddml} short-stacking weights.{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}
