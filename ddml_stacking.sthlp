{smcl}
{* *! version 26jun2023}{...}
{hline}
{cmd:help ddml stacking}{right: v1.2}
{hline}

{title:ddml - Stata package for Double Debiased Machine Learning}

{p2colset 5 19 21 2}{...}
{p2col:{hi: ddml} {hline 2}}Stata package for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}Stacking regression is a simple and powerful method for 
combining predictions from multiple learners.
{helpb pystacked} is the recommended way to specify multiple learners in {opt ddml},
and {opt ddml} has integrated support for various features provided by {helpb pystacked}.
This help file provides an overview of how to implement stacking
when estimating using {opt ddml}.


{marker section_stacking}{...}
{title:Stacking}

{pstd}
Stacking regression ({help ddml stacking##Wolpert1992:Wolpert, 1992})
is a way of combining predictions from multiple base ("level-0") learners
into a final prediction.
A final estimator ("level-1") is used to combine the base predictions.
A common approach is to use the cross-validated (out-of-sample, OOS) predictions
of the base learners to obtain the weights for combining the learners.

{pstd}
Three versions of stacking are supported by {opt ddml}:
standard stacking, provided via the {help pystacked} package;
pooled stacking, a variant of standard stacking;
and short-stacking, a version of stacking specific to double-debiased machine learning.

{pstd}{helpb pystacked} is the recommended way to specify multiple learners in {opt ddml}.
{helpb pystacked} provides a fast way of estimating all the learners in a single call to one program,
and {opt ddml} has integrated support for various features provided by {helpb pystacked}.
{opt ddml} will store the predicted values of the specified base learners as well as the combined ("stacked") predicted values.
It also stores the standard and pooled stacking weights used by {help pystacked}
along with the {opt ddml} short-stacking weights.

{pstd}{bf:Important}: For these features to be available, {helpb pystacked} needs to be the only learner for each conditional expectation.
Multiple learners must be specified in the call to {helpb pystacked}; see the examples below.
{helpb pystacked} can be provided directly to {opt ddml} as one of several learners for a conditional expectation,
but in this case the extra features for {helpb pystacked} will not be availabe.


{title:Standard stacking}

{pstd}
Standard stacking is implemented via the {help pystacked} package.
This is done by specifying {help pystacked} as the learner for a conditional expectation;
{help pystacked} in turn estimates using the user-specified base learners,
and stacks them to get the stacked (ensemble) prediction.
This is done in the context of {opt ddml}'s cross-fitting algorithm,
meaning that for k-fold cross-fitting, stacking is done k times,
once for each of the cross-fit folds.

{pstd}
The {help pystacked} base learners are specified at the {help ddml eq} stage,
when the supervised ML learners are added.
This is also where other {help pystacked} options can be specified, e.g.,
the final ("level-1") estimator used to combined the predictions.
The {help pystacked} default final predictor for stacking
regession is non-negative least squares (NNLS) without an intercept
and with the constraint that weights sum to one.
See {help pystacked} for alternative final estimators.


{title:Pooled stacking}

{pstd}
Pooled stacking is a variant of standard stacking that implements additional regularization
via the {help ddml crossfit:cross-fitting} step of {opt ddml}.
Pooled stacking is done once, after all cross-fitting has been done
and a full set of all cross-validated OOS predictions has been obtained.
This means that a single set of stacking weight is used to obtain all OOS cross-fit predictions.
This is in contrast to standard stacking and k-fold cross-fitting,
where k different sets of stacking weights are estimated and used to obtain the k OOS cross-fit predictions.
Pooled stacking is specified at the {help ddml crossfit} stage using the {opt poolstack} option.
The default final estimator is the same as with {help pystacked},
and can be changed using the {opt psfinalest(estimator)} option.


{title:Short-stacking}

{pstd}
Short-stacking is a form of stacking specific to double-debiased machine learning and cross-fitting.
Short-stacking uses the cross-fitted predicted values to obtain
the stacked (weighted average) of the multiple base learner predictions.
It is computationally faster (often much faster) than
either standard stacking or pooled stacking available via {help pystacked}.
Short-stacking also does not require use of {help pystacked};
the predictions of any base learners specified by the user can be short-stacked.
Short-stacking is specified at the {help ddml crossfit} stage using the {opt shortstack} option.
The default final estimator is the same as with {help pystacked},
and can be changed using the {opt ssfinalest(estimator)} option.


{title:Stacking weights}

{pstd}
The weights used for standard stacking, pooled stacking and short-stacking
can be inspected after estimation using {help ddml extract}
with the {opt show(stweights)}, {opt show(psweights)} and/or {opt show(ssweights)}, respectively.
In the case of standard stacking, the mean weights across cross-fit folds are displayed;
to display all k standard stacking weights along with the separate learner MSEs,
use the {opt show(pystacked)} option.


{marker examples}{...}
{title:Examples}

For more examples of usage see {help ddml##examples:help ddml}.
See {help ddml init:help ddml init} for details of model initialization and learner specification options.

{pstd}Note: the additional support provided by {opt ddml} for {helpb pystacked} (see {help ddml##pystacked:above})
is available only if, as in this example, {help pystacked} is the only learner for each conditional expectation.
Mutliple learners are provided to {help pystacked}, not directly to {opt ddml}.{p_end}

{pstd}Preparation: load the data, define global macros, set the seed and initialize the model.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}

{pstd}Add supervised machine learners for estimating conditional expectations.
For simplicity, we use {help pystacked}'s default learners: OLS, cross-validated lasso, and gradient boosting.{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}

{pstd} Cross-fitting: The learners are iteratively fitted on the training data.
In addition to the standard stacking done by {helpb pystacked},
also request short-stacking and pooled-stacking to be done by {opt ddml}.{p_end}
{phang2}. {stata "ddml crossfit, shortstack poolstack"}{p_end}

{pstd}Estimate the coefficients of interest.
Specify heteroskedastic-consistent SEs.{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Examine the standard ({cmd:pystacked}) stacking weights as well as
the {opt ddml} short-stacking and pooled-stacking weights.{p_end}
{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}
{phang2}. {stata "ddml extract, show(psweights)"}{p_end}

{pstd}Shorthand for displaying all weights:{p_end}
{phang2}. {stata "ddml extract, show(weights)"}{p_end}

{pstd}Examine the full set of {help pystacked} stacking weights
by cross-fit folds plus other {help pystacked} results:{p_end}
{phang2}. {stata "ddml extract, show(weights)"}{p_end}

{pstd}As above, but specify the base learners explicitly.
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


{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

{marker Hastie2009}{...}
{pstd}
Hastie, T., Tibshirani, R., & Friedman, J. (2009). 
The elements of statistical learning: data mining, inference,
and prediction. Springer Science & Business Media.

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}


{title:Authors}

{pstd}
Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland  {break}
achim.ahrens@gess.ethz.ch

{pstd}
Christian B. Hansen, University of Chicago, USA {break}
Christian.Hansen@chicagobooth.edu

{pstd}
Mark E Schaffer, Heriot-Watt University, UK {break}
m.e.schaffer@hw.ac.uk	

{pstd}
Thomas Wiemann, University of Chicago, USA {break}
wiemann@uchicago.edu


{title:Also see (if installed)}

{pstd}
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}, {helpb pystacked}.{p_end}
