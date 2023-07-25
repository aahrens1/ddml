{smcl}
{* *! version 25jul2023}{...}
{viewerjumpto "Stacking" "ddml_stacking##stacking"}{...}
{viewerjumpto "Standard stacking with pystacked" "ddml_stacking##std_stack"}{...}
{viewerjumpto "Pooled stacking" "ddml_stacking##pool_stack"}{...}
{viewerjumpto "Short-stacking" "ddml_stacking##short_stack"}{...}
{viewerjumpto "Re-stacking" "ddml_stacking##restack"}{...}
{viewerjumpto "Retrieving stacking weights" "ddml_stacking##stack_weights"}{...}
{viewerjumpto "Examples" "ddml_stacking##examples"}{...}
{viewerjumpto "Installation" "ddml_stacking##installation"}{...}
{viewerjumpto "References" "ddml_stacking##references"}{...}
{viewerjumpto "Authors" "ddml_stacking##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_stacking##also_see"}{...}
{hline}
{cmd:help ddml stacking}{right: v1.2}
{hline}

{title:ddml - Stata package for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}Stacking regression is a simple and powerful method for 
combining predictions from multiple learners.
{help pystacked} is the recommended way to specify multiple learners in {opt ddml},
and {opt ddml} has integrated support for various features provided by {help pystacked}.
This help file provides an overview of how to implement stacking
when estimating using {opt ddml}.


{marker stacking}{...}
{title:Stacking}

{pstd}
Stacking regression ({help ddml stacking##Wolpert1992:Wolpert, 1992})
is a way of combining predictions from multiple base ("level-0") learners
into a final prediction.
A final estimator ("level-1") is used to combine the base predictions.
A common approach is to use the cross-validated (out-of-sample, OOS) predictions
of the base learners to obtain the weights for combining the learners.

{pstd}
Three ways of pairing stacking with DDML are supported:
{it:standard stacking}, provided via the {help pystacked} package;
{it:pooled stacking}, a variant of standard stacking;
and {it:short-stacking}, a version of stacking specific to double-debiased machine learning.

{marker pystacked}{...}
{pstd}{help pystacked} is the recommended way to specify multiple learners in {opt ddml}.
{help pystacked} provides a fast way of estimating all the learners in a single call to one program,
and {opt ddml} has integrated support for various features provided by {help pystacked}.
{opt ddml} will store the predicted values of the specified base learners as well as the combined ("stacked") predicted values.
It also stores the standard and pooled stacking weights used by {help pystacked}
along with the {opt ddml} short-stacking weights.

{pstd}{bf:Important}: For these features to be available, {help pystacked} needs to be the only learner for each conditional expectation.
Multiple learners must be specified in the call to {help pystacked}; see the examples below.
{help pystacked} can be provided directly to {opt ddml} as one of several learners for a conditional expectation,
but in this case the extra features for {help pystacked} will not be availabe.

{pstd}Note: some of the {opt ddml} stacking options available via {help pystacked} integration
are not available for the flexible IV model.
See this {help ddml_example_flexiv_anylearner_detailed:help file} for examples and discussion
of how to stack and short-stack when using the flexible IV model).{p_end}


{marker std_stack}{...}
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


{marker pool_stack}{...}
{title:Pooled stacking}

{pstd}
Pooled stacking is a variant of standard stacking that implements additional regularization
via the {help ddml crossfit:cross-fitting} step of {opt ddml}.
Pooled stacking is done once, after all cross-fitting has been done
and a full set of all cross-validated OOS predictions has been obtained.
This means that a single set of stacking weights is used to obtain all OOS cross-fit predictions.
This is in contrast to standard stacking and k-fold cross-fitting,
where k different sets of stacking weights are estimated and used to obtain the k OOS cross-fit predictions.
Pooled stacking is specified at the {help ddml crossfit} stage using the {opt poolstack} option.
Pooled stacking is available only in conjunction with standard stacking.
The default final estimator is the same as with {help pystacked},
and can be changed using the {opt psfinalest(estimator)} option.

{pstd}Note: all final estimators available with {help pystacked} are also available for pooled stacking.
However, the current version of {help pystacked} generates the necessary cross-validated OOS predicted values
only if the standard stacking final estimator used by {help pystacked}
is either {opt nnls1} (the default), {opt ls1}, {opt ols}, {opt ridge} or {opt singlebest}.
Hence when using pooled stacking,
the standard stacking final estimator specified with {help pystacked} needs to be one of these.


{marker short_stack}{...}
{title:Short-stacking}

{pstd}
Short-stacking is a form of stacking specific to double debiased machine learning and cross-fitting.
Short-stacking uses the cross-fitted predicted values to obtain
the stacked (weighted average) of the multiple base learner predictions.
It is computationally faster (often much faster) than
either standard stacking or pooled stacking available via {help pystacked}.
Short-stacking also does not require use of {help pystacked};
the predictions of any base learners specified by the user can be short-stacked.
Short-stacking is specified at the {help ddml crossfit} stage using the {opt shortstack} option.
The default final estimator is the same as with {help pystacked},
and can be changed using the {opt finalest(estimator)} option.

{pstd}
Because short-stacking is typically much faster than standard or pooled stacking,
users may wish to use short-stacking as the only stacking method.
This can be done efficiently in combination with {help pystacked}.
To do this, (1) use {help pystacked} as the single learner in each equation;
(2) at the cross-fitting stage, specify the {opt shortstack} and {cmdab:nostd:stack} options.
This causes {help pystacked} to estimate the base learners
without the computationally-costly stacking step in each cross-fit fold.


{marker restack}{...}
{title:Re-stacking after cross-fitting}

{pstd}
Users have the option of re-stacking the base learner predictions using a different final estimator
without having to re-cross-fit/re-estimate the entire model.
This is done by specifying the stacking method and final estimator
at the {help ddml estimate:ddml estimate} step.
This feature is available only if {help pystacked} is the single learner in every equation.


{marker stack_weights}{...}
{title:Stacking weights}

{pstd}
The weights used for standard stacking, pooled stacking and short-stacking
can be inspected after estimation using {help ddml extract}
with the {opt show(stweights)}, {opt show(psweights)} and/or {opt show(ssweights)}, respectively.
In the case of standard stacking, the mean weights across cross-fit folds are displayed;
to display the standard stacking weights for all k folds along with the separate learner MSEs,
use the {opt show(pystacked)} option.


{marker examples}{...}
{title:Examples}

{pstd}
See {help ddml init:help ddml init} for details of model initialization and learner specification options.

{pstd}Note: the additional support provided by {opt ddml} for {help pystacked} (see above)
is available only if {help pystacked} is the sole learner for each conditional expectation.
Mutliple learners are provided to {help pystacked}, not directly to {opt ddml}.{p_end}


{smcl}
INCLUDE help ddml_example_stacking.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
