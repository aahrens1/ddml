{smcl}
{* *! version 25july2023}{...}
{viewerjumpto "Syntax" "ddml_estimate##syntax"}{...}
{viewerjumpto "Cross-fit options" "ddml_estimate##crossfit"}{...}
{viewerjumpto "Estimation options" "ddml_estimate##estimation"}{...}
{viewerjumpto "Replay options" "ddml_estimate##replay"}{...}
{viewerjumpto "User-specified variables" "ddml_estimate##userspec"}{...}
{viewerjumpto "Installation" "ddml_estimate##installation"}{...}
{viewerjumpto "References" "ddml_estimate##references"}{...}
{viewerjumpto "Authors" "ddml_estimate##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_estimate##also_see"}{...}
{hline}
{cmd:help ddml crossfit, ddml estimate}{right: v1.2}
{hline}

{title:ddml crossfit and estimate commands for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml crossfit} implements the cross-fitting algorithm.
Each learner is fitted iteratively on training folds and out-of-sample predicted values are obtained.

{pstd}
{opt ddml estimate} estimates the model using the conditional expectations obtained from the cross-fitting step.

{pstd}
Before cross-fitting, the model must be defined using {help ddml init} and the learners specified using {help ddml eq}.

{pstd}
See the help for {help ddml stacking} for a detailed discussion and examples of stacking with {opt ddml}.


{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml crossfit} [ , {opt mname(name)} {opt shortstack} {opt poolstack} {cmdab:NOSTD:stack} {opt finalest(name)}{bind: ]} 

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)}
{opt atet} {opt ateu} {opt trim(real)}
{opt mname(name)} {opt shortstack} {opt poolstack} {opt stdstack} {opt finalest(name)}{bind: ]} 

{p 8 14}Replay options (available after model estimation):

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {opt spec(integer or string)} {opt rep(integer or string)} {opt allcombos} {opt not:able} {opt replay}{bind: ]} 

{p 8 14}Using a user-specified combination of {help ddml crossfit} conditional expectations (models {opt partial}, {opt partialiv}, {opt fiv}):

{p 8 14}{cmd:ddml estimate} , {opt y(varname)} {opt d(varlist)} [ {opt z(varlist)} {opt dh(varname)} {opt mname(name)} {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)}{bind: ]} 

{p 8 14}Using a user-specified combination of {help ddml crossfit} conditional expectations (models {opt interactive}, {opt interactiviv}):

{p 8 14}{cmd:ddml estimate} , {opt y0(varname)} {opt y1(varname)} [ {opt d(varname)} {opt d0(varname)} {opt d1(varname)} {opt z(varname)} {opt mname(name)}  {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)}{bind: ]}


{marker crossfit}{...}
{synoptset 20}{...}
{synopthdr:Cross-fitting}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt shortstack}} asks for short-stacking to be used.
Short-stacking uses the
cross-fitted predicted values to obtain a weighted average
of the multiple base learners.
It is computationally faster (often much faster) than standard stacking
(implemented via {help pystacked}) is used to specify the base learners
{p_end}
{synopt:{opt poolstack}} is available as an alternative to standard stacking
when {help pystacked} is used to specify the base learners.
Pooled-stacking adds additional regularization
by obtaining a single set of stacking weights
from the full set of out-of-sample base learner predicted values
(in contrast to {help pystacked}, which stacks each cross-fit fold separately).
{p_end}
{synopt:{opt nostdstack}} is used in conjunction with short-stacking and {help pystacked}.
It tells {help pystacked} to generate the base learner predictions without
the computationally-expensive additional step of obtaining the stacking weights.
This option should be used if short-stacking is the only stacking method needed.
{p_end}
{synopt:{opt finalest(name)}} sets the final estimator for all stacking methods;
the default is the {help pystacked} default of non-negative nonlinear least squares.
See {help pystacked} for alternative stacking final estimators.
NB: use of this option is incompatible with use of the {opt finalest(.)} option
when {help pystacked} is the learner specified in an equation using {help ddml eq};
use {opt finalest} in one or the other, or neither (the default), but not both.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{marker estimation}{...}
{synoptset 20}{...}
{synopthdr:Estimation}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{cmdab:r:obust}}
report SEs that are robust to the
presence of arbitrary heteroskedasticity.
{p_end}
{synopt:{opt cluster(varname)}}
select cluster-robust variance-covariance estimator, e.g. {cmd:vce(hc3)} or {cmd:vce(cluster id)}.
{p_end}
{synopt:{opt vce(type)}}
select variance-covariance estimator; see {help regress##vcetype:here}.
{p_end}
{synopt:{cmdab:noc:onstant}}
suppress constant term ({it:partial}, {it:iv}, {it:fiv} models only). Since the residualized outcome 
and treatment may not be exactly mean-zero in finite samples, {cmd:ddml} includes the constant by 
default in the estimation stage of partially linear models.
{p_end}
{synopt:{cmdab:showc:onstant}}
display constant term in summary estimation output table ({it:partial}, {it:iv}, {it:fiv} models only).
{p_end}
{synopt:{opt atet}}
report average treatment effect of the treated (default is ATE).
{p_end}
{synopt:{opt ateu}}
report average treatment effect of the untreated (default is ATE).
{p_end}
{synopt:{opt trim(real)}}
trimming of propensity scores for the Interactive and Interactive IV models. The default is 0.01
(that is, values below 0.01 and above 0.99 are set 
to 0.01 and 0.99, respectively).
{p_end}
{synopt:{opt shortstack}} requests re-stacking of the short-stacking results
using the final estimator specified with {opt finalest(.)};
this option is available only if {help pystacked} is the single learner for each equation.
Re-stacking is fast because it doesn't require re-cross-fitting.
{p_end}
{synopt:{opt poolstack}} requests re-stacking of the pooled stacking results
using the final estimator specified with {opt finalest(.)};
this option is available only if {help pystacked} is the single learner for each equation.
Re-stacking is fast because it doesn't require re-cross-fitting.
{p_end}
{synopt:{opt stdstack}} requests re-stacking of the standard stacking results
using the final estimator specified with {opt finalest(.)};
this option is available only if {help pystacked} is the single learner for each equation.
Re-stacking is fast because it doesn't require re-cross-fitting.
{p_end}
{synopt:{opt finalest(name)}} sets the final estimator for all stacking methods;
the default is the {help pystacked} default of non-negative nonlinear least squares.
See {help pystacked} for alternative stacking final estimators.
{p_end}
{p2colreset}{...}
{pstd}

{marker replay}{...}
{synoptset 20}{...}
{synopthdr:Replay}
{synoptline}
{synopt:{opt spec(integer/string)}}
select specification. This can either be the specification number, {it:mse} for minimum-MSE specification (the default) or {it:ss} for short-stacking. 
{p_end}
{synopt:{opt rep(integer/string)}}
select resampling iteration. This can either be the cross-fit repetition number, {it:mn} for mean aggregation or {it:md} for median aggregation (the default).
{p_end}
{synopt:{opt allcombos}}
estimates all possible specifications. By default, only the min-MSE, short-stacking or or pooled-stacking
specification is estimated and displayed.
{p_end}
{synopt:{opt replay}}
used in combination with {opt spec()} and {opt rep()} to display and return estimation results.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{marker userspec}{...}
{synoptset 20}{...}
{synopthdr:User-specified vars}
{synoptline}
{synopt:{opt y(varname)}}
estimated conditional expectation of dependent variable (models {opt partial}, {opt partialiv}, {opt fiv})
{p_end}
{synopt:{opt d(varname)}}
estimated conditional expectation of causal variable of interest (models {opt partial}, {opt partialiv}, {opt fiv})
{p_end}
{synopt:{opt z(varname)}}
estimated conditional expectation of instrumental variables (models {opt partialiv}, {opt interactiv})
{p_end}
{synopt:{opt dh(varname)}}
estimated optimal IV = E[D|X,Z] - E[D^|X] (model {opt fiv} only)
{p_end}
{synopt:{opt y0(varname)}}
estimated E[Y|X,D=0] (model {opt interactive}) or E[Y|X,Z=0] (model {opt interactiveiv})
{p_end}
{synopt:{opt y1(varname)}}
estimated E[Y|X,D=1] (model {opt interactive}) or E[Y|X,Z=1] (model {opt interactiveiv})
{p_end}
{synopt:{opt d0(varname)}}
estimated E[D|X,Z=0] (model {opt interactiveiv})
{p_end}
{synopt:{opt d1(varname)}}
estimated E[D|X,Z=1] (model {opt interactiveiv})
{p_end}


{smcl}
INCLUDE help ddml_install_ref_auth
