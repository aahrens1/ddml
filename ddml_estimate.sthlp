{smcl}
{* *! version 11july2023}{...}
{hline}
{cmd:help ddml crossfit, ddml estimate}{right: v1.3}
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
Before cross-fitting, the model must be defined using {helpb ddml init} and the learners specified using {helpb ddml eq};
see the help files for details.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml crossfit} [ , {opt mname(name)} {opt shortstack} {opt poolstack}
{opt ssfinalest(name)} {opt psfinalest(name)} {opt finalest(name)}{bind: ]} 

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)} {opt atet} {opt ateu} {opt trim(real)}{bind: ]} 

{p 8 14}Replay options:

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {opt spec(integer or string)} {opt rep(integer or string)} {opt allcombos} {opt not:able} {opt replay} {bind: ]} 


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
available when {help pystacked} is used to specify the base learners
{p_end}
{synopt:{opt ssfinalest(name)}} specifies the estimator used
to obtain the short-stacking weights; the default is constrained non-negative least squares.
For the list of available final estimators, see the help for {help pystacked}.
{p_end}
{synopt:{opt poolstack}} is available as an alternative to standard stacking
when {help pystacked} is used to specify the base learners.
Pooled-stacking adds additional regularization
by obtaining a single set of stacking weights
from the full set of out-of-sample base learner predicted values
(in contrast to {help pystacked}, which stacks each cross-fit fold separately).
{p_end}
{synopt:{opt psfinalest(name)}} specifies the estimator used
to obtain the pooled-stacking weights; the default is constrained non-negative least squares.
For the list of available final estimators, see the help for {help pystacked}.
{p_end}
{synopt:{opt finalest(name)}} sets the final estimator for both
short-stacking and pooled-stacking.
NB: to set the final estimator for standard stacking using {help pystacked},
use the {help pystacked} {opt finalest} option when specifying the base learners using {help ddml eq}.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

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
select variance-covariance estimator; see {helpb regress##vcetype:here}.
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
{synoptline}
{p2colreset}{...}
{pstd}

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


{title:Examples}

{pstd}
For more examples of usage see the links via the main {help ddml##examples:ddml help file}.
See {help ddml init:help ddml init} for details of model initialization and learner specification options.

{pstd}Note: the additional support provided by {opt ddml} for {helpb pystacked} (see {help ddml stacking:help ddml stacking})
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


{marker installation}{title:Installation}

{pstd}
To get the latest stable version of {cmd:ddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.


{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}


{pstd}
To verify that {cmd:ddml} is correctly installed, 
click on or type {stata "whichpkg ddml"} 
(which requires {helpb whichpkg} 
to be installed; {stata "ssc install whichpkg"}).


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
Help: {helpb pystacked}, {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}.{p_end}
