{smcl}
{* *! version 15jun2023}{...}
{hline}
{cmd:help ddml crossfit, ddml estimate}{right: v1.2}
{hline}

{title:ddml crossfit and estimate commands for Double Debiased Machine Learning}

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

{p 8 14}{cmd:ddml crossfit} [ , {opt mname(name)} {opt shortstack}{bind: ]} 

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
Short-stacking runs constrained non-negative least squares on the
cross-fitted predicted values to obtain a weighted average
of several base learners.
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
estimates all possible specifications. By default, only the min-MSE (or short-stacking)
specification is estimated and displayed.
{p_end}
{synopt:{opt replay}}
used in combination with {opt spec()} and {opt rep()} to display and return estimation results.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{title:Examples}

For examples of usage see {help ddml##examples:help ddml}.


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

{marker installation}{title:Installation}

{pstd}
To get the latest stable version of {cmd:ddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.

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
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}, {helpb pystacked}.{p_end}
