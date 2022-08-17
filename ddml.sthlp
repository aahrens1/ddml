{smcl}
{* *! version 18sep2020}{...}
{hline}
{cmd:help ddml}{right: v0.1.2}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{hi: ddml} {hline 2}}Stata package for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continous treatment variables and endogeneity. 
{opt ddml} supports a variety of different ML programs, including
but not limited to {helpb lassopack} and {helpb pystacked}. 

{marker syntax}{...}
{title:Syntax}

The estimation proceeds in fours steps.

{pstd}
{ul:Step 1:} Initialise {cmd:ddml} and select model:

{p 8 14}{cmd:ddml init}
{it:model} 
[, {opt mname(string)} {opt kfolds(integer)}
{opt foldvar(varlist)} {opt reps(integer)} 
{opt tabfold} {opt vars(varlist)}]

{pstd}
where {it:model} is either {it:partial}, 
{it:iv}, {it:interactive}, {it:optimaliv}, {it:late};
see model descriptions below.

{pstd}
{ul:Step 2:} Add supervised ML programs for estimating conditional expectations:

{p 8 14}{cmd:ddml} {it:eq} 
[, {opt mname(string)} {opt vname(string)} {opt vtilde(string)}
{opt vtype(string)}
{opt predopt(string)}]:
{it:command} {it:depvar} {it:vars} [, {it:cmdopt}]

{pstd}
where, depending on the chosen model,
{it:eq} is either 
{it:E[Y|X]} {it:E[Y|D,X]} {it:E[Y|X,Z]} {it:E[D|X]} {it:E[D|X,Z]} {it:E[Z|X]}.
{it:command} is a ML program that supports the standard {it:reg y x}-type syntax,
e.g. {helpb pystacked}. 
{it:cmdopt} are specific to that program.
See compatibility below.

{pstd}
{ul:Step 3:} Cross-fitting

{p 8 14}{cmd:ddml crossfit} [, {opt mname(string)} {opt shortstack}] 

{pstd}
{ul:Step 4:} Estimate causal effects

{p 8 14}{cmd:ddml estimate} [, {opt mname(string)} {opt spec(integer)} {opt rep(integer)} {cmdab:r:obust} {opt clustervar(varname)} {opt vce(type)}] 

{pstd}
{ul:Auxiliary sub-programs:}
     
{pstd} 
Update {cmd:ddml}:

{p 8 14}{cmd:ddml update} 

{pstd}
Describe model:

{p 8 14}{cmd:ddml desc} [, {opt mname()}]

{pstd}
Set-up Mata library:

{p 8 14}{cmd:ddml setup} 

{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:init options}
{synoptline}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt kfolds(integer)}}
number of cross-fitting folds. The default is 5.
{p_end}
{synopt:{opt foldvar(varname)}}
integer variable with user-specified cross-fitting folds.
{p_end}
{synopt:{opt reps(integer)}}
number of re-sampling iterations, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds. 
{p_end}
{synopt:{opt tabfold}}
prints a table with frequency of observations by default.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:Equation options}
{synoptline}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt vname(string)}}
name of the dependent variable in the reduced form estimation. 
This is usually inferred from the command line but is mandatory
for the {it:ivhd} model.
{p_end}
{synopt:{opt vtilde(string)}}
name of the variable to be created. 
{p_end}
{synopt:{opt vtype(string)}}
variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank.
{p_end}
{synopt:{opt predopt(string)}}
{cmd:predict} option to be used to get predicted values. 
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:Cross-fitting}
{synoptline}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt shortstack}} asks for short-stacking to be used.
Short-stacking runs contrained non-negative least squares on the
cross-fitted predicted values to obtain a weighted average
of several base learners.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:Estimation}
{synoptline}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt spec(integer)}}
select specification
{p_end}
{synopt:{opt spec(integer)}}
select resampling iteration
{p_end}
{synopt:{opt vce(type)}}
select variance-covariance estimator, see {helpb regress##vcetype:here}
{p_end}
{synopt:{cmdab:r:obust}}
report SEs that are robust to the
presence of arbitrary heteroskedasticity.
{p_end}
{synopt:{opt clustervar(varname)}}
select cluster-robust variance-covariance estimator.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{marker models}{...}
{title:Models}

{pstd}
The following models are implemented: 

{pstd}
{ul:Partial linear model} [{it:partial}]

	Y = {it:a}.D + g(X) + U
        D = m(X) + V

{pstd}
where the aim is to estimate {it:a} while controlling for X. To this end, 
we estimate the conditional expectations
E[Y|X] and E[D|X] using a supervised machine learner.

{pstd}
{ul:Interactive model} [{it:interactive}]

	Y = g(X,D) + U
        D = m(X) + V

{pstd}
which relaxes the assumption that X and D are separable. 
D is a binary treatment variable. 
We estimate the conditional expectations 
E[Y|X,D=0] and E[Y|X,D=1] as well as 
E[D|X] using a supervised machine learner.

{pstd}
{ul:Partial linear IV model} [{it:iv}]

	Y = {it:a}.D + g(X) + U
        Z = m(X) + V

{pstd}
where the aim is to estimate {it:a}. 
We estimate the conditional expectations E[Y|X], 
E[D|X] and E[Z|X] using a supervised machine
learner.

{pstd}
{ul:Interactive IV model}  [{it:late}]

	Y = g(Z,X) + U
        D = h(Z,X) + V
        Z = m(X) + E

{pstd}
where the aim is to estimate the local average treatment effect.
We estimate, using a supervised machine
learner, the following conditional expectations:
E[Y|X,Z=0] and E[Y|X,Z=1] (jointly referred to as "{opt yeq}");
E[D|X,Z=0] and E[D|X,Z=1] (jointly referred to as "{opt deq}");
E[Z|X] ("zeq").

{pstd}
{ul:High-dimensional IV model} [{it:ivhd}]

	Y = {it:a}.D + g(X) + U
        D = m(Z) + g(X) + V 

{pstd}
where the estimand of interest is {it:a}. 
We estimate the conditional expectations
E[Y|X], 
E[D^|X] and D^:=E[D|Z,X] using a supervised machine
learner. The instrument is then formed as D^-E^[D^|X] where E^[D^|X] denotes
the estimate of E[D^|X]. 

{marker compatibility}{...}
{title:Compatible programs}

{pstd}
{opt ddml} is compatible with a large set of user-written Stata commands. 
It has been tested with 

{p 7 9 0} 
- {opt lassopack} for regularized regression (see {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}).

{p 7 9 0} 
- the {helpb pylearn} package by Michael Droste (see {helpb pytree}, {helpb pyforest}, {helpb pymlp}, {helpb pyadaboost}, {helpb pygradboost}). 
Note that {helpb pylearn} requires Stata 16.

{p 7 9 0} 
- the {helpb pystacked} package (see {helpb pystacked}. 
Note that {helpb pystacked} requires Stata 16.

{p 7 9 0} 
- {helpb rforest} by Zou & Schonlau.

{p 7 9 0} 
- {helpb svmachines} by Guenther & Schonlau.

{pstd}
Beyond these, it is compatible with almost any Stata program that uses the standard {it:reg y x}-type Syntax,
supports {it:if}-conditions and comes with predict post-estimation programs.

{pstd}
If you are aware of a program that is not compatible with {opt ddml}, but think it should be, please
do not hesitate to contact us.

{marker examples}{...}
{title:Examples}

{pstd}
To be added.

{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

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
Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland
achim.ahrens@gess.ethz.ch

{pstd}
Christian B. Hansen, University of Chicago, USA
Christian.Hansen@chicagobooth.edu

{pstd}
Mark E Schaffer, Heriot-Watt University, UK
m.e.schaffer@hw.ac.uk	

{title:Also see (if installed)}

{pstd}
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}, {helpb pystacked}.{p_end}
