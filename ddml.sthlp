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
(Economics Letters, 2018). Five different models are supported; allowing for 
binary or continous treatment variables and endogeneity. 
{opt ddml} supports a variety of different ML programs, including
but not limited to {helpb lassopack} and {helpb pylearn}. 

{marker syntax}{...}
{title:Syntax}

{pstd}
{ul:Step 1:} Initialise {cmd:ddml} and select model:

        {cmd:ddml init} {it:model}

{pstd}
where {it:model} is either {it:partial}, 
{it:iv}, {it:interactive}, {it:optimaliv}, {it:late};
see model descriptions below.

{pstd}
{ul:Step 2:} Add supervisd ML programs for estimating conditional expectations:

        {cmd:ddml} {it:eq} {it:newvarname} [, {it:eqopt}]: {it:command} {it:depvar} {it:vars} [, {it:cmdopt}]

{pstd}
where {it:eq} is either {it:yeq}, {it:deq} or {it:zeq}. {it:command} is a
ML program that support the standard {it:reg y x}-type syntax. 
{it:cmdopt} are specific to that program.
See compatibility below.

{pstd}
{ul:Step 3:} Cross-fitting

        {cmd:ddml crossfit} [, {it:crossfitopt}] 

{pstd}
{ul:Step 4:} Estimate causal effects

        {cmd:ddml estimate} [, {it:estopt}] 

{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:eqopt}
{synoptline}
{synopt:{cmdab:postl:ogit}}
???
{p_end}
{synopt:{cmdab:nocon:stant}}
???
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:crossfitopt}
{synoptline}
{synopt:{opt k(integer)}}
number of folds for cross-fitting / sample splitting. 
The default is 2. The theory of DDML does not
depend on the number of folds; yet, we recommend 
to consider test higher number of folds (e.g., 5, 10)
to check robustness of your results.
{p_end}
{synopt:{opt absorb(varlist)}}
partial out fixed effects
{p_end}
{synopt:{opt tabf:old}}
show number of observations per fold
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:estopt}
{synoptline}
{synopt:{opt robust}}
report SEs that are robust to the
presence of arbitrary heteroskedasticity
{p_end}
{synopt:{opt show(string)}}
all if all combinations should be estimated. default is opt which 
only shows the optimal combatination. 
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{marker models}{...}
{title:Models}

{pstd}
The following models are implemented: 

{pstd}
{ul:{it:Partial linear model:}}

	Y = {it:a}.D + g(X) + U
        D = m(X) + V

{pstd}
where the aim is to estimate {it:a} while controlling for X.

{pstd}
{ul:{it:Interactive model:}}

	Y = g(X,D) + U
        D = m(X) + V

{pstd}
where we are, as in the Partial Linear Model interested in the ATE, but do not 
assume that X and D are separable.

{pstd}
{ul:{it:Partial linear IV model:}}

	Y = {it:a}.D + g(X) + U
        Z = m(X) + V

{pstd}
where the aim is to estimate the average treatment effect or local average treatment effect.

{pstd}
{ul:{it:LATE model:}}

	Y = g(Z,X) + U
        D = m(Z,X) + V
        Z = m(X) + V

{pstd}
where the aim is to estimate the average treatment effect or local average treatment effect.

{pstd}
{ul:{it:Optimal IV model:}}

	Y = {it:a}.D + U
        D = m(Z) + V 

{pstd}
where the aim is to estimate the average treatment effect or local average treatment effect.


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
- {helpb rforest} by Zou & Schonlau.

{pstd}
Beyond these, it is compatible with almost any Stata program that uses the standard {it:reg y x}-type Syntax
and comes with predict post-estimation programs.

{pstd}
If you are aware of a program that is not compatible with {opt ddml}, but think it should be, please
do not hesitate to contact us.

{marker examples}{...}
{title:Examples}

{pstd}
This help file is under construction.

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

{title:Also see}

{pstd}
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso}, {helpb pdslasso}, {helpb pylearn}.{p_end}