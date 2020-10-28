{smcl}
{* *! version 18sep2020}{...}
{hline}
{cmd:help ddml}{right: v0.1.2}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{hi: ddml} {hline 2}}Stata package for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters.}
 (Economics Letters, 2018).

{marker syntax}{...}
{title:Syntax}

{marker models}{...}
{title:Models}

The following models are implemented: 

{pstd}
{it:Partial linear model:}

	Y = {it:a}.D + g(X) + U
        D = m(X) + V

{pstd}
where the aim is to estimate {it:a} while controlling for X.

{pstd}
{it:Interactive model:}

	Y = g(X,D) + U
        D = m(X) + V

{pstd}
where we are, as in the Partial Linear Model interested in the ATE, but do not 
assume that X and D are separable.

{pstd}
{it:IV model:}

	Y = g(Z,X) + U
        D = m(Z,X) + V
        Z = m(X) + V

{pstd}
where the aim is to estimate the average treatment effect or local average treatment effect.

{marker compatibility}{...}
{title:Compatible programs}

{opt ddml} is compatible with a large set of user-written Stata commands. 
It has been tested with 

{pstd} 
{opt lassopack} for regularized regression (see {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}).

{pstd}  
the {helpb pylearn} package by Michael Droste (see {helpb pytree}, {helpb pyforest}, {helpb pymlp}, {helpb pyadaboost}, {helpb pygradboost}). 
Note that {helpb pylearn} requires Stata 16.

{pstd} 
{helpb rforest} by Zou & Schonlau.

Beyond these, it is compatible with almost any Stata program that uses the standard {it:reg y x}-type Syntax
and comes with predict post-estimation programs.

If you are aware of a program that is not compatible with {opt ddml}, but think it should be, please
do not hesitate to contact us.

{marker examples}{...}
{title:Examples}

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

{p 7 14 2}
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso}, {helpb pdslasso} (if installed).{p_end}
