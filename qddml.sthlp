{smcl}
{* *! version 18sep2020}{...}
{hline}
{cmd:help ddml}{right: v0.1.2}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{hi: qddml} {hline 2}}Stata program for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
{opt qddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported; allowing for 
binary or continous treatment variables and endogeneity. 
{opt qddml} uses stacking regression ({helpb pystacked}) as the default machine learning algorithm. 
{opt qddml} is a wrapper for {opt ddml}, which offers more flexibility.

{p 8 14 2}
{cmd:qddml}
{it:depvar} {it:regressors} [{cmd:(}{it:hd_controls}{cmd:)}]
{cmd:(}{it:endog}{cmd:=}{it:instruments}{cmd:)}
[{cmd:if} {it:exp}] [{cmd:in} {it:range}]
{opt model(name)}
{bind:[ {cmd:,}}
{opt cmd(varlist)}
{opt cmdopt(varlist)}
{opt mname(string)}
{bind:{cmdab:noc:onstant} ]}

{pstd}
Since {opt qddml} uses {helpb pystacked} per default, 
it requires Stata 16 and Python in the default setting. See HERE for how to set up
Python your system.
If you don't have Stata 16,
you need to change the ML program used for estimating conditional expectations, e.g., 
using the option {opt cmd(rlasso)}.

{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:Option}
{synoptline}
{synopt:{opt model(name)}}
the model to be estimated; allows for {it:partial}, {it:interactive},
{it:iv}, {it:optimaliv}, {it:late}. See {helpb ddml##models:here} for a description.
{p_end}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:ddml1}.
{p_end}

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
 {helpb pdslasso}, {helpb pylearn}.{p_end}
