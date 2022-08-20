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
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 
{opt ddml} supports a variety of different ML programs, including
but not limited to {helpb lassopack} and {helpb pystacked}. 

{pstd}
{opt qddml} is a wrapper program of {cmd:ddml}. It provides a convenient 
one-line syntax with almost the full flexibility of {cmd:ddml}.

{pstd}
{opt qddml} uses stacking regression ({helpb pystacked}) as the default machine learning algorithm. 

{p 8 14 2}
{cmd:qddml}
{it:depvar} {it:regressors} [{cmd:(}{it:hd_controls}{cmd:)}]
{cmd:(}{it:endog}{cmd:=}{it:instruments}{cmd:)}
[{cmd:if} {it:exp}] [{cmd:in} {it:range}]
{opt model(name)}
{bind:[ {cmd:,}}
{opt cmd(string)}
{opt cmdopt(string)}
{opt mname(string)}
{opt ...} ]}

{pstd}
Since {opt qddml} uses {helpb pystacked} per default, 
it requires Stata 16 or higher, Python 3.x and at least scikit-learn 0.24. See 
{helpb python:this help file}, {browse "https://blog.stata.com/2020/08/18/stata-python-integration-part-1-setting-up-stata-to-use-python/":this Stata blog entry}
and 
{browse "https://www.youtube.com/watch?v=4WxMAGNhcuE":this Youtube video}
for how to set up
Python on your system.
In short, install Python 3.x (we recommend Anaconda) 
and set the appropriate Python path using {cmd:python set exec}.
If you don't have Stata 16+,
you can still use {cmd:pystacked} with programs that don't rely on Python, 
e.g., using the option {opt cmd(rlasso)}.

{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:General}
{synoptline}
{synopt:{opt model(name)}}
the model to be estimated; allows for {it:partial}, {it:interactive},
{it:iv}, {it:ivhd}, {it:late}. See {helpb ddml##models:here} for an overview.
{p_end}
{synopt:{opt mname(string)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt kfolds(integer)}}
number of cross-fitting folds. The default is 5.
{p_end}
{synopt:{opt fcluster(varname)}}
cluster identifiers for cluster randomization of random folds.
{p_end}
{synopt:{opt foldvar(varname)}}
integer variable with user-specified cross-fitting folds.
{p_end}
{synopt:{opt reps(integer)}}
number of re-sampling iterations, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds. 
{p_end}
{synopt:{opt shortstack}} asks for short-stacking to be used.
Short-stacking runs contrained non-negative least squares on the
cross-fitted predicted values to obtain a weighted average
of several base learners.
{p_end}
{synopt:{cmdab:r:obust}}
report SEs that are robust to the
presence of arbitrary heteroskedasticity.
{p_end}
{synopt:{opt vce(type)}}
select variance-covariance estimator, see {helpb regress##vcetype:here}
{p_end}
{synopt:{opt cluster(varname)}}
select cluster-robust variance-covariance estimator.
{p_end}

{synoptset 20}{...}
{synopthdr:Learners}
{synoptline}
{synopt:{opt cmd(string)}}
ML program used for estimating conditional expectations. 
Defaults to {helpb pystacked}. 
See {helpb ddml##compatibility:here} for 
other supported programs.
{p_end}
{synopt:{opt ycmd(string)}}
ML program used for estimating the conditional expectations of the outcome {it:Y}. 
Defaults to {opt cmd(string)}. 
{p_end}
{synopt:{opt dcmd(string)}}
ML program used for estimating the conditional expectations of the treatment variable(s) {it:D}. 
Defaults to {opt cmd(string)}. 
{p_end}
{synopt:{opt zcmd(string)}}
ML program used for estimating conditional expectations of instrumental variable(s) {it:Z}. 
Defaults to {opt cmd(string)}. 
{p_end}
{synopt:{opt cmdopt(string)}}
options that are passed on to ML program
{p_end}
{synopt:{opt ycmdopt(string)}}
options that are passed on to ML program used for 
conditional expectations of the outcome {it:Y}. 
{p_end}
{synopt:{opt dcmdopt(string)}}
options that are passed on to ML program used for 
conditional expectations of the treatment variable(s) {it:D}. 
{p_end}
{synopt:{opt zcmdopt(string)}}
options that are passed on to ML program used for 
conditional expectations of instrumental variable(s) {it:Z}. 
{p_end}
{synopt:{opt *vtype(string)}}
variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank 
(this is required when using {cmd:ddml} with {helpb rforest}.)
Replace {cmd:*} with either {cmd:y}, {cmd:d} or {cmd:z} to pass
option to learner of conditional expectations of {it:Y}, 
{it:D} or {it:Z}.
{p_end}
{synopt:{opt *predopt(string)}}
{cmd:predict} option to be used to get predicted values. 
Typical values could be {opt xb} or {opt pr}. Default is 
blank. 
Replace {cmd:*} with either {cmd:y}, {cmd:d} or {cmd:z} to pass
option to learner of conditional expectations of {it:Y}, 
{it:D} or {it:Z}.
{p_end}

{synoptset 20}{...}
{synopthdr:Output}
{synoptline}
{synopt:{opt verb:ose}}
show detailed output
{p_end}
{synopt:{opt vverb:ose}}
show even more output
{p_end}

{marker models}{...}
{title:Models}

{pstd} 
See {helpb ddml##models:here}.

{marker compatibility}{...}
{title:Compatible programs}

{pstd} 
See {helpb ddml##compatibility:here}.

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
