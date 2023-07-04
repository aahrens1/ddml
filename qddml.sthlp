{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help qddml}{right: v1.2}
{hline}

{title:ddml, qddml - Stata package for Double Debiased Machine Learning}

{pstd}
{help ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 
{help ddml} supports a variety of different ML programs, including
but not limited to {help pystacked} and {help lassopack}. 

{pstd}
{opt qddml} is a wrapper program for {help ddml}. It provides a convenient 
one-line syntax with almost the full flexibility of {help ddml}.
The main restriction of {cmd:qddml} is that it only allows to be used 
with one machine learning program at the time, while {help ddml} 
allow for multiple learners per reduced form equation.

{pstd}
{opt qddml} uses stacking regression implemented via {help pystacked}
as the default machine learning program.
{opt qddml} has various options specific to {help pystacked}; see below.

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

{pstd}
Please check the {helpb qddml##examples:examples} provided at the end of this help file.


{marker syntax}{...}
{title:Syntax}

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
{bind:... ]}


{title:Options}

{synoptset 20}{...}
{synopthdr:General}
{synoptline}
{synopt:{opt model(name)}}
the model to be estimated; allows for {it:partial}, {it:interactive},
{it:iv}, {it:fiv}, {it:late}. See {helpb ddml##models:here} for an overview.
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
{synopthdr:Stacking}
{synoptline}
{synopt:{opt shortstack}} asks for short-stacking to be used.
Short-stacking combines the cross-fitted predicted values
to obtain a weighted average of several base learners.
The default behavior of {opt qddml} is to report short-stacked results;
the {opt shortstack} option is used when the user wants to report
more than one type of stacking results.
{p_end}
{synopt:{opt ssfinalest(name)}}
specified the final (ensemble) estimator used to combine
the base learner short-stacked cross-fitted predicted values.
The default is constrained non-negative least squares;
for alternatives, see the choices available for standard stacking in {help pystacked}.
{p_end}
{synopt:{opt stdstack}} requests results based on standard stacking via {help pystacked}.
Stacking is done k times when cross-fitting.
Note that the final estimator for combining the base learner predicted values
is set via the options provided directly to {help pystacked};
the {help pystacked} default is constrained non-negative least squares.
Standard stacking is available only in combination with {help pystacked}.
{p_end}
{synopt:{opt poolstack}} is an alternative to standard stacking via {help pystacked}.
Pooled-stacking combines the out-of-sample predicted values of the {help pystacked}
using a single final estimation and set of weights.
{p_end}
{synopt:{opt psfinalest(name)}}
specified the final (ensemble) estimator used to combine
the base learner cross-fitted predicted values when pooled-stacking.
The default is constrained non-negative least squares;
for alternatives, see the choices available for standard stacking in {help pystacked}.
Pooled-stacking is available only in combination with {help pystacked}.
{p_end}

{pstd}
For more on stacking, see {help ddml stacking:help ddml stacking}.
{p_end}

{synoptset 20}{...}
{synopthdr:Learners - pystacked}
{synoptline}
{synopt:{opt pystacked(string)}}
{help pystacked} options;
applies to all estimated conditional expectations. 
If no options are specified, {help pystacked}'s defaults are used.
{p_end}
{synopt:{opt pystacked_y(string)}}
{help pystacked} options specific to the estimation of the conditional expectations of the outcome {it:Y}. 
{p_end}
{synopt:{opt pystacked_d(string)}}
{help pystacked} options specific to the estimation of the conditional expectations of the treatment variable(s) {it:D}.
{p_end}
{synopt:{opt pystacked_z(string)}}
{help pystacked} options specific to the estimation of the conditional expectations of instrumental variable(s) {it:Z}. 
{p_end}

{synoptset 20}{...}
{synopthdr:Learners - other}
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
{synopt:{opt *cmdopt(string)}}
options that are passed on to ML program. 
The asterisk {cmd:*} can be replaced with either nothing 
(setting the default for all reduced form equations), 
{cmd:y} (setting the default for the conditional expectation of {it:Y}), 
{cmd:d} (setting the default for {it:D})
or {cmd:z} (setting the default for {it:Z}).
{p_end}
{synopt:{opt *vtype(string)}}
variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank 
(this is required when using {help ddml} with {helpb rforest}.)
The asterisk {cmd:*} can be replaced with either nothing 
(setting the default for all reduced form equations), 
{cmd:y} (setting the default for the conditional expectation of {it:Y}), 
{cmd:d} (setting the default for {it:D})
or {cmd:z} (setting the default for {it:Z}).
{p_end}
{synopt:{opt *predopt(string)}}
{cmd:predict} option to be used to get predicted values. 
Typical values could be {opt xb} or {opt pr}. Default is 
blank. The asterisk {cmd:*} can be replaced with either nothing 
(setting the default for all reduced form equations), 
{cmd:y} (setting the default for the conditional expectation of {it:Y}), 
{cmd:d} (setting the default for {it:D})
or {cmd:z} (setting the default for {it:Z}).
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
See {helpb ddml##help:here}.


{marker compatibility}{...}
{title:Compatible programs}

{pstd} 
See {helpb ddml##compatibility:here}.


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:qddml} for each of the 5 models supported. 
Note that estimation models are chosen for demonstration purposes only and 
kept simple to allow you to run the code quickly.
Please also see the examples via the links in the {helpb ddml##examples:ddml help file}

{pstd}{ul:Partially linear model.} 

{pstd}Preparations: we load the data and define global macros.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}

{pstd}The options {cmd:model(partial)} selects the partially linear model
and {cmd:kfolds(2)} selects two cross-fitting folds.
We use the options {cmd:cmd()} and {cmd:cmdopt()} to select
random forest and cross-validated lasso for estimating the conditional expectations.{p_end}

{pstd}Note that we set the number of random folds to 2, so that 
the model runs quickly. The default is {opt kfolds(5)}. We recommend 
to consider at least 5-10 folds and even more if your sample size is small.{p_end}

{pstd}Note also that we recommend to re-run the model multiple time on 
different random folds; see options {opt reps(integer)}.{p_end}

{phang2}. {stata "qddml $Y $D ($X), kfolds(2) model(partial) cmd(pystacked) cmdopt(type(reg) method(rf lassocv))"}{p_end}

{pstd}Postestimation options for {helpb ddml} also work after {opt qddml}.
Here we request display of the {helpb pystacked} stacking weights and MSEs.{p_end}

{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}

{pstd}{ul:Equivalence of ddml and qddml.} 

{pstd}Here we illustrate how to replicate the output of qddml using separate ddml command.
Note that to guarantee exact replication, we need to set Stata's random-number seed
and also specify two {help pystacked}-specific options.
The {opt noshuffle} option is relevant to the cross-validation estimators supported by {help pystacked};
it stops the data being shuffled before splitting into cross-validation folds.
The {opt randomstate{int)} option is specific to the {help pystacked} random forest estimator;
it controls the randomness of the bootstrap samples and sampling of features.
We also need to ensure that we use the same cross-fold split in both estimations.
{p_end}

{pstd}We specify 5 base learners: OLS, cross-validated lasso and ridge, and two random forests.
We will use the same specification for both conditional expectations E[Y|X} and E[D|X].{p_end}
{phang2}. {stata "global rflow max_features(5) min_samples_leaf(1) max_samples(.7) random_state(0)"}{p_end}
{phang2}. {stata "global rfhigh max_features(5) min_samples_leaf(10) max_samples(.7) random_state(0)"}{p_end}
{phang2}. {stata "global psoptions method(ols lassocv ridgecv rf rf) cmdopt4($rflow) cmdopt5($rfhigh) type(reg) noshuffle"}{p_end}
pstd}Estimation using {opt ddml}. Note that initialization also generates the (random) cross-fit split and identifier.{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, $psoptions"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X, $psoptions"}{p_end}
{pstd}Reset the seed prior to cross-fitting and estimation:{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{pstd}Cross-fit and estimate:{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}
{pstd}Estimation using {opt qddml}. Note that we re-use the same cross-fit split via the {opt foldvar(varname)} option.{p_end}
{phang2}. {stata "qddml $Y $D ($X), mname(m1) model(partial) foldvar(m0_fid_1) pystacked($psoptions)"}{p_end}

{pstd}{ul:Partially linear IV model.} 

{pstd}Preparations: we load the data and define global macros.{p_end}
{phang2}. {stata "use https://statalasso.github.io/dta/AJR.dta, clear"}{p_end}
{phang2}. {stata "global Y logpgp95"}{p_end}
{phang2}. {stata "global D avexpr"}{p_end}
{phang2}. {stata "global Z logem4"}{p_end}
{phang2}. {stata "global X lat_abst edes1975 avelf temp* humid* steplow-oilres"}{p_end}

{pstd}Since the data set is very small, we consider 30 cross-fitting folds.{p_end} 
{pstd}We need to add the option {opt vtype(none)} for {helpb rforest} to 
work with {help ddml} since {helpb rforests}'s {cmd:predict} command doesn't
support variable types.{p_end}

{phang2}. {stata "qddml $Y ($X) ($D=$Z), kfolds(30) model(iv) cmd(rforest) cmdopt(type(reg)) vtype(none) robust"}{p_end}

{pstd}{ul:Interactive model--ATE and ATET estimation.} 

{pstd}Preparations: we load the data and define global macros.{p_end}
{phang2}. {stata "webuse cattaneo2, clear"}{p_end}
{phang2}. {stata "global Y bweight"}{p_end}
{phang2}. {stata "global D mbsmoke"}{p_end}
{phang2}. {stata "global X mage prenatal1 mmarried fbaby mage medu"}{p_end}

{pstd}
Note that we use gradient boosted regression trees for E[Y|X,D] (see {opt ycmdopt()}),
but gradient boosted classification trees for E[D|X] (see {opt dcmdopt()}).
{p_end} 
{phang2}. {stata "qddml $Y $D ($X), kfolds(5) reps(5) model(interactive) cmd(pystacked) ycmdopt(type(reg) method(gradboost)) dcmdopt(type(class) method(gradboost))"}{p_end}

{pstd}{cmd:qddml} reports the ATE effect by default. The option {cmd:atet}
returns the ATET estimate.{p_end}

{pstd}If we want retrieve the ATET estimate after estimation, 
we can simply use {ddml estimate}.{p_end}
{phang2}. {stata "ddml estimate, atet"}{p_end}

{pstd}{ul:Interactive IV model--LATE estimation.} 

{pstd}Preparations: we load the data and define global macros.{p_end}
{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta,clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}

{phang2}. {stata "qddml $Y (c.($X)# #c($X)) ($D=$Z), kfolds(5) model(interactiveiv) cmd(pystacked) ycmdopt(type(reg) m(lassocv)) dcmdopt(type(class) m(lassocv)) zcmdopt(type(class) m(lassocv))"}{p_end}

{pstd}{ul:Flexible Partially Linear IV model.} 

{pstd}Preparations: we load the data and define global macros.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear"}{p_end}
{phang2}. {stata "global Y share"}{p_end}
{phang2}. {stata "global D price"}{p_end}
{phang2}. {stata "global X hpwt air mpd space"}{p_end}
{phang2}. {stata "global Z sum*"}{p_end}

{pstd}The syntax is the same as in the Partially Linear IV model, 
but we now estimate the optimal instrument flexibly.{p_end}
{phang2}. {stata "qddml $Y ($X) ($D=$Z), model(fiv)"}{p_end}


{marker installation}{title:Installation}

{pstd}
To get the latest stable version of {help ddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.

{pstd}
To verify that {help ddml} is correctly installed, 
click on or type {stata "whichpkg ddml"} 
(which requires {helpb whichpkg} 
to be installed; {stata "ssc install whichpkg"}).


{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}


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
