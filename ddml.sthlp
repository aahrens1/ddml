{smcl}
{* *! version 11july2023}{...}
{hline}
{cmd:help ddml}{right: v1.3}
{hline}

{title:ddml - Stata package for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml} supports a variety of different ML programs, including
but not limited to {helpb pystacked} and {helpb lassopack}. 
{helpb pystacked} is the recommended way to specify multiple learners in {opt ddml},
and {opt ddml} has integrated support for various features provided by {helpb pystacked}.

{pstd}
The {opt ddml} package includes the wrapper program {helpb qddml},
which uses a simplified one-line syntax, 
but offers less flexibility.

{pstd}
{opt ddml} and {opt qddml} rely on {helpb crossfit}, which can be used as a standalone program.


{title:Contents}

{p 2}{help ddml##help:Links to detailed help files}{p_end}
{p 2}{help ddml##overview:Overview of main steps to estimation using ddml}{p_end}
{p 2}{help ddml##compatibility:Compatible programs and pystacked integration}{p_end}
{p 2}{help ddml##installation:Installation}{p_end}
{p 2}{help ddml##references:References}{p_end}
{p 2}{help ddml##authors:Authors}{p_end}


{marker help}{...}
{title:Follow links below to detailed help}

{p2colset 5 25 25 0}{...}
{p 2}Help files: qddml{p_end}
{p2col:{help qddml}} One-step DDML estimation using {help qddml}.{p_end}

{p 2}Help files: main steps in using {opt ddml}{p_end}
{p2col:{help ddml overview}} Overview of main steps to estimation using ddml{p_end}
{p2col:{help ddml stacking}} Overview of stacking in {opt ddml} and {help pystacked}{p_end}

{p 2}Help files: detailed help files for main steps in using {opt ddml}{p_end}
{p2col:{help ddml init}} 1. Initialize {opt ddml} and select model.{p_end}
{p2col:{help ddml eq}} 2. Add supervised ML programs for estimating conditional expectations.{p_end}
{p2col:{help ddml crossfit}} 3. Cross-fitting to estimate conditional expectations.{p_end}
{p2col:{help ddml estimate}} 4. Estimate causal model and report/post results.{p_end}

{marker examples}{...}
{p 2}Help files: detailed help files and examples for models supported by {opt ddml}{p_end}
{p2col:{help ddml partial}} Partially-linear model{p_end}
{p2col:{help ddml iv}} Partially-linear IV model{p_end}
{p2col:{help ddml fiv}} Flexible partially-linear IV model{p_end}
{p2col:{help ddml interactive}} Interactive model - ATE and ATET estimation{p_end}
{p2col:{help ddml interactiveiv}} Interactive IV model - LATE estimation{p_end}

{p 2}Help files: auxiliary programs{p_end}
{p2col:{help ddml describe}} Report information about the model setup and/or results.{p_end}
{p2col:{help ddml extract}} Report information about saved results e.g. stacking weigths.{p_end}
{p2col:{help ddml sample}} Report information about the estimation sample, folds, etc.{p_end}
{p2col:{help ddml export}} Save the {opt ddml} estimated conditional expectations to a csv file.{p_end}
{p2col:{help ddml overlap}} (interactive models only) Generate overlap plots for propensity-score-based models{p_end}
{p2col:{help crossfit}} Use {opt crossfit} as a standalone program for cross-fitting and cross-validation.{p_end}


{marker overview}{...}
{title:Overview of main steps when estimating with ddml}

{pstd}Estimation with {cmd:ddml} proceeds in four steps. 

{pstd}
{ul:Step 1.} Initialize {cmd:ddml} and select model:

{p 8 14}{cmd:ddml init}
{it:model} [if] [in]
[ , {opt mname(name)} {opt kfolds(integer)}
{opt fcluster(varname)}
{opt foldvar(varlist)} {opt reps(integer)} 
{opt norandom} {opt tabfold} {opt vars(varlist)}{bind: ]}

{pstd}
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, {it:interactiveiv};
see {helpb ddml##models:model descriptions}.

{pstd}
{ul:Step 2.} Add supervised ML programs for estimating conditional expectations:

{p 8 14}{cmd:ddml} {it:eq} 
[ , {opt mname(name)} {opt vname(varname)} {opt l:earner(varname)}
{opt vtype(string)}
{opt predopt(string)}{bind: ] :}
{it:command} {it:depvar} {it:vars} [ , {it:cmdopt}{bind: ]}

{pstd}
where, depending on model chosen in Step 1,
{it:eq} is either 
{it:E[Y|X]} {it:E[Y|D,X]} {it:E[Y|X,Z]} {it:E[D|X]} {it:E[D|X,Z]} {it:E[Z|X]}.
{it:command} is a supported supervised ML program (e.g. {helpb pystacked} or {helpb cvlasso}). 
See {helpb ddml##compatibility:supported programs}.

{pstd}
Note: Options before ":" and after the first comma refer to {cmd:ddml}. 
Options that come after the final comma refer to the estimation command. 
{p_end}

{pstd}
{ul:Step 3.} Cross-fitting:

{p 8 14}{cmd:ddml crossfit} [ , {opt mname(name)} {opt shortstack}{bind: ]} 

{pstd}
This step implements the cross-fitting algorithm. Each learner is fitted iteratively on training folds and out-of-sample predicted values are obtained.

{pstd}
{ul:Step 4.} Estimate causal effects:

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)} {opt atet} {opt ateu} {opt trim(real)}{bind: ]} 

{pstd}
The {cmd:ddml estimate} command returns treatment effect estimates for all combination of learners 
added in Step 2.

{pstd}
{ul:Optional.} Report/post selected results:

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {opt spec(integer or string)} {opt rep(integer or string)} {opt allcombos} {opt not:able} {opt replay} {bind: ]} 

{pstd}
{marker auxiliary}{...}
{ul:Optional.} Retrieve information from {cmd:ddml}:

{p 8 14}{cmd:ddml extract} [ {it:object_name} , {opt mname(name)} {opt show(display_item)} {opt ename(name)} {opt vname(varname)}
{opt stata} {opt keys} {opt key1(string)} {opt key2(string)} {opt key3(string)} {opt subkey1(string)}
{opt subkey2(string)}{bind: ]}

{pstd}
{it:display_item} can be {it:stweights}, {it:ssweights}, {it:psweights}, {it:weights}, {it:mse}, {it:n}, or {it:pystacked}.
{cmd:ddml} stores many internal results on associative arrays.
See {helpb ddml extract} for details.

{pstd}
For full details and further options, follow the links to the detailed help files {helpb ddml##help:above}.


{marker compatibility}{...}
{title:Compatible programs}

{pstd}
{marker general}{...}
{opt ddml} is compatible with a large set of user-written Stata commands. 
It has been tested with 

{p 7 9 0} 
- the {helpb pystacked} package (see {helpb pystacked} and {help ddml##pystacked:below}). 
Note that {helpb pystacked} requires Stata 16.

{p 7 9 0} 
- {helpb lassopack} for regularized regression (see {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}).

{p 7 9 0} 
- {helpb rforest} by Zou & Schonlau. Note that {cmd:rforest} requires the option 
{cmd:vtype(none)}. 

{p 7 9 0} 
- {helpb svmachines} by Guenther & Schonlau.

{pstd}
Beyond these, it is compatible with any Stata program that 

{p 7 9 0} 
- uses the standard "{it:reg y x}" syntax,

{p 7 9 0} 
- supports {it:if}-conditions,

{p 7 9 0} 
- and comes with {helpb predict} post-estimation programs.

{marker pystacked}{...}
{pstd}
{helpb pystacked} implements stacking regression ({helpb pystacked##Wolpert1992:Wolpert, 1992})
via Stata's Python integration in combination with
{browse "https://scikit-learn.org/stable/index.html":scikit-learn}'s 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingRegressor.html":sklearn.ensemble.StackingRegressor} and 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingClassifier.html":sklearn.ensemble.StackingClassifier}. 
Stacking is a way of combining multiple supervised
machine learners (the "base" or "level-0" learners) into
a meta learner.

{pstd}{helpb pystacked} is the recommended way to specify multiple learners in {opt ddml}.
{helpb pystacked} provides a fast way of estimating all the learners in a single call to one program,
and {opt ddml} has integrated support for various features provided by {helpb pystacked}.
{opt ddml} will store the predicted values of the specified base learners as well as the combined ("stacked") predicted values.
It also stores the stacking weights used by {help pystacked} along with the {opt ddml} short-stacking weights.

{pstd}{bf:Important}: For these features to be available, {helpb pystacked} needs to be the only learner for each conditional expectation.
Multiple learners must be specified in the call to {helpb pystacked}; see the examples below.
{helpb pystacked} can be provided directly to {opt ddml} as one of several learners for a conditional expectation,
but in this case the extra features for {helpb pystacked} will not be availabe.


{marker installation}{title:Installation}

{pstd}
To get the latest stable version of {cmd:ddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.
Alternatively, just use the {opt update} subcommand or click here: {stata "ddml update"}.

{pstd}
To verify that {cmd:ddml} is correctly installed, 
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

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}


{marker authors}{title:Authors}

{pstd}
Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland  {break}
achim.ahrens@gess.ethz.ch

{pstd}
Christian B. Hansen, University of Chicago, USA {break}
Christian.Hansen@chicagobooth.edu

{pstd}
Mark E. Schaffer, Heriot-Watt University, UK {break}
m.e.schaffer@hw.ac.uk	

{pstd}
Thomas Wiemann, University of Chicago, USA {break}
wiemann@uchicago.edu


{title:Also see (if installed)}

{pstd}
Help: {helpb pystacked}, {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}.{p_end}
