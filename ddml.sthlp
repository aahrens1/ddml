{smcl}
{* *! version 21feb2023}{...}
{hline}
{cmd:help ddml}{right: v1.2}
{hline}

{title:ddml - Stata package for Double Debiased Machine Learning}

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
{p 2}{help ddml##examples:Examples}{p_end}
{p 4}{help ddml##plm_i:Partially linear model I - multiple learners with pystacked}{p_end}
{p 4}{help ddml##plm_ii:Partially linear model II - multiple learners with ddml}{p_end}
{p 4}{help ddml##plm_iii:Partially linear model III - multiple treatments}{p_end}
{p 4}{help ddml##pliv:Partially linear IV model}{p_end}
{p 4}{help ddml##interactive:Interactive model - ATE and ATET estimation}{p_end}
{p 4}{help ddml##late:Interactive IV model - LATE estimation}{p_end}
{p 4}{help ddml##fiv:Flexible Partially Linear IV model}{p_end}
{p 2}{help ddml##installation:Installation}{p_end}
{p 2}{help ddml##references:References}{p_end}


{marker help}{...}
{title:Follow links below to detailed help}

{p2colset 5 20 20 0}{...}
{p 2}Help files: qddml{p_end}
{p2col:{help qddml}} One-step DDML estimation using {help qddml}.{p_end}

{p 2}Help files: main steps in using {opt ddml}{p_end}
{p2col:{help ddml overview}} Overview of main steps to estimation using ddml{p_end}
{p2col:{help ddml stacking}} Overview of stacking in ddml and {help pystacked}{p_end}

{p 2}Help files: detailed help files for main steps in using {opt ddml}{p_end}
{p2col:{help ddml init}} 1. Initialize {opt ddml} and select model.{p_end}
{p2col:{help ddml eq}} 2. Add supervised ML programs for estimating conditional expectations.{p_end}
{p2col:{help ddml crossfit}} 3. Cross-fitting to estimate conditional expectations.{p_end}
{p2col:{help ddml estimate}} 4. Estimate causal model and report/post results.{p_end}

{p 2}Help files: auxiliary programs{p_end}
{p2col:{help ddml describe}} Report information about the model setup and/or results.{p_end}
{p2col:{help ddml extract}} Report information about saved results e.g. stacking weigths.{p_end}
{p2col:{help ddml sample}} Report information about the estimation sample, folds, etc.{p_end}
{p2col:{help ddml export}} Save the {opt ddml} estimated conditional expectations to a csv file.{p_end}
{p2col:{help ddml overlap}} (interactive models only) Generate overlap plots for propensity-score-based models{p_end}
{p2col:{help crossfit}} Use {opt crossfit} as a standalone program for cross-fitting and cross-validation.{p_end}


{marker overview}{...}
{title:Overview of main step when estimating with ddml}

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


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:ddml} for each of the 5 models supported. 
Note that estimation models are chosen for demonstration purposes only and 
kept simple to allow you to run the code quickly.

{marker plm_i}{...}
{pstd}{ul:Partially linear model I. Stacking regression using {helpb pystacked}.}

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We next initialize the ddml estimation and select the model.
{it:partial} refers to the partially linear model.
The model will be stored on a Mata object with the default name "m0"
unless otherwise specified using the {opt mname(name)} option.{p_end}

{pstd}We set the number of random folds to 2 so that 
the model runs quickly. The default is {opt kfolds(5)}. We recommend 
considering at least 5-10 folds and even more if your sample size is small.{p_end}

{pstd}We recommend re-running the model multiple times on 
different random folds; see options {opt reps(integer)}.
Here we set the number of repetions to 2, again only so that the model runs quickly.{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}

{pstd}Stacking regression is a simple and powerful method for 
combining predictions from multiple learners.
Here we use {helpb pystacked} with the partially linear model,
but it can be used with any model supported by {cmd:ddml}.{p_end}

{pstd}Note: the additional support provided by {opt ddml} for {helpb pystacked} (see {help ddml##pystacked:above})
is available only if, as in this example, {helpb pystacked} is the only learner for each conditional expectation.
Mutliple learners are provided to {helpb pystacked}, not directly to {opt ddml}.

{pstd}Add supervised machine learners for estimating conditional expectations.
The first learner in the stacked ensemble is OLS.
We also use cross-validated lasso, ridge and two random forests with different settings, 
which we save in the following macros:{p_end}
{phang2}. {stata "global rflow max_features(5) min_samples_leaf(1) max_samples(.7)"}{p_end}
{phang2}. {stata "global rfhigh max_features(5) min_samples_leaf(10) max_samples(.7)"}{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X || method(ols) || method(lassocv) || method(ridgecv) || method(rf) opt($rflow) || method(rf) opt($rfhigh), type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X || method(ols) || method(lassocv) || method(ridgecv) || method(rf) opt($rflow) || method(rf) opt($rfhigh), type(reg)"}{p_end}

{pstd}Note: Options before ":" and after the first comma refer to {cmd:ddml}. 
Options that come after the final comma refer to the estimation command. 
Make sure to not confuse the two types of options.{p_end}

{pstd}Check if learners were correctly added:{p_end}
{phang2}. {stata "ddml desc, learners"}{p_end}

{pstd} Cross-fitting: The learners are iteratively fitted on the training data.
This step may take a while, depending on the number of learners, repetitions, folds, etc.
In addition to the standard stacking done by {helpb pystacked},
also request short-stacking to be done by {opt ddml}.
Whereas stacking relies on (out-of-sample) cross-validated predicted values
to obtain the relative weights for the base learners,
short-stacking uses the (out-of-sample) cross-fitted predicted values.{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}

{pstd}Finally, we estimate the coefficients of interest.{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Examine the standard ({cmd:pystacked}) stacking weights as well as the {opt ddml} short-stacking weights.{p_end}
{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}

{marker plm_ii}{...}
{pstd}{ul:Partially linear model II - multiple learners with ddml.} 

{pstd}Here we used {opt ddml} to add learners. This allows use of learners not supported by,
or as alternatives to, those available via {helpb pystacked}.
It is also possible to use {helpb pystacked} as a standalone learner in this way.{p_end}

{pstd}Preparation: use the data and globals as above.
Use the name {cmd:m1} for this new estimation, 
to distinguish it from the previous example that uses the default name {cmd:m0}.
This enables having multiple estimations available for comparison.
Also specify 5 resamplings.{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(5) mname(m1)"}{p_end}

{pstd}We add supervised machine learners for estimating the conditional 
expectation E[Y|X].
In each step, we add the {cmd:mname(m1)} option to ensure that the learners
are not added to the {cmd:m0} model which is still in memory.
We also specify the names of the variables containing the estimated conditional
expectations using the {opt learner(varname)} option.
This avoids overwriting the variables created for the {cmd:m0} model using default naming.{p_end}

{pstd} We first add simple linear regression.{p_end}
{phang2}. {stata "ddml E[Y|X], mname(m1) learner(Y_reg_m1): reg $Y $X"}{p_end}

{pstd}We can add more than one learner per reduced form equation. Here, we 
add a random forest learner. We do this using {helpb pystacked} to implement a single learner.{p_end}
{phang2}. {stata "ddml E[Y|X], mname(m1) learner(Y_pys_m1): pystacked $Y $X, type(reg) method(rf)"}{p_end}

{pstd}We do the same for the conditional expectation E[D|X].{p_end}
{phang2}. {stata "ddml E[D|X], mname(m1) learner(D_reg_m1): reg $D $X"}{p_end}
{phang2}. {stata "ddml E[D|X], mname(m1) learner(D_pys_m1): pystacked $D $X, type(reg) method(rf)"}{p_end}

{pstd}Check if learners were correctly added:{p_end}
{phang2}. {stata "ddml desc, mname(m1) learners"}{p_end}

{pstd}Cross-fitting and estimation.
Since we added two learners for each of our two reduced form equations, 
there are four possible specifications. 
By default, the result shown corresponds to the specification 
with the lowest out-of-sample MSPE:
{p_end}
{phang2}. {stata "ddml crossfit, mname(m1)"}{p_end}
{phang2}. {stata "ddml estimate, mname(m1) robust"}{p_end}

{pstd}To estimate all four specifications, we use the {cmd:allcombos} option:
{p_end}
{phang2}. {stata "ddml estimate, mname(m1) robust allcombos"}{p_end}

{pstd}After having estimated all specifications, we can retrieve 
specific results. Here we use the specification relying on OLS for both
estimating both E[Y|X] and E[D|X], from the 4th cross-fit split ({opt rep(4)}.
(Note: Working interactively, the simplest way to do this
is to click on the hyperlink in the summary table in the {opt ddml estimate} output above.)
The {opt notable} option suppresses the summary table:
{p_end}
{phang2}. {stata "ddml estimate, mname(m1) spec(1) rep(4) replay notable"}{p_end}

{pstd}You could manually retrieve the same point estimate by 
typing:
{p_end}
{phang2}. {stata "reg Y_reg_m1_4 D_reg_m1_4, robust"}{p_end}
{pstd}or graphically:
{p_end}
{phang2}. {stata "twoway (scatter Y_reg_m1_4 D_reg_m1_4) (lfit Y_reg_m1_4 D_reg_m1_4)"}{p_end}

{pstd}where {opt Y_reg_m1_4} and {opt D_reg_m1_4} are the orthogonalized
versions of {opt net_tfa} and {opt e401} from the 4th cross-fit estimation.
{p_end}

{pstd}To describe the ddml model setup or results in detail,
you can use {cmd: ddml describe} with the relevant option ({opt sample}, {opt learners}, {opt crossfit}, {opt estimates}),
or just describe them all with the {opt all} option:
{p_end}
{phang2}. {stata "ddml describe, mname(m1) all"}{p_end}

{pstd}We can compare the effects with the first {cmd:ddml} model 
(if you have run the first example above).{p_end}
{phang2}. {stata "ddml estimate, mname(m0) replay"}{p_end}

{marker plm_iii}{...}
{pstd}{ul:Partially linear model III. Multiple treatments.}

{pstd}We can also run the partially linear model with multiple treatments. 
In this simple example, we estimate the effect of both 401k elligibility 
{cmd:e401} and education {cmd:educ}. 
Note that we remove {cmd:educ} from the set of controls.
We again use {helpb pystacked} as the single learner provided to {opt ddml};
the two base learners, OLS and random forest, are provided via {opt pystacked}.
This time we use the alternative simplified syntax supported by {helpb pystacked}.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D1 e401"}{p_end}
{phang2}. {stata "global D2 educ"}{p_end}
{phang2}. {stata "global X tw age inc fsize db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Initialize the model.{p_end}
{phang2}. {stata "ddml init partial, kfolds(2)"}{p_end}

{pstd}Add learners. Note that we add leaners with both {cmd:$D1} and
{cmd:$D2} as the dependent variable.{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg) methods(ols rf)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D1 $X, type(reg) methods(ols rf)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D2 $X, type(reg) methods(ols rf)"}{p_end}

{pstd}Cross-fitting.{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}

{pstd}Estimation.{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Because we have used {pystacked as the single {opt ddml} learner,
we can access the saved {opt pystacked} information.
Here we use the {opt pystacked} option to get the stacking weights and MSEs by cross-fit fold:{p_end}
{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}

{marker pliv}{...}
{pstd}{ul:Partially linear IV model.} 

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}
{phang2}. {stata "use https://statalasso.github.io/dta/AJR.dta, clear"}{p_end}
{phang2}. {stata "global Y logpgp95"}{p_end}
{phang2}. {stata "global D avexpr"}{p_end}
{phang2}. {stata "global Z logem4"}{p_end}
{phang2}. {stata "global X lat_abst edes1975 avelf temp* humid* steplow-oilres"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Preparation: we load the data, define global macros and set the seed. Since the
data set is very small, we consider 30 cross-fitting folds.{p_end}
{phang2}. {stata "ddml init iv, kfolds(30)"}{p_end}

{pstd}The partially linear IV model has three conditional expectations: 
E[Y|X], E[D|X] and E[Z|X]. For each reduced form equation, we use two learners:
OLS and random forest.
To illustrate how {opt ddml} works with other packages,
instead of a single call to {opt pystacked} specifying two base learners
we specify Stata's {helpb regress} and {helpb rforest} by Zou and Schonlau as the two learners.
We need to add the option {opt vtype(none)} for {helpb rforest} to 
work with {cmd:ddml} since {helpb rforest}'s {cmd:predict} command doesn't
support variable types.{p_end}
{phang2}. {stata "ddml E[Y|X]: reg $Y $X"}{p_end}
{phang2}. {stata "ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: reg $D $X"}{p_end}
{phang2}. {stata "ddml E[D|X], vtype(none): rforest $D $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: reg $Z $X"}{p_end}
{phang2}. {stata "ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)"}{p_end}

{pstd}Cross-fitting and estimation.{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}If you are curious about what {cmd:ddml} does in the background:{p_end}
{phang2}. {stata "ddml estimate, allcombos robust"}{p_end}
{phang2}. {stata "ddml estimate, spec(8) rep(1) replay notable"}{p_end}
{phang2}. {stata "ivreg Y2_rf (D2_rf = Z2_rf), robust"}{p_end}

{marker interactive}{...}
{pstd}{ul:Interactive model--ATE and ATET estimation.} 

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}
{phang2}. {stata "webuse cattaneo2, clear"}{p_end}
{phang2}. {stata "global Y bweight"}{p_end}
{phang2}. {stata "global D mbsmoke"}{p_end}
{phang2}. {stata "global X prenatal1 mmarried fbaby mage medu"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We use 5 folds and 5 resamplings; that is, 
we estimate the model 5 times using randomly chosen folds.{p_end}
{phang2}. {stata "ddml init interactive, kfolds(5) reps(5)"}{p_end}

{pstd}We need to estimate the conditional expectations of E[Y|X,D=0], 
E[Y|X,D=1] and E[D|X]. The first two conditional expectations 
are added jointly.{p_end} 
{pstd}We consider two supervised learners: linear regression and gradient boosted
trees, stacked using {helpb pystacked}.
Note that we use gradient boosted regression trees for E[Y|X,D], but
gradient boosted classification trees for E[D|X].
{p_end} 
{phang2}. {stata "ddml E[Y|X,D]: pystacked $Y $X, type(reg) methods(ols gradboost)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X, type(class) methods(logit gradboost)"}{p_end}

{pstd}Cross-fitting and short-stacking:{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}

{pstd}In the final estimation step, we can estimate
the average treatment effect (the default),
the average treatment effect of the treated ({opt atet}),
or the average treatment effect of the untreated ({opt ateu}).{p_end}
{phang2}. {stata "ddml estimate"}{p_end}
{phang2}. {stata "ddml estimate, atet"}{p_end}

{pstd}Recall that we have specified 5 resampling iterations ({opt reps(5)})
By default, the median over short-stacked resampling iterations is shown.
At the bottom, a table of summary statistics over resampling iterations is shown. 
To display the mean over standard stacking results, i.e.,
the results where the weights derive from {helpb pystacked} and vary by cross-fit fold,
we use [opt ddml estimate, replay} with {opt spec(mse) and {opt rep(mn)}
(because {opt pystacked} is the only learner, it is also the "minimum MSE learner").{p_end}
{phang2}. {stata "ddml estimate, spec(mse) rep(mn) replay notable"}{p_end}

{marker late}{...}
{pstd}{ul:Interactive IV model--LATE estimation.} 

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}
{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We initialize the model.{p_end}
{phang2}. {stata "ddml init interactiveiv, kfolds(5)"}{p_end}

{pstd}We use {helpb pystacked} with two base learners for each reduced form equation.{p_end}
{phang2}. {stata "ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(ols lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}

{pstd}Cross-fitting and estimation, with short-stacking implemented via {opt ddml}.{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Compare the short-stacking estimation above with standard (within-cross-fit-fold) stacking:{p_end}
{phang2}. {stata "ddml estimate, spec(mse) rep(1) replay notable"}{p_end}

{marker fiv}{...}
{pstd}{ul:Flexible Partially Linear IV model.} 

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear"}{p_end}
{phang2}. {stata "global Y share"}{p_end}
{phang2}. {stata "global D price"}{p_end}
{phang2}. {stata "global X hpwt air mpd space"}{p_end}
{phang2}. {stata "global Z sum*"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We initialize the model.{p_end}
{phang2}. {stata "ddml init fiv"}{p_end}

{pstd}We add learners for E[Y|X] in the usual way. Here we use {helpb pystacked}'s default base learners.{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg)"}{p_end}

{pstd}There are some pecularities that we need to bear in mind
when adding learners for E[D|Z,X] and E[D|X].
The reason for this is that the estimation of E[D|X]
depends on the estimation of E[D|X,Z].
More precisely, we first obtain the fitted values D^=E[D|X,Z] and 
fit these against X to estimate E[D^|X].{p_end}

{pstd}
When adding learners for E[D|Z,X],
we need to provide a name
for each learners using {opt learner(name)}.{p_end}
{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)"}{p_end}

{pstd}
When adding learners for E[D|X], we explicitly refer to the learner from 
the previous step (e.g., {cmd:learner(Dhat_pystacked)}) and
also provide the name of the treatment variable ({cmd:vname($D)}).
Finally, we use the placeholder {cmd:{D}} in place of the dependent variable. 
{p_end}
{phang2}. {stata "ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)"}{p_end}
 
{pstd}That's it. Now we can move to cross-fitting and estimation.{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}If you are curious about what {cmd:ddml} does in the background:{p_end}
{phang2}. {stata "gen Dtilde = $D - Dhat_pystacked_h_1"}{p_end}
{phang2}. {stata "gen Zopt = Dhat_pystacked_1 - Dhat_pystacked_h_1"}{p_end}
{phang2}. {stata "ivreg Y1_pystacked_1 (Dtilde=Zopt), robust"}{p_end}


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
