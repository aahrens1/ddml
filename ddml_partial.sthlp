{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help ddml partial}{right: v1.2}
{hline}

{title:ddml - estimation examples for the partially-linear model in Double Debiased Machine Learning}

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
The {opt ddml} package also includes the wrapper program {helpb qddml},
which uses a simplified one-line syntax, but offers less flexibility.

{pstd}
This help file illustrates usage of the {ul:partially-linear model}.
For examples of other models,
follow the links in the main {help ddml:ddml help file}.

{pstd}
We use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Partially-linear model} [{it:partial}]

	Y = {it:a}.D + g(X) + U
        D = m(X) + V

{pstd}
where the aim is to estimate {it:a} while controlling for X. To this end, 
we estimate the conditional expectations
E[Y|X] and E[D|X] using a supervised machine learner.


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:ddml} for the partially-linear model. 
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

{pstd}Because we have used {help pystacked} as the single {opt ddml} learner,
we can access the saved {opt pystacked} information.
Here we use the {opt pystacked} option to get the stacking weights and MSEs by cross-fit fold:{p_end}
{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}


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
