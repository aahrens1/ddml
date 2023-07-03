{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help ddml init, ddml eq, ddml sample}{right: v1.2}
{hline}

{title:ddml init, eq and sample commands for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml init} {it:model} initializes the model,
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, or {it:interactiveiv}.

{pstd}
{cmd: ddml eq: command} adds supervised ML programs for estimating conditional expectations,
where {it:eq} is the conditional expectation to be estimated (e.g., {it:E[Y|X]})
and {it:command} is a supported supervised ML program.

{pstd}
{opt ddml sample} adds cross-fitting repetitions to an existing and possibly already-estimated model.


{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml init}
{it:model} [if] [in]
[ , {opt mname(name)} {opt kfolds(integer)}
{opt fcluster(varname)}
{opt foldvar(varlist)} {opt reps(integer)} 
{opt norandom} {opt tabfold} {opt vars(varlist)}{bind: ]}

{pstd}
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, {it:interactiveiv}.
See {help ddml_init##models:Models} below.

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
Options that come after ":" and the final comma refer to the estimation command. 
{p_end}

{p 8 14}{cmd:ddml sample} [ , {opt append}[{cmd:(}{it:integer}{cmd:)}] {opt foldvar(varlist)} {bind: ]}

{pstd}
adds cross-fitting repetitions to an existing and possibly already-estimated model,
where the additional repetitions is indicated either by {opt append(#)}
or by {opt append} and the cross-fit fold identifiers in {opt foldvar(varlist)}.


{synoptset 20}{...}
{synopthdr:init options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt kfolds(integer)}}
number of cross-fitting folds. The default is 5.
{p_end}
{synopt:{opt fcluster(varname)}}
cluster identifiers for cluster randomization of random folds.
{p_end}
{synopt:{opt foldvar(varlist)}}
integer variable with user-specified cross-fitting folds (one per cross-fitting repetition).
{p_end}
{synopt:{opt norandom}}
use observations in existing order instead of randomizing before splitting into folds;
if multiple resamples, applies to first resample only;
ignored if user-defined fold variables are provided in {opt foldvar(varlist)}.
{p_end}
{synopt:{opt reps(integer)}}
cross-fitting repetitions, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds. 
{p_end}
{synopt:{opt tabfold}}
prints a table with frequency of observations by fold.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:equation options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt vname(varname)}}
name of the dependent variable in the reduced form estimation. 
This is usually inferred from the command line but is mandatory
for the {it:fiv} model.
{p_end}
{synopt:{opt l:earner(varname)}}
optional name of the variable to be created. 
{p_end}
{synopt:{opt vtype(string)}}
(rarely used) optional variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank 
(required when using {cmd:ddml} with {helpb rforest}.)
{p_end}
{synopt:{opt predopt(string)}}
(rarely used) {cmd:predict} option to be used to get predicted values. 
Typical values could be {opt xb} or {opt pr}. Default is 
blank. 
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:sample options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt append(#)}}
number of additional resamples to cross-fit.
{p_end}
{synopt:{opt append}}
when no number of resamples to append is provided, this is based on the list fold IDs in {opt foldvar(varlist)}.
{p_end}
{synopt:{opt foldvar(varlist)}}
integer variable with user-specified cross-fitting folds (one per cross-fitting repetition).
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{marker models}{...}
{title:Models}

{pstd}
Throughout we use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, 
{it:Z} to denote instrumental variable(s), and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Partially-linear model} [{it:partial}]

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
D is a binary treatment variable,
and we aim to estimate the average treatment effect (ATE)
or average treatment effect on the treated (ATET).
We estimate the conditional expectations E[D|X], as well as 
E[Y|X,D=0] and E[Y|X,D=1] (jointly added using {cmd:ddml E[Y|X,D]}).

{pstd}
{ul:Partially-linear IV model} [{it:iv}]

	Y = {it:a}.D + g(X) + U
        Z = m(X) + V

{pstd}
where the aim is to estimate {it:a}. 
We estimate the conditional expectations E[Y|X], 
E[D|X] and E[Z|X] using a supervised machine
learner.

{pstd}
{ul:Interactive IV model}  [{it:interactiveiv}]

	Y = g(Z,X) + U
        D = h(Z,X) + V
        Z = m(X) + E

{pstd}
where the aim is to estimate the local average treatment effect (LATE).
We estimate, using a supervised machine
learner, the following conditional expectations:
E[Y|X,Z=0] and E[Y|X,Z=1] (jointly added using {cmd:ddml E[Y|X,Z]});
E[D|X,Z=0] and E[D|X,Z=1] (jointly added using {cmd:ddml E[D|X,Z]});
E[Z|X].

{pstd}
{ul:Flexible partially-linear IV model} [{it:fiv}]

	Y = {it:a}.D + g(X) + U
        D = m(Z) + g(X) + V 

{pstd}
where the estimand of interest is {it:a}. 
We estimate the conditional expectations
E[Y|X], 
E[D^|X] and D^:=E[D|Z,X] using a supervised machine
learner. The instrument is then formed as D^-E^[D^|X] where E^[D^|X] denotes
the estimate of E[D^|X]. 

{pstd}
Note: "{D}" is a placeholder that is used because last step (estimation of E[D|X]) 
uses the fitted values from estimating E[D|X,Z].


{title:Examples}

{pstd}
For more examples of usage see the links via the main {help ddml##examples:ddml help file}.
See {help ddml estimate:help ddml estimate} for details of cross-fitting and estimation options.

{pstd}Partially-linear model: load the data, define global macros, set the seed and initialize the model.
Use 2-fold cross-fitting with two repetitions (resamples)
Use {help pystacked}'s default learners as the supervised learners: OLS, cross-validated lasso, and gradient boosting.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Interactive model with 5 folds and 2 resamplings.
We need to estimate the conditional expectations of E[Y|X,D=0], E[Y|X,D=1] and E[D|X].
The first two conditional expectations are added jointly.
Two supervised learners: linear regression and gradient boosted
trees, stacked using {help pystacked}.
We use {help pystacked}'s 2nd syntax and stack using the single-best learner
(rather than the default constrained least squares).
Note that we use gradient boosted regression trees for E[Y|X,D], but
gradient boosted classification trees for E[D|X].
{p_end}
{phang2}. {stata "webuse cattaneo2, clear"}{p_end}
{phang2}. {stata "global Y bweight"}{p_end}
{phang2}. {stata "global D mbsmoke"}{p_end}
{phang2}. {stata "global X prenatal1 mmarried fbaby mage medu"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init interactive, kfolds(5) reps(2)"}{p_end}
{phang2}. {stata "ddml E[Y|X,D]: pystacked $Y $X || method(ols) || method(gradboost) || , type(reg) finalest(singlebest)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X || method(logit) || method(gradboost) || , type(class) finalest(singlebest)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Partially linear IV model.
The model has three conditional expectations: 
E[Y|X], E[D|X] and E[Z|X]. For each reduced form equation, we use two learners:
OLS and random forest.
To illustrate how {opt ddml} works with other packages,
instead of a single call to {opt pystacked} specifying two base learners
we specify Stata's {helpb regress} and {helpb rforest} by Zou and Schonlau as the two learners.
We need to add the option {opt vtype(none)} for {helpb rforest} to 
work with {cmd:ddml} since {helpb rforest}'s {cmd:predict} command doesn't
support variable types.
Since the data set is very small, we consider 30 cross-fitting folds.{p_end}
{phang2}. {stata "use https://statalasso.github.io/dta/AJR.dta, clear"}{p_end}
{phang2}. {stata "global Y logpgp95"}{p_end}
{phang2}. {stata "global D avexpr"}{p_end}
{phang2}. {stata "global Z logem4"}{p_end}
{phang2}. {stata "global X lat_abst edes1975 avelf temp* humid* steplow-oilres"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init iv, kfolds(30)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: reg $Y $X"}{p_end}
{phang2}. {stata "ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: reg $D $X"}{p_end}
{phang2}. {stata "ddml E[D|X], vtype(none): rforest $D $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: reg $Z $X"}{p_end}
{phang2}. {stata "ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Interactive IV model--LATE estimation.
We use {help pystacked} with two base learners for each reduced form equation.{p_end}
{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init interactiveiv, kfolds(5)"}{p_end}
{phang2}. {stata "ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(ols lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Flexible partially-linear IV model: first load the data, define global macros, set the seed and initialize the model.
We add learners for E[Y|X] in the usual way. Here we use {helpb pystacked}'s default base learners.{p_end}
{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear"}{p_end}
{phang2}. {stata "global Y share"}{p_end}
{phang2}. {stata "global D price"}{p_end}
{phang2}. {stata "global X hpwt air mpd space"}{p_end}
{phang2}. {stata "global Z sum*"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init fiv"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg)"}{p_end}
{pstd}Adding learners for E[D|Z,X] and E[D|X] in the Flexible Partially-Linear IV Model is different.
The reason for this is that the estimation of E[D|X]
depends on the estimation of E[D|X,Z].
When adding learners for E[D|Z,X],
we need to provide a name for each learners using {opt learner(name)}.
When adding learners for E[D|X], we explicitly refer to the learner from 
the previous step (e.g., {cmd:learner(Dhat_pystacked)}) and
also provide the name of the treatment variable ({cmd:vname($D)}),
and we use the placeholder {cmd:{D}} in place of the dependent variable. 
{p_end}
{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}


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
Help:  {helpb pystacked}, {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}.{p_end}
