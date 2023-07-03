{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help ddml iv, help ddml fiv}{right: v1.2}
{hline}

{title:ddml - estimation examples for the partially-linear IV models in Double Debiased Machine Learning}

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
This help file illustrates usage of the {ul:partially-linear IV model}
and the {ul:flexible partially-linear IV model}.
For examples of other models,
follow the links in the main {help ddml:ddml help file}.

{pstd}
We use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, 
{it:Z} to denote instrumental variable(s), and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Partially-linear IV model} [{it:iv}]

	Y = {it:a}.D + g(X) + U
        Z = m(X) + V

{pstd}
where the aim is to estimate {it:a}. 
We estimate the conditional expectations E[Y|X], 
E[D|X] and E[Z|X] using supervised machine learners.
Note that the instrument set Z is low-dimensional.

{pstd}
{ul:Flexible partially-linear IV model} [{it:fiv}]

	Y = {it:a}.D + g(X) + U
        D = m(Z) + g(X) + V 

{pstd}
where the estimand of interest is {it:a}. 
We estimate the conditional expectations
E[Y|X], 
E[D^|X] and D^:=E[D|Z,X] using supervised machine learnerd.
The instrument is then formed as D^-E^[D^|X] where E^[D^|X] denotes
the estimate of E[D^|X]. 

{pstd}
Note: "{D}" is a placeholder that is used because last step (estimation of E[D|X]) 
uses the fitted values from estimating E[D|X,Z].

{pstd}
{ul:Which IV model?}

{pstd}
The flexible partially-linear IV Model allows for approximation of optimal instruments
as in Belloni et al. ({help ddml iv##BCCH2012:2012}),
but relies on a stronger independence assumption than the partially-linear IV Model.
Specifically, the partially-linear IV model uses an orthogonality condition,

	E[Cov(U,Z|X)] = 0

{pstd}
whereas the flexible partially-linear IV model uses the conditional mean independence condition

	E[U|Z,X] = 0

{pstd}
Note that the generated instruments are generally valid,
Also note that (unlike the standard partially-linear IV model above),
the flexible partially-linear IV model can accommodate both low- and high-dimensional instrument sets Z.


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:ddml} for partially-linear IV models.
Note that estimation models are chosen for demonstration purposes only and 
kept simple to allow you to run the code quickly.


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


{marker references}{title:References}

{marker BCCH2012}{...}
{pstd}
Belloni, A., Chen, D., Chernozhukov, V. and Hansen, C. 2012.
Sparse models and methods for optimal instruments with an application to eminent domain.
{it:Econometrica} 80(6):2369-2429.
{browse "http://onlinelibrary.wiley.com/doi/10.3982/ECTA9626/abstract"}

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
