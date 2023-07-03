{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help ddml partial}{right: v1.2}
{hline}

{title:ddml - estimation examples for the interactive (ATE, ATET) model in Double Debiased Machine Learning}

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
The {opt ddml} package also includes the wrapper program {helpb qddml},
which uses a simplified one-line syntax, but offers less flexibility.

{pstd}
This help file illustrates usage of the {ul:interactive model}
used to obtain estimates of the ATE (average treatment effect)
and ATET (average treatment effect on the treated).
For examples of other models,
follow the links in the main {help ddml:ddml help file}.

{pstd}
We use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Interactive model} [{it:interactive}]

	Y = g(X,D) + U
        D = m(X) + V

{pstd}
which (compared to the {help ddml partial:partially-linear model}
relaxes the assumption that X and D are separable. 
D is a binary treatment variable. 
We estimate, using a supervised machine
learner, the following conditional expectations:
{p_end}
{phang2}1. E[Y|X,D=0] and E[Y|X,D=1], jointly added using {cmd:ddml E[Y|X,D]}{p_end}
{phang2}2. E[D|X], added using {cmd:ddml E[D|X]}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:ddml} for the interactive model. 
Note that estimation models are chosen for demonstration purposes only and 
kept simple to allow you to run the code quickly.

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
the average treatment effect on the treated ({opt atet}),
or the average treatment effect on the untreated ({opt ateu}).{p_end}
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
