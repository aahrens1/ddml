{smcl}
{* *! version 3jul2023}{...}
{hline}
{cmd:help ddml partial}{right: v1.2}
{hline}

{title:ddml - estimation examples for the interactive IV (LATE) model in Double Debiased Machine Learning}

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
This help file illustrates usage of the {ul:interactive IV model}
used to obtain estimates of the LATE (local average treatment effect).
For examples of other models,
follow the links in the main {help ddml:ddml help file}.

{pstd}
We use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Interactive IV model}  [{it:interactiveiv}]

	Y = g(Z,X) + U
        D = h(Z,X) + V
        Z = m(X) + E

{pstd}
where the aim is to estimate the local average treatment effect (LATE).
We estimate, using a supervised machine
learner, the following conditional expectations:
{p_end}
{phang2}1. E[Y|X,Z=0] and E[Y|X,Z=1], jointly added using {cmd:ddml E[Y|X,Z]}{p_end}
{phang2}2. E[D|X,Z=0] and E[D|X,Z=1], jointly added using {cmd:ddml E[D|X,Z]}{p_end}
{phang2}3. E[Z|X], added using {cmd:ddml E[Z|X]}{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
Below we demonstrate the use of {cmd:ddml} for the interactive IV model. 
Note that estimation models are chosen for demonstration purposes only and 
kept simple to allow you to run the code quickly.

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
