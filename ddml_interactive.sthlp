{smcl}
{* *! version 24jul2023}{...}
{hline}
{cmd:help ddml interactive}{right: v1.2}
{hline}

{title:ddml - estimation of the interactive (ATE, ATET) model in Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables.

{pstd}
{opt ddml} supports a variety of different ML programs, including
but not limited to {help pystacked} and {help lassopack}. 
{help pystacked} is the recommended way to specify multiple learners in {opt ddml},
and {opt ddml} has integrated support for various features provided by {help pystacked}.

{pstd}
The {opt ddml} package also includes the wrapper program {help qddml},
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

{marker interactive_i}{...}
{pstd}{help ddml_example_interactive_pystacked_basic:Basic example of the interactive model (ATE, ATET) with pystacked}

{marker interactive_ii}{...}
{pstd}{help ddml_example_interactive_pystacked_detailed:Detailed example of the interactive model (ATE, ATET) with pystacked}


{smcl}
INCLUDE help ddml_install_ref_auth
