{smcl}
{* *! version 24jul2023}{...}
{hline}
{cmd:help ddml partial}{right: v1.2}
{hline}

{title:ddml - estimation of the partially-linear model in Double Debiased Machine Learning}

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
{pstd}{help ddml_example_partial_pystacked_basic:Basic example of the partially-linear model with pystacked}

{marker plm_iii}{...}
{pstd}{help ddml_example_partial_pystacked_detailed:Detailed example of the partially-linear model with pystacked}

{marker plm_iii}{...}
{pstd}{help ddml_example_partial_anylearner_detailed:Detailed general example of the partially-linear model with any learner(s)}

{marker plm_iii}{...}
{pstd}{help ddml_example_partial_pystacked_multitreat:Estimating the partially-linear model with multiple treatments}


{smcl}
INCLUDE help ddml_install_ref_auth

