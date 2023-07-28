{smcl}
{* *! version 27july2023}{...}
{viewerjumpto "Examples" "ddml_interactiveiv##examples"}{...}
{viewerjumpto "Installation" "ddml_interactiveiv##installation"}{...}
{viewerjumpto "References" "ddml_interactiveiv##references"}{...}
{viewerjumpto "Authors" "ddml_interactiveiv##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_interactiveiv##also_see"}{...}
{hline}
{cmd:help ddml interactiveiv}{right: v1.4.1}
{hline}

{title:ddml - estimation of the interactive IV (LATE) model in Double Debiased Machine Learning}

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
may be kept simple to allow you to run the code quickly.

{pstd}{help ddml interactiveiv##pystacked_basic:1. Basic example of the interactive IV model (LATE) with pystacked}{p_end}
{pstd}{help ddml interactiveiv##pystacked_detailed:2. Detailed example of the interactive IV model (LATE) with pystacked}{p_end}

{marker pystacked_basic}
{smcl}
INCLUDE help ddml_example_interactiveiv_pystacked_basic.sthlp

{marker pystacked_detailed}
{smcl}
INCLUDE help ddml_example_interactiveiv_pystacked_detailed.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
