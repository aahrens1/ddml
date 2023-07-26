{smcl}
{* *! version 26july2023}{...}
{viewerjumpto "Examples" "ddml_iv##examples"}{...}
{viewerjumpto "Installation" "ddml_iv##installation"}{...}
{viewerjumpto "References" "ddml_iv##references"}{...}
{viewerjumpto "Authors" "ddml_iv##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_iv##also_see"}{...}
{hline}
{cmd:help ddml iv, help ddml fiv}{right: v1.4}
{hline}

{title:ddml - estimation of partially-linear IV models in Double Debiased Machine Learning}

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
may be kept simple to allow you to run the code quickly.

{pstd}{help ddml iv##iv_pystacked_basic:1. Basic example of the partially-linear IV model with pystacked}{p_end}
{pstd}{help ddml iv##iv_anylearner_basic:2. Basic example of the partially-linear IV model with any learner(s)}{p_end}
{pstd}{help ddml iv##fiv_anylearner_basic:3. Basic example of the flexible partially-linear IV model with any learner(s)}{p_end}
{pstd}{help ddml iv##fiv_anylearner_detailed:4. Detailed example of the flexible partially-linear IV model with any learner(s)}{p_end}

{marker iv_pystacked_basic}
{smcl}
INCLUDE help ddml_example_partialiv_pystacked_basic.sthlp

{marker iv_anylearner_basic}
{smcl}
INCLUDE help ddml_example_partialiv_anylearner_basic.sthlp

{marker fiv_anylearner_basic}
{smcl}
INCLUDE help ddml_example_flexiv_anylearner_basic.sthlp

{marker fiv_anylearner_detailed}
{smcl}
INCLUDE help ddml_example_flexiv_anylearner_detailed.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
