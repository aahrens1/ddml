{smcl}
{* *! version 6aug2023}{...}
{viewerjumpto "Links to detailed help files" "ddml##help"}{...}
{viewerjumpto "Models" "ddml##models"}{...}
{viewerjumpto "Estimation steps" "ddml##estimation"}{...}
{viewerjumpto "Compatible programs and pystacked" "ddml##compatibility"}{...}
{viewerjumpto "Basic example" "ddml##example"}{...}
{viewerjumpto "Installation" "ddml##installation"}{...}
{viewerjumpto "References" "ddml##references"}{...}
{viewerjumpto "Authors" "ddml##authors"}{...}
{vieweralsosee "Also see" "ddml##also_see"}{...}
{hline}
{cmd:help ddml}{right: v1.4.3}
{hline}

{title:ddml - Stata package for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
({helpb ddml##Chern2018:Chernozhukov et al., Econometrics Journal, 2018}). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 
{opt ddml} also implements the stacking approaches discussed in {helpb ddml##Ahrens2024stacking:Ahrens et al. (2024a)}. 
For a companion paper, see {helpb ddml##Ahrens2024ddml:Ahrens et al. (2024b)}.

{pstd}
{opt ddml} supports a variety of different ML programs, including
but not limited to {help pystacked} and {help lassopack}. 
{help pystacked} is the recommended way to specify multiple learners in {opt ddml},
and {opt ddml} has integrated support for various features provided by {help pystacked}.

{pstd}
The {opt ddml} package includes the wrapper program {help qddml},
which uses a simplified one-line syntax, 
but offers less flexibility.

{pstd}
{opt ddml} and {opt qddml} rely on {help crossfit}, which can be used as a standalone program.


{title:Contents}

{p 2}{help ddml##help:Links to detailed help files}{p_end}
{p 2}{help ddml##models:Supported models in ddml}{p_end}
{p 2}{help ddml##estimation:Main steps in estimation using ddml}{p_end}
{p 2}{help ddml##compatibility:Compatible programs and pystacked integration}{p_end}
{p 2}{help ddml##example:Basic example: the partially-linear model with pystacked}{p_end}
{p 2}{help ddml##installation:Installation}{p_end}
{p 2}{help ddml##references:References}{p_end}
{p 2}{help ddml##authors:Authors}{p_end}


{marker help}{...}
{title:Follow links below to detailed help}

{p2colset 5 25 25 0}{...}
{p 2}Help files: qddml{p_end}
{p2col:{help qddml}} One-step DDML estimation using {help qddml}.{p_end}

{p 2}Help files: detailed help files including syntax/options for main steps in using {opt ddml}{p_end}
{p2col:{help ddml init}} 1. Initialize {opt ddml} and select model.{p_end}
{p2col:{help ddml eq}} 2. Add supervised ML programs for estimating conditional expectations.{p_end}
{p2col:{help ddml crossfit}} 3. Cross-fitting to estimate conditional expectations.{p_end}
{p2col:{help ddml estimate}} 4. Estimate causal model and report/post results.{p_end}
{p2col:{help ddml stacking}} Overview of stacking in {opt ddml} and {help pystacked}{p_end}

{marker examples}{...}
{p 2}Help files: detailed help files and examples for models supported by {opt ddml}{p_end}
{p2col:{help ddml partial}} Partially-linear model{p_end}
{p2col:{help ddml iv}} Partially-linear IV model{p_end}
{p2col:{help ddml fiv}} Flexible partially-linear IV model{p_end}
{p2col:{help ddml interactive}} Interactive model - ATE and ATET estimation{p_end}
{p2col:{help ddml interactiveiv}} Interactive IV model - LATE estimation{p_end}
{p2col:{help ddml examples}} Clickable list of all ddml examples{p_end}

{p 2}Help files: auxiliary programs{p_end}
{p2col:{help ddml describe}} Report information about the model setup and/or results.{p_end}
{p2col:{help ddml extract}} Report information about saved results e.g. stacking weigths.{p_end}
{p2col:{help ddml sample}} Add cross-fitting repetitions to an existing model.{p_end}
{p2col:{help ddml export}} Save the {opt ddml} estimated conditional expectations to a csv file.{p_end}
{p2col:{help ddml overlap}} (interactive models only) Generate overlap plots for propensity-score-based models{p_end}
{p2col:{help crossfit}} Use {opt crossfit} as a standalone program for cross-fitting and cross-validation.{p_end}

{marker overview}
{smcl}
INCLUDE help ddml_overview.sthlp


{marker example}{...}
{title:Example}

{pstd}A basic example of how to use {opt ddml} is below.
For a clickable list of all examples in the package, see {help ddml examples:help ddml examples}.{p_end}

{smcl}
INCLUDE help ddml_example_partial_pystacked_basic.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth.ihlp

