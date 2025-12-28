{smcl}
{* *! version 28dec2025}{...}
{viewerjumpto "Syntax" "ddml_describe##syntax"}{...}
{viewerjumpto "Examples" "ddml_describe##examples"}{...}
{viewerjumpto "Installation" "ddml_describe##installation"}{...}
{viewerjumpto "References" "ddml_describe##references"}{...}
{viewerjumpto "Authors" "ddml_describe##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_describe##also_see"}{...}
{hline}
{cmd:help ddml describe}{right: v1.5.0}
{hline}

{title:ddml describe utility for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml describe} provides information about the model setup and/or results in detail.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml describe}
[ , {opt mname(name)}
{opt sample}
{opt learn:ers}
{opt cross:fit} 
{opt est:imates}
{opt stack:ing}
{opt stweights}
{opt ssweights}
{opt psweights}
{opt cvc}
{opt mse}
{opt rmse}
{opt rsq}
{opt n}
{opt all}

{synoptset 20}{...}
{synopthdr:options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt sample}}
information about the estimation sample, folds, etc.
{p_end}
{synopt:{opt learn:ers}}
information about the different learners used to estimate conditional expectations.
{p_end}
{synopt:{opt cross:fit}}
information about results of the cross-fitting step.
{p_end}
{synopt:{opt est:imates}}
information about the estimation estimation results.
{p_end}
{synopt:{opt stack:ing}}
information about stacking weights.
{p_end}
{synopt:{opt stweights}}
information about standard stacking weights ({opt pystacked} only).
{p_end}
{synopt:{opt ssweights}}
information about short-stacking weights).
{p_end}
{synopt:{opt psweights}}
information about pooled-stacking weights ({opt pystacked} only).
{p_end}
{synopt:{opt cvc}}
results of the {helpb ddml##Lei2020:Lei (2020)} CVC (cross-validation with confidence) test.
{p_end}
{synopt:{opt mse}}
MSE by learner and resample.
{p_end}
{synopt:{opt rmse}}
RMSE by learner and resample.
{p_end}
{synopt:{opt rsq}}
R-sqs by learner and resample.
{p_end}
{synopt:{opt n}}
sample sizes by learner and resample.
{p_end}
{synopt:{opt all}}
report information relating to all options.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{marker examples}{...}
{title:Examples}

{smcl}
INCLUDE help ddml_example_describe.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
