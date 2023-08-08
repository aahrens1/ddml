{smcl}
{* *! version 8aug2023}{...}
{viewerjumpto "Syntax" "ddml_describe##syntax"}{...}
{viewerjumpto "Examples" "ddml_describe##examples"}{...}
{viewerjumpto "Installation" "ddml_describe##installation"}{...}
{viewerjumpto "References" "ddml_describe##references"}{...}
{viewerjumpto "Authors" "ddml_describe##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_describe##also_see"}{...}
{hline}
{cmd:help ddml describe}{right: v1.4.2}
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
{opt learners}
{opt crossfit} 
{opt estimates}
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
{synopt:{opt learners}}
information about the differ learners used to estimate conditional expectations.
{p_end}
{synopt:{opt crossfit}}
information about results of the cross-fitting step.
{p_end}
{synopt:{opt estimates}}
information about the estimation estimation results.
{p_end}
{synopt:{opt all}}
equivalent to {opt sample} + {opt learners} + {opt crossfit} + {opt estiamtes}.
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
