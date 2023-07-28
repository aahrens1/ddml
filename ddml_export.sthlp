{smcl}
{* *! version 28july2023}{...}
{viewerjumpto "Syntax" "ddml_export##syntax"}{...}
{viewerjumpto "Examples" "ddml_export##examples"}{...}
{viewerjumpto "Installation" "ddml_export##installation"}{...}
{viewerjumpto "References" "ddml_export##references"}{...}
{viewerjumpto "Authors" "ddml_export##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_export##also_see"}{...}
{hline}
{cmd:help ddml export}{right: v1.4.1}
{hline}

{title:ddml export utility for Double Debiased Machine Learning}

{p2colset 5 19 21 2}{...}
{p2col:{hi: ddml} {hline 2}}Stata package for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml export} saves the estimated conditional expectations, cross-fold identifers, etc.
to a CSV file.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml export}
[using {it:filename} , {opt mname(name)}
{opt addvars(varlist)}

{synoptset 20}{...}
{synopthdr:options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt addvars(varlist)}}
additional Stata variables to include with {opt ddml} variables.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{marker examples}{...}
{title:Examples}

{smcl}
INCLUDE help ddml_example_export.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
