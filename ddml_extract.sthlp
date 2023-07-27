{smcl}
{* *! version 27july2023}{...}
{viewerjumpto "Syntax" "ddml_extract##syntax"}{...}
{viewerjumpto "Options" "ddml_extract##options"}{...}
{viewerjumpto "Examples" "ddml_extract##examples"}{...}
{viewerjumpto "Installation" "ddml_extract##installation"}{...}
{viewerjumpto "References" "ddml_extract##references"}{...}
{viewerjumpto "Authors" "ddml_extract##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_extract##also_see"}{...}
{hline}
{cmd:help ddml extract}{right: v1.4}
{hline}

{title:ddml extract utility for Double Debiased Machine Learning}

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

{p2colset 5 19 21 2}{...}
{p2col:{hi: ddml extract}} Stata extract utility for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
Please check the {helpb ddml extract##examples:examples} provided at the end of the help file.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml extract} [ {it:object_name} , {opt mname(name)} {opt show(display_item)} {opt ename(name)} {opt vname(varname)}
{opt stata} {opt keys} {opt key1(string)} {opt key2(string)} {opt key3(string)} {opt subkey1(string)} {opt subkey2(string)}{bind: ]}

{pstd}
{it:display_item} can be {it:mse}, {it:n} or {it:pystacked}.
{cmd:ddml} stores many internal results on associative arrays.
These can be retrieved using the different key options.

{marker options}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:main options}
{synoptline}
{synopt:{opt mname(name)}}
Name of the DDML model; a Mata object. Defaults to {it:m0}.
{p_end}
{synopt:{opt vname(name)}}
Name of a Y, D or Z variable corresponding to a DDML equation.
{p_end}
{synopt:{opt ename(name)}}
Name of a DDML equation struct; a Mata object.
Use with {helpb crossfit} or with a DDML eStruct that has been separately extracted.
{p_end}
{synopt:{opt stata}}
Saves extracted {it:object_name} in a Stata r(.) macro (default is to leave it as Mata object).
NB: does not apply to {opt show(display_item)} (see below).
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:show options}
{synoptline}
{synopt:{opt show(stweights)}}
Extracts standard stacking ({opt pystacked}) weights.
{p_end}
{synopt:{opt show(ssweights)}}
Extracts {opt shortstack} weights.
{p_end}
{synopt:{opt show(sweights)}}
Extracts {opt poolstack} weights.
{p_end}
{synopt:{opt show(weights)}}
Extracts all available weights: standard, short-stacked, pool-stacked.
{p_end}
{synopt:{opt show(pystacked)}}
Extracts detailed {opt pystacked} weights and learner MSEs, including a breakdown by cross-fit fold.
The MSEs are cross-validation MSEs and correspond to the predictions used to obtain the stacking weights;
see {helpb pystacked:help pystacked}.
{p_end}
{synopt:{opt show(mse)}}
Extracts OOS MSEs by crossfitting fold.
{p_end}
{synopt:{opt show(n)}}
Extracts sample size by crossfitting fold.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:key options}
{synoptline}
{synopt:{opt keys}}
List all keys on the relevant associative array.
{p_end}
{synopt:{opt key1(string)}}
Associative array key #1.
{p_end}
{synopt:{opt key2(string)}}
Associative array key #2.
{p_end}
{synopt:{opt key3(string)}}
Associative array key #3.
{p_end}
{synopt:{opt subkey1(string)}}
Associative array subkey #1.
{p_end}
{synopt:{opt subkey2(string)}}
Associative array subkey #2.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{marker examples}{...}
{title:Examples}

{smcl}
INCLUDE help ddml_example_extract.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
