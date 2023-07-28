{smcl}
{* *! version 28july2023}{...}
{viewerjumpto "Syntax" "ddml_export##syntax"}{...}
{viewerjumpto "Examples" "ddml_export##examples"}{...}
{viewerjumpto "Installation" "ddml_export##installation"}{...}
{viewerjumpto "References" "ddml_export##references"}{...}
{viewerjumpto "Authors" "ddml_export##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "ddml interactive" "ddml interactive"}{...}
{vieweralsosee "Other" "ddml_export##also_see"}{...}
{hline}
{cmd:help ddml overlap}{right: v1.4.1}
{hline}

{title:ddml overlap commands for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{cmd:ddml overlap} reports overlap plots following estimation of the {opt ddml} { {it:interactive} and {it:interactiveiv} models.
One overlap (line) plot of propensity scores is reported for each treatment variable learner;
by default, propensity scores for all crossfit samples are plotted.
Overlap plots for the treatment variables are combined using {helpb graph combine}.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml overlap} [ {opt mname(name)} {opt replist(numlist)} {opt pslist(namelist)} {opt n(integer)} {opt kernel(name)}
{opt name(name [, replace])} {opt title(string)} {opt subtitle(string)} {opt lopt0(string)}
{opt lopt1(string)}{bind: ]}

{synoptset 20}{...}
{synopthdr:Options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt replist(numlist)}}
list of crossfitting resamples to plot. Defaults to all.
{p_end}
{synopt:{opt pslist(namelist)}}
varnames of propensity scores to plot (excluding the resample number). Defaults to all.
{p_end}
{synopt:{opt n(integer)}}
see {helpb teffects overlap}.
{p_end}
{synopt:{opt kernel(name)}}
see {helpb teffects overlap}.
{p_end}
{synopt:{opt name(name)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt title(string)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt subtitle(string)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt lopt0(string)}}
options for line plot of untreated; default is solid/navy; see {helpb line}.
{p_end}
{synopt:{opt lopt0(string)}}
options for line plot of treated; default is short dash/dark orange; see {helpb line}.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{title:Examples}

{smcl}
INCLUDE help ddml_example_overlap.sthlp


{smcl}
INCLUDE help ddml_install_ref_auth
