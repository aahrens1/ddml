{smcl}
{* *! version 8aug2023}{...}
{viewerjumpto "Syntax" "ddml_init##syntax"}{...}
{viewerjumpto "Options" "ddml_init##options"}{...}
{viewerjumpto "Installation" "ddml_init##installation"}{...}
{viewerjumpto "References" "ddml_init##references"}{...}
{viewerjumpto "Authors" "ddml_init##authors"}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_init##also_see"}{...}
{hline}
{cmd:help ddml init, ddml eq, ddml sample}{right: v1.4.2}
{hline}

{title:ddml init, eq and sample commands for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{opt ddml init} {it:model} initializes the model,
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, or {it:interactiveiv}.

{pstd}
{cmd: ddml eq: command} adds supervised ML programs for estimating conditional expectations,
where {it:eq} is the conditional expectation to be estimated (e.g., {it:E[Y|X]})
and {it:command} is a supported supervised ML program.

{pstd}
{opt ddml sample} adds cross-fitting repetitions to an existing and possibly already-estimated model.


{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml init}
{it:model} [if] [in]
[ , {opt mname(name)}
{opt prefix}
{opt kfolds(integer)}
{opt fcluster(varname)}
{opt foldvar(varlist)}
{opt reps(integer)} 
{opt norandom}
{opt tabfold}
{opt vars(varlist)}{bind: ]}

{pstd}
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, {it:interactiveiv}.

{p 8 14}{cmd:ddml} {it:eq} 
[ , {opt mname(name)}
{opt vname(varname)}
{opt l:earner(varname)}
{opt vtype(string)}
{opt predopt(string)}{bind: ] :}
{it:command} {it:depvar} {it:vars} [ , {it:cmdopt}{bind: ]}

{pstd}
where, depending on model chosen in Step 1,
{it:eq} is either 
{it:E[Y|X]} {it:E[Y|D,X]} {it:E[Y|X,Z]} {it:E[D|X]} {it:E[D|X,Z]} {it:E[Z|X]}.
{it:command} is a supported supervised ML program (e.g. {helpb pystacked} or {helpb cvlasso}).

{pstd}
Note: Options before ":" and after the first comma refer to {cmd:ddml}. 
Options that come after ":" and the final comma refer to the estimation command. 
{p_end}

{p 8 14}{cmd:ddml sample} [ , {opt append}[{cmd:(}{it:integer}{cmd:)}] {opt foldvar(varlist)} {bind: ]}

{pstd}
adds cross-fitting repetitions to an existing and possibly already-estimated model,
where the additional repetitions is indicated either by {opt append(#)}
or by {opt append} and the cross-fit fold identifiers in {opt foldvar(varlist)}.


{marker options}{...}
{synoptset 20}{...}
{synopthdr:init options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Allows to run multiple DDML
models simultaneously. Defaults to {it:m0}.
{p_end}
{synopt:{opt prefix}}
tells {opt ddml} to prefix the names of all created variables
with name of the DDML model.
Default is to prefix only the created sample and fold ID variables.
{p_end}
{synopt:{opt kfolds(integer)}}
number of cross-fitting folds. The default is 5.
{p_end}
{synopt:{opt fcluster(varname)}}
cluster identifiers for cluster randomization of random folds.
{p_end}
{synopt:{opt foldvar(varlist)}}
integer variable with user-specified cross-fitting folds (one per cross-fitting repetition).
{p_end}
{synopt:{opt norandom}}
use observations in existing order instead of randomizing before splitting into folds;
if multiple resamples, applies to first resample only;
ignored if user-defined fold variables are provided in {opt foldvar(varlist)}.
{p_end}
{synopt:{opt reps(integer)}}
cross-fitting repetitions, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds. 
{p_end}
{synopt:{opt tabfold}}
prints a table with frequency of observations by fold.
{p_end}
{synopt:{opt vars(varlist)}}
tells {opt ddml} that the variables in {it:varlist} are used in the estimation.
Useful if you want the fold split to take account of
observations dropped because of missing values.
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:equation options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt vname(varname)}}
name of the dependent variable in the reduced form estimation. 
This is usually inferred from the command line but is mandatory
for the {it:fiv} model.
{p_end}
{synopt:{opt l:earner(varname)}}
optional name of the variable to be created. 
{p_end}
{synopt:{opt vtype(string)}}
(rarely used) optional variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank 
(required when using {cmd:ddml} with {helpb rforest}.)
{p_end}
{synopt:{opt predopt(string)}}
(rarely used) {cmd:predict} option to be used to get predicted values. 
Typical values could be {opt xb} or {opt pr}. Default is 
blank. 
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}

{synoptset 20}{...}
{synopthdr:sample options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt append(#)}}
number of additional resamples to cross-fit.
{p_end}
{synopt:{opt append}}
when no number of resamples to append is provided, this is based on the list fold IDs in {opt foldvar(varlist)}.
{p_end}
{synopt:{opt foldvar(varlist)}}
integer variable with user-specified cross-fitting folds (one per cross-fitting repetition).
{p_end}
{synoptline}
{p2colreset}{...}
{pstd}


{smcl}
INCLUDE help ddml_install_ref_auth
