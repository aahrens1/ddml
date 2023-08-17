{smcl}
{* *! version 17aug2023}{...}
{marker models}{...}
{title:Models}

{pstd}
Throughout we use {it:Y} to denote the outcome variable, 
{it:X} to denote confounders, 
{it:Z} to denote instrumental variable(s), and
{it:D} to denote the treatment variable(s) of interest.

{pstd}
{ul:Partially-linear model} [{it:partial}]

	Y = {it:a}.D + g(X) + U
        D = m(X) + V

{pstd}
where the aim is to estimate {it:a} while controlling for X. To this end, 
we estimate the conditional expectations
E[Y|X] and E[D|X] using a supervised machine learner.

{pstd}
{ul:Interactive model} [{it:interactive}]

	Y = g(X,D) + U
        D = m(X) + V

{pstd}
which relaxes the assumption that X and D are separable. 
D is a binary treatment variable,
and we aim to estimate the average treatment effect (ATE)
or average treatment effect on the treated (ATET).
We estimate the conditional expectations E[D|X], as well as 
E[Y|X,D=0] and E[Y|X,D=1] (jointly added using {cmd:ddml E[Y|X,D]}).

{pstd}
{ul:Partially-linear IV model} [{it:iv}]

	Y = {it:a}.D + g(X) + U
        Z = m(X) + V

{pstd}
where the aim is to estimate {it:a}. 
We estimate the conditional expectations E[Y|X], 
E[D|X] and E[Z|X] using a supervised machine
learner.

{pstd}
{ul:Flexible partially-linear IV model} [{it:fiv}]

	Y = {it:a}.D + g(X) + U
        D = m(Z) + g(X) + V 

{pstd}
where the estimand of interest is {it:a}. 
We estimate the conditional expectations
E[Y|X], 
E[D^|X] and D^:=E[D|Z,X] using a supervised machine
learner. The instrument is then formed as D^-E^[D^|X] where E^[D^|X] denotes
the estimate of E[D^|X]. 

{pstd}
Note: "{D}" is a placeholder that is used because last step (estimation of E[D|X]) 
uses the fitted values from estimating E[D|X,Z].

{pstd}
{ul:Interactive IV model}  [{it:interactiveiv}]

	Y = g(Z,X) + U
        D = h(Z,X) + V
        Z = m(X) + E

{pstd}
where the aim is to estimate the local average treatment effect (LATE).
We estimate, using a supervised machine
learner, the following conditional expectations:
E[Y|X,Z=0] and E[Y|X,Z=1] (jointly added using {cmd:ddml E[Y|X,Z]});
E[D|X,Z=0] and E[D|X,Z=1] (jointly added using {cmd:ddml E[D|X,Z]});
and E[Z|X].


{marker estimation}{...}
{title:Main steps when estimating with ddml}

{pstd}Estimation with {cmd:ddml} proceeds in four steps. 

{pstd}
{ul:Step 1.} Initialize {cmd:ddml} and select model:

{p 8 14}{cmd:ddml init}
{it:model} [if] [in]
[ , {opt mname(name)} {opt kfolds(integer)}
{opt fcluster(varname)}
{opt foldvar(varlist)} {opt reps(integer)} 
{opt norandom} {opt tabfold} {opt vars(varlist)}{bind: ]}

{pstd}
where {it:model} is either {it:partial}, {it:iv}, {it:interactive}, {it:fiv}, {it:interactiveiv};
see {help ddml##models:model descriptions}.

{pstd}
{ul:Step 2.} Add supervised ML programs for estimating conditional expectations:

{p 8 14}{cmd:ddml} {it:eq} 
[ , {opt mname(name)} {opt vname(varname)} {opt l:earner(varname)}
{opt vtype(string)}
{opt predopt(string)}{bind: ] :}
{it:command} {it:depvar} {it:vars} [ , {it:cmdopt}{bind: ]}

{pstd}
where, depending on model chosen in Step 1,
{it:eq} is either 
{it:E[Y|X]} {it:E[Y|D,X]} {it:E[Y|X,Z]} {it:E[D|X]} {it:E[D|X,Z]} {it:E[Z|X]}.
{it:command} is a supported supervised ML program (e.g. {help pystacked} or {help cvlasso}). 
See {help ddml##compatibility:supported programs}.

{pstd}
Note: Options before ":" and after the first comma refer to {cmd:ddml}. 
Options that come after the final comma refer to the estimation command. 
{p_end}

{pstd}
{ul:Step 3.} Cross-fitting:

{p 8 14}{cmd:ddml crossfit} [ , {opt mname(name)} {opt shortstack}{bind: ]} 

{pstd}
This step implements the cross-fitting algorithm. Each learner is fitted iteratively on training folds and out-of-sample predicted values are obtained.

{pstd}
{ul:Step 4.} Estimate causal effects:

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {cmdab:r:obust} {opt cluster(varname)} {opt vce(type)} {opt atet} {opt ateu} {opt trim(real)}{bind: ]} 

{pstd}
The {cmd:ddml estimate} command returns treatment effect estimates for all combination of learners 
added in Step 2.

{pstd}
{ul:Optional.} Report/post selected results:

{p 8 14}{cmd:ddml estimate} [ , {opt mname(name)} {opt spec(integer or string)} {opt rep(integer or string)} {opt allcombos} {opt not:able} {opt replay} {bind: ]} 

{pstd}
{marker auxiliary}{...}
{ul:Optional.} Retrieve information from {cmd:ddml}:

{p 8 14}{cmd:ddml extract} [ {it:object_name} , {opt mname(name)} {opt show(display_item)} {opt ename(name)} {opt vname(varname)}
{opt stata} {opt keys} {opt key1(string)} {opt key2(string)} {opt key3(string)} {opt subkey1(string)}
{opt subkey2(string)}{bind: ]}

{pstd}
{it:display_item} can be {it:stweights}, {it:ssweights}, {it:psweights}, {it:weights}, {it:mse}, {it:n}, or {it:pystacked}.
{cmd:ddml} stores many internal results on associative arrays.
See {help ddml extract} for details.

{pstd}
For full details and further options, follow the links to the detailed help files {help ddml##help:above}.


{marker compatibility}{...}
{title:Compatible programs}

{pstd}
{marker general}{...}
{opt ddml} is compatible with a large set of user-written Stata commands. 
It has been tested with 

{p 7 9 0} 
- the {help pystacked} package (see {help pystacked} and {help ddml##pystacked:below}). 
Note that {help pystacked} requires Stata 16.

{p 7 9 0} 
- {help lassopack} for regularized regression (see {help lasso2}, {help cvlasso}, {help rlasso}).

{p 7 9 0} 
- {help rforest} by Zou & Schonlau. Note that {cmd:rforest} requires the option 
{cmd:vtype(none)}. 

{p 7 9 0} 
- {help svmachines} by Guenther & Schonlau.

{pstd}
Beyond these, it is compatible with any Stata program that 

{p 7 9 0} 
- uses the standard "{it:reg y x}" syntax,

{p 7 9 0} 
- supports {it:if}-conditions,

{p 7 9 0} 
- and comes with {help predict} post-estimation programs.

{marker pystacked}{...}
{pstd}
{help pystacked} implements stacking regression ({help pystacked##Wolpert1992:Wolpert, 1992})
via Stata's Python integration in combination with
{browse "https://scikit-learn.org/stable/index.html":scikit-learn}'s 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingRegressor.html":sklearn.ensemble.StackingRegressor} and 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingClassifier.html":sklearn.ensemble.StackingClassifier}. 
Stacking is a way of combining multiple supervised
machine learners (the "base" or "level-0" learners) into
a meta learner.

{pstd}{help pystacked} is the recommended way to specify multiple learners in {opt ddml}.
{help pystacked} provides a fast way of estimating all the learners in a single call to one program,
and {opt ddml} has integrated support for various features provided by {help pystacked}.
{opt ddml} will store the predicted values of the specified base learners as well as the combined ("stacked") predicted values.
It also stores the stacking weights used by {help pystacked} along with the {opt ddml} short-stacking weights.
See {help ddml stacking} for more details.
