{smcl}
{* *! version 8aug2023}{...}
{smcl}
{pstd}{ul:Partially-linear model - Multiple treatments with {help pystacked}}

{pstd}We can also run the partially-linear model with multiple treatments. 
In this simple example, we estimate the effect of both 401k elligibility 
{cmd:e401} and education {cmd:educ}. 
Note that we remove {cmd:educ} from the set of controls.
We again use {help pystacked} as the single learner provided to {opt ddml};
the two base learners, OLS and random forest, are provided via {help pystacked}.
We use the simplified syntax supported by {help pystacked}.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D1 e401"}{p_end}
{phang2}. {stata "global D2 educ"}{p_end}
{phang2}. {stata "global X tw age inc fsize db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Initialize the model.{p_end}

{phang2}. {stata "ddml init partial, kfolds(2)"}{p_end}

{pstd}Add learners. Note that we add learners with both {cmd:$D1} and
{cmd:$D2} as the dependent variable.{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg) methods(ols rf)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D1 $X, type(reg) methods(ols rf)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D2 $X, type(reg) methods(ols rf)"}{p_end}

{pstd}Cross-fitting.{p_end}

{phang2}. {stata "ddml crossfit"}{p_end}

{pstd}Estimation.{p_end}

{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Because we have used {help pystacked} as the single {opt ddml} learner,
we can access the saved {opt pystacked} information.
Here we use the {opt pystacked} option to get the stacking weights and MSEs by cross-fit fold:{p_end}

{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}
