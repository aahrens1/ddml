{smcl}
{* *! version 27july2023}{...}
{smcl}
{pstd}{ul:Interactive IV model (LATE) - Detailed example with {help pystacked}}

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}

{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We initialize the model.{p_end}

{phang2}. {stata "ddml init interactiveiv, kfolds(5)"}{p_end}

{pstd}We use {helpb pystacked} with two base learners for each reduced form equation.
Note that E[Y|X,Z] is a regression problem,
whereas E[D|X,Z] and E[Z|X] are classification problems.{p_end}

{phang2}. {stata "ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(ols lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}

{pstd}Cross-fitting and estimation, with short-stacking implemented via {opt ddml}.{p_end}

{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}Compare the short-stacking estimation above with standard (within-cross-fit-fold) stacking:{p_end}

{phang2}. {stata "ddml estimate, spec(st) rep(1) replay notable"}{p_end}

{pstd}Short-stacking is typically considerably faster than standard stacking.
We can estimate using short-stacking only by specifying the {opt nostd} option when cross-fitting.
We re-set the seed for comparability.{p_end}

{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml crossfit, shortstack nostdstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}
