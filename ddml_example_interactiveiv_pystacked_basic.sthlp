{smcl}
{* *! version 25july2023}{...}

{pstd}We use {help pystacked} with two base learners for each reduced form equation.{p_end}

{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init interactiveiv"}{p_end}
{phang2}. {stata "ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(ols lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}
