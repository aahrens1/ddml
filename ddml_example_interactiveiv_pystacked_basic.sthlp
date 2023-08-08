{smcl}
{* *! version 3august2023}{...}
{smcl}
{pstd}{ul:Interactive IV model (LATE) - Basic example with {help pystacked}}

{pstd}We use {help pystacked} with two base learners for each reduced form equation.{p_end}

{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global D training"}{p_end}
{phang2}. {stata "global Z assignmt"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}

{pstd}Drop observations where treatment=1 even though assignment=0.
(Up to the user how to handle such observations;
{opt ddml} can handle these cases,
and we drop them here only to illustrate how this is reflected in the stacking weights.){p_end}

{phang2}. {stata "tab $D $Z"}{p_end}
{phang2}. {stata "drop if $D==1 & $Z==0"}{p_end}

{pstd}Set the seed, initialize, cross-fit and estimate:{p_end}

{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init interactiveiv"}{p_end}
{phang2}. {stata "ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(ols lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(logit lassocv)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Report the stacking weights.
Note that the weights for the non-existent case (treatment=1, assignment=0) are missing.{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
