{smcl}
{* *! version 26july2023}{...}
{smcl}
{pstd}{ul:Flexible partially-linear IV model - Basic example with {help pystacked}}

{pstd}First load the data, define global macros, set the seed and initialize the model.
We add learners for E[Y|X] in the usual way.
We illustrate with single {help pystacked} estimations,
but the procedure applies to all learners.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear"}{p_end}
{phang2}. {stata "global Y share"}{p_end}
{phang2}. {stata "global D price"}{p_end}
{phang2}. {stata "global X hpwt air mpd space"}{p_end}
{phang2}. {stata "global Z sum*"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init fiv"}{p_end}

{pstd}Adding learners for E[Y|X] is the same as for other {opt ddml} linear models:{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg)"}{p_end}

{pstd}Adding learners for E[D|Z,X] and E[D|X] in the {opt fiv} model is different
from how it's done in the {opt partialiv} model.
The reason for this is that the estimation of E[D|X]
depends on the estimation of E[D|X,Z].{p_end}

{pstd}When adding learners for E[D|Z,X],
we need to provide a name for each learners using {opt learner(name)}.
Here we use the name "Dhat_pys".{p_end}

{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_pys): pystacked $D $X $Z, type(reg)"}{p_end}

{pstd}When adding learners for E[D|X], we explicitly refer to the name of the learner from 
the previous step (here, "Dhat_pys").
We also provide the name of the treatment variable ({cmd:vname($D)}),
and we use the placeholder {cmd:{D}} in place of the dependent variable.{p_end}

{phang2}. {stata "ddml E[D|X], learner(Dhat_pys) vname($D): pystacked {D} $X, type(reg)"}{p_end}

{pstd}The crossfit and estimation commands with the {opt fiv} model are standard.{p_end}

{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}
