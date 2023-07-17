{smcl}
{pstd}{ul:Flexible partially-linear IV model - Detailed example with {help pystacked}}

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear"}{p_end}
{phang2}. {stata "global Y share"}{p_end}
{phang2}. {stata "global D price"}{p_end}
{phang2}. {stata "global X hpwt air mpd space"}{p_end}
{phang2}. {stata "global Z sum*"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We initialize the model.{p_end}

{phang2}. {stata "ddml init fiv"}{p_end}

{pstd}We add learners for E[Y|X] in the usual way. Here we use {helpb pystacked}'s default base learners.{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg)"}{p_end}

{pstd}Adding learners for E[D|Z,X] and E[D|X] in the {opt fiv} model is different
from how it's done in the {opt partialiv} model.
The reason for this is that the estimation of E[D|X]
depends on the estimation of E[D|X,Z].
More precisely, we first obtain the fitted values D^=E[D|X,Z] and 
fit these against X to estimate E[D^|X].{p_end}

{pstd}
When adding learners for E[D|Z,X],
we need to provide a name for each learner using {opt learner(name)}
Here we use "Dhat_pys".{p_end}

{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_pys): pystacked $D $X $Z, type(reg)"}{p_end}

{pstd}When adding learners for E[D|X], we explicitly refer to the learner from 
the previous step (here, "Dhat_pys")
and also provide the name of the treatment variable ({cmd:vname($D)}).
Finally, we use the placeholder {cmd:{D}} in place of the dependent variable.{p_end}

{phang2}. {stata "ddml E[D|X], learner(Dhat_pys) vname($D): pystacked {D} $X, type(reg)"}{p_end}
 
{pstd}That's it. Now we can move to cross-fitting and estimation.{p_end}

{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}To replicate what {cmd:ddml} does in the background:{p_end}

{phang2}. {stata "cap drop Ytilde"}{p_end}
{phang2}. {stata "cap drop Dtilde"}{p_end}
{phang2}. {stata "cap drop Ztilde"}{p_end}
{phang2}. {stata "gen double Ytilde = $Y - Y1_pystacked_1"}{p_end}
{phang2}. {stata "gen Dtilde = $D - Dhat_pys_h_1"}{p_end}
{phang2}. {stata "gen Zopt = Dhat_pys_1 - Dhat_pys_h_1"}{p_end}
{phang2}. {stata "ivreg Ytilde (Dtilde=Zopt), robust"}{p_end}
