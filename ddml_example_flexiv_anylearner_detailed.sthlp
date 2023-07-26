{smcl}
{pstd}{ul:Flexible partially-linear IV model - Detailed example with {help pystacked}}

{pstd}Here will illustrate how to do standard- and short-stacking with the flexible IV model.{p_end}

{pstd}Note: Support for {help pystacked} integration is relatively limited for the flexible IV model.
In particular, short-stacking requires that individual learners appear in separate {help pystacked} commands,
pooled stacking is not available, and re-stacking with {opt ddml estimate} is also not available.{p_end}

{pstd}First we illustrate how to do standard stacking with {help pystacked}.
To start, we load the data, define global macros, set the seed and initialize the model.{p_end}

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
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}The stacking weights and other results from {help pystacked} are available via {help ddml extract}.
Note this is only the case if {help pystacked} is the single learner in each equation.{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}

{pstd}To replicate what {cmd:ddml} does in the background:{p_end}

{phang2}. {stata "cap drop Ytilde"}{p_end}
{phang2}. {stata "cap drop Dtilde"}{p_end}
{phang2}. {stata "cap drop Ztilde"}{p_end}
{phang2}. {stata "gen double Ytilde = $Y - Y1_pystacked_1"}{p_end}
{phang2}. {stata "gen Dtilde = $D - Dhat_pys_h_1"}{p_end}
{phang2}. {stata "gen Zopt = Dhat_pys_1 - Dhat_pys_h_1"}{p_end}
{phang2}. {stata "ivreg Ytilde (Dtilde=Zopt), robust"}{p_end}

{pstd}Next we illustrate how to do short-stacking with the flexible IV model.
We again use {help pystacked}, but the procedure applies to any set of learners.
Here we use the same learners as in the standard stacking estimation with {help pystacked} above,
in order to facilitate direct comparison of the two sets of results.
We begin by re-initializing the model.{p_end}

{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init fiv"}{p_end}

{pstd}We add learners for E[Y|X] in the usual way,
but we need to specify each {help pystacked} learner in a separate equation.{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg) m(ols)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg) m(lassocv)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X, type(reg) m(gradboost)"}{p_end}

{pstd}
As above, when adding learners for E[D|Z,X],
we need to provide a name for each learner using {opt learner(name)}.{p_end}

{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_ols): pystacked $D $X $Z, type(reg) m(ols)"}{p_end}
{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_lassocv): pystacked $D $X $Z, type(reg) m(lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|Z,X], learner(Dhat_gradboost): pystacked $D $X $Z, type(reg) m(gradboost)"}{p_end}

{pstd}Again as above, when adding learners for E[D|X],
we explicitly refer to the learner from the previous step,
the name of the treatment variable ({cmd:vname($D)}),
and the placeholder {cmd:{D}} in place of the dependent variable.{p_end}

{phang2}. {stata "ddml E[D|X], learner(Dhat_ols) vname($D): pystacked {D} $X, type(reg) m(ols)"}{p_end}
{phang2}. {stata "ddml E[D|X], learner(Dhat_lassocv) vname($D): pystacked {D} $X, type(reg) m(lassocv)"}{p_end}
{phang2}. {stata "ddml E[D|X], learner(Dhat_gradboost) vname($D): pystacked {D} $X, type(reg) m(gradboost)"}{p_end}
 
{pstd}Short-stacking is requested when cross-fitting.
Short-stacking weights can be examined using {help ddml extract}.{p_end}

{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}
