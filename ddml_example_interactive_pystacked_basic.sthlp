{smcl}
{* *! version 28july2023}{...}
{smcl}
{pstd}{ul:Interactive model - Basic example with {help pystacked}}{p_end}

{pstd}We need to estimate the conditional expectations of E[Y|X,D=0], E[Y|X,D=1] and E[D|X].
The first two conditional expectations are added jointly.
We use 5 cross-fit folds and 2 resamplings
(more resamplings would be advisable; we use 2 in this example so the code runs faster).
We specify two supervised learners: linear regression and gradient boosted
trees, stacked using {help pystacked}.
We use {help pystacked}'s 2nd syntax and stack using the single-best learner
(rather than the default constrained least squares).
Note that we use gradient boosted regression trees for E[Y|X,D],
but gradient boosted classification trees for E[D|X].{p_end}

{phang2}. {stata "webuse cattaneo2, clear"}{p_end}
{phang2}. {stata "global Y bweight"}{p_end}
{phang2}. {stata "global D mbsmoke"}{p_end}
{phang2}. {stata "global X prenatal1 mmarried fbaby mage medu"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init interactive, kfolds(5) reps(2)"}{p_end}
{phang2}. {stata "ddml E[Y|X,D]: pystacked $Y $X || method(ols) || method(gradboost) || , type(reg) finalest(singlebest)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X || method(logit) || method(gradboost) || , type(class) finalest(singlebest)"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}

{pstd}{opt ddml estimate} reports the ATE (average treatment effect) by default:{p_end}

{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Request the ATET (average treatment effect on the treated) instead:{p_end}

{phang2}. {stata "ddml estimate, atet"}{p_end}
