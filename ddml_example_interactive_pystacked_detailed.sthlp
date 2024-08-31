{smcl}
{* *! version 30aug2024}{...}
{smcl}
{pstd}{ul:Interactive model - Detailed example with {help pystacked}}{p_end}

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}

{phang2}. {stata "webuse cattaneo2, clear"}{p_end}
{phang2}. {stata "global Y bweight"}{p_end}
{phang2}. {stata "global D mbsmoke"}{p_end}
{phang2}. {stata "global X prenatal1 mmarried fbaby mage medu"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}We use 5 folds and 5 resamplings; that is, 
we estimate the model 5 times using randomly chosen folds.{p_end}

{phang2}. {stata "ddml init interactive, kfolds(5) reps(5)"}{p_end}

{pstd}We need to estimate the conditional expectations of E[Y|X,D=0], 
E[Y|X,D=1] and E[D|X]. The first two conditional expectations 
are added jointly.{p_end} 
{pstd}We consider two supervised learners: linear regression and gradient boosted
trees, stacked using {helpb pystacked}.
Note that we use gradient boosted regression trees for E[Y|X,D], but
gradient boosted classification trees for E[D|X].{p_end}

{phang2}. {stata "ddml E[Y|X,D]: pystacked $Y $X, type(reg) methods(ols gradboost)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X, type(class) methods(logit gradboost)"}{p_end}

{pstd}Cross-fitting and short-stacking:{p_end}

{phang2}. {stata "ddml crossfit, shortstack"}{p_end}

{pstd}In the final estimation step, we can estimate
the average treatment effect (the default),
the average treatment effect on the treated ({opt atet}),
or the average treatment effect on the untreated ({opt ateu}).{p_end}

{phang2}. {stata "ddml estimate"}{p_end}
{phang2}. {stata "ddml estimate, atet"}{p_end}
{phang2}. {stata "ddml estimate, ateu"}{p_end}

{pstd}Recall that we have specified 5 resampling iterations ({opt reps(5)})
By default, the median over short-stacked resampling iterations is shown.
At the bottom, a table of summary statistics over resampling iterations is shown. 
To display the mean over standard stacking results, i.e.,
the results where the weights derive from {helpb pystacked} and vary by cross-fit fold,
we use {opt ddml estimate, replay} with {opt spec(st)} and {opt rep(mn)}.{p_end}

{phang2}. {stata "ddml estimate, spec(st) rep(mn) notable replay"}{p_end}

{pstd}Generate an overlap plot using {opt ddml overlap}:{p_end}

{phang2}. {stata "ddml overlap"}{p_end}

{pstd}Report the standard stacking and short-stacking weights:{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}
