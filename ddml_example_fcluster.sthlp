{smcl}
{* *! version 17aug2023}{...}
{smcl}
{pstd}{ul:Cluster sampling with cross-fit folds - Basic example with {help pystacked}}{p_end}

{pstd}Load the data, define global macros and set the seed.{p_end}

{phang2}. {stata "webuse nlsw88, clear"}{p_end}
{phang2}. {stata "gen lwage = ln(wage)"}{p_end}
{phang2}. {stata "global Y lwage"}{p_end}
{phang2}. {stata "global D union"}{p_end}
{phang2}. {stata "global X age-c_city hours-tenure"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Initialize the model.
The {opt fcluster(industry)} ("fold-cluster") option tells {opt ddml}
to ensure that clusters (here, identified by the variable {opt industry})
are not split across cross-fit folds, i.e., each cluster appears in only one cross-fit fold.
Here we specify 2 cross-fit folds,
so all observations for each cluster will appear in either fold 1 or in fold 2.
NB: This example is somewhat artificial, because there are only 12 clusters (industries).{p_end}

{phang2}. {stata "ddml init partial, kfolds(2) fcluster(industry)"}{p_end}
{phang2}. {stata "tab industry m0_fid_1"}{p_end}

{pstd}Since there are 12 clusters defined by {opt industry},
we could achieve the same cross-fit split either by specifying {opt fcluster(industry)},
or by using {opt fcluster(industry)} as the fold identifier and specifying {opt foldvar(industry)}.
(NB: The split is the same but the fold numbering is different.){p_end}

{phang2}. {stata "ddml init partial, foldvar(industry)"}{p_end}
{phang2}. {stata "tab industry m0_fid_1"}{p_end}

{phang2}. {stata "ddml init partial, kfolds(12) fcluster(industry)"}{p_end}
{phang2}. {stata "tab industry m0_fid_1"}{p_end}

{pstd}Estimation is standard,
but to obtain cluster-robust SEs the covariance estimator
needs to be requested with {opt ddml estimate}:{p_end}

{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate, cluster(industry)"}{p_end}
