{smcl}
{* *! version 17aug2023}{...}
{smcl}
{pstd}{ul:Partially-linear IV model - Basic example with {help pystacked}} 

{pstd}The model has three conditional expectations: E[Y|X], E[D|X] and E[Z|X].
For each reduced form equation, we use {help pystacked}'s default learners: 
OLS, cross-validated lasso, and gradient boosting.
Since the data set is very small, we consider 30 cross-fitting folds.
NB: The model specification and results will be stored on a Mata object
with the default name "m0".{p_end}

{phang2}. {stata "use https://statalasso.github.io/dta/AJR.dta, clear"}{p_end}
{phang2}. {stata "global Y logpgp95"}{p_end}
{phang2}. {stata "global D avexpr"}{p_end}
{phang2}. {stata "global Z logem4"}{p_end}
{phang2}. {stata "global X lat_abst edes1975 avelf temp* humid* steplow-oilres"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init iv, kfolds(30)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}
{phang2}. {stata "ddml E[Z|X]: pystacked $Z $X"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Replicate the {opt ddml estimate} results for the 1st cross-fit estimation (resample 1) by hand,
using the estimated conditional expectations generated by {opt ddml} and {help pystacked};
"_1" means resample 1.
Compare using {opt ddml estimate, replay}.{p_end}

{phang2}. {stata "cap drop Yresid"}{p_end}
{phang2}. {stata "cap drop Dresid"}{p_end}
{phang2}. {stata "cap drop Zresid"}{p_end}
{phang2}. {stata "gen double Yresid = $Y - Y1_pystacked_1"}{p_end}
{phang2}. {stata "gen double Dresid = $D - D1_pystacked_1"}{p_end}
{phang2}. {stata "gen double Zresid = $Z - Z1_pystacked_1"}{p_end}
{phang2}. {stata "ivreg Yresid (Dresid=Zresid)"}{p_end}
{phang2}. {stata "ddml estimate, mname(m0) spec(st) rep(1) notable replay"}{p_end}

