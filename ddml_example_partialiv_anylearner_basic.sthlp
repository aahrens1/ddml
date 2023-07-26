{smcl}
{pstd}{ul:Partially-linear IV model - Basic example with various learners} 

{pstd}Preparation: we load the data, define global macros and set the seed.{p_end}

{phang2}. {stata "use https://statalasso.github.io/dta/AJR.dta, clear"}{p_end}
{phang2}. {stata "global Y logpgp95"}{p_end}
{phang2}. {stata "global D avexpr"}{p_end}
{phang2}. {stata "global Z logem4"}{p_end}
{phang2}. {stata "global X lat_abst edes1975 avelf temp* humid* steplow-oilres"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Preparation: we load the data, define global macros and set the seed. Since the
data set is very small, we consider 30 cross-fitting folds.{p_end}

{phang2}. {stata "ddml init iv, kfolds(30)"}{p_end}

{pstd}The partially linear IV model has three conditional expectations: 
E[Y|X], E[D|X] and E[Z|X]. For each reduced form equation, we use two learners:
OLS and random forest.
To illustrate how {opt ddml} works with other packages,
instead of a single call to {opt pystacked} specifying two base learners
we specify Stata's {help regress} and {help rforest} by Zou and Schonlau as the two learners.
We need to add the option {opt vtype(none)} for {help rforest} to 
work with {cmd:ddml} since {help rforest}'s {cmd:predict} command doesn't
support variable types.{p_end}

{phang2}. {stata "ddml E[Y|X]: reg $Y $X"}{p_end}
{phang2}. {stata "ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: reg $D $X"}{p_end}
{phang2}. {stata "ddml E[D|X], vtype(none): rforest $D $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[Z|X]: reg $Z $X"}{p_end}
{phang2}. {stata "ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)"}{p_end}

{pstd}Cross-fitting and estimation; report all combinations of estimated conditional expectations.{p_end}

{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust allcombos"}{p_end}
	
