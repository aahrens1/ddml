{smcl}
{* *! version 28july2023}{...}
{smcl}
{pstd}{ul:ddml describe utility - Basic example with {help pystacked}}{p_end}

{pstd}Load the data, define global macros, set the seed and initialize the model.
Use 2-fold cross-fitting with two repetitions (resamples)
Use {help pystacked}'s default learners as the supervised learners.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked $D $X"}{p_end}
{phang2}. {stata "ddml crossfit"}{p_end}
{phang2}. {stata "ddml estimate"}{p_end}

{pstd}Default of {opt ddml describe} is to report a brief summary.{p_end}

{phang2}. {stata "ddml describe"}{p_end}

{pstd}Options: report details of the total and cross-fit samples,
learners (including the estimation strings),
cross-fit results,
and estimation results.{p_end}

{phang2}. {stata "ddml describe, sample"}{p_end}
{phang2}. {stata "ddml describe, learners"}{p_end}
{phang2}. {stata "ddml describe, crossfit"}{p_end}
{phang2}. {stata "ddml describe, estimates"}{p_end}

{pstd}The {opt all} option is equivalent to specifying all 4 options.{p_end}

{phang2}. {stata "ddml describe, all"}{p_end}
{phang2}. {stata "ddml describe, sample learners crossfit estimates"}{p_end}
