{smcl}
{* *! version 26july2023}{...}
{smcl}
{pstd}{ul:Partially-linear model - Detailed general example with multiple learners} 

{pstd}Here we used {opt ddml} to add learners. This allows use of learners not supported by,
or as alternatives to, those available via {help pystacked}.
It is also possible to use {help pystacked} as a standalone learner in this way.{p_end}

{pstd}Preparation: load the data and define the globals.
Use the name "m1" for this new estimation, 
to distinguish it from any model estimated previously that uses the default name "m0".
This enables having multiple estimations available for comparison.
Also specify 5 resamplings.{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(5) mname(m1)"}{p_end}

{pstd}We add supervised machine learners for estimating the conditional 
expectation E[Y|X].
In each step, we add the {opt mname(m1)} option to ensure that the learners
are added to correct model.
We also specify the names of the variables containing the estimated conditional
expectations using the {opt learner(varname)} option.
This avoids overwriting any variables created some other model using default naming.{p_end}

{pstd} We first add simple linear regression.{p_end}

{phang2}. {stata "ddml E[Y|X], mname(m1) learner(Y_m1_reg): reg $Y $X"}{p_end}

{pstd}We can add more than one learner per reduced form equation. Here, we 
add a random forest learner. We do this using {help pystacked} to implement a single learner.{p_end}

{phang2}. {stata "ddml E[Y|X], mname(m1) learner(Y_m1_pys): pystacked $Y $X, type(reg) method(rf)"}{p_end}

{pstd}We do the same for the conditional expectation E[D|X].{p_end}

{phang2}. {stata "ddml E[D|X], mname(m1) learner(D_m1_reg): reg $D $X"}{p_end}
{phang2}. {stata "ddml E[D|X], mname(m1) learner(D_m1_pys): pystacked $D $X, type(reg) method(rf)"}{p_end}

{pstd}Check if learners were correctly added:{p_end}

{phang2}. {stata "ddml desc, mname(m1) learners"}{p_end}

{pstd}Cross-fitting and estimation.
Since we added two learners for each of our two reduced form equations, 
there are four possible specifications. 
By default, the result shown corresponds to the specification 
with the lowest out-of-sample MSPE:{p_end}

{phang2}. {stata "ddml crossfit, mname(m1)"}{p_end}
{phang2}. {stata "ddml estimate, mname(m1) robust"}{p_end}

{pstd}To estimate all four specifications, we use the {cmd:allcombos} option:{p_end}

{phang2}. {stata "ddml estimate, mname(m1) robust allcombos"}{p_end}

{pstd}After having estimated all specifications, we can retrieve 
specific results. Here we use the specification relying on OLS for both
estimating both E[Y|X] and E[D|X], from the 4th cross-fit split ({opt rep(4))}.
(Note: Working interactively, the simplest way to do this
is to click on the hyperlink in the summary table in the {opt ddml estimate} output above.)
The {opt notable} option suppresses the summary table:{p_end}

{phang2}. {stata "ddml estimate, mname(m1) spec(1) rep(4) replay notable"}{p_end}

{pstd}You could manually retrieve the same point estimate by 
cacluating the orthogonalized versions of {opt net_tfa} and {opt e401}
from the 4th cross-fit estimation and then using {help regress}:{p_end}

{phang2}. {stata "cap drop Yresid"}{p_end}
{phang2}. {stata "cap drop Dresid"}{p_end}
{phang2}. {stata "gen double Yresid = $Y - Y_m1_reg_4"}{p_end}
{phang2}. {stata "gen double Dresid = $D - D_m1_reg_4"}{p_end}
{phang2}. {stata "regress Yresid Dresid"}{p_end}

{pstd}You can also compare the estimated conditional expectations graphically:{p_end}

{phang2}. {stata "twoway (scatter $Y Y_m1_pys_4) "}{p_end}

{pstd}To describe the ddml model setup or results in detail,
you can use {cmd: ddml describe} with the relevant option ({opt sample}, {opt learners}, {opt crossfit}, {opt estimates}),
or just describe them all with the {opt all} option:{p_end}

{phang2}. {stata "ddml describe, mname(m1) all"}{p_end}

{pstd}If there is a previously-estimated {opt ddml} model called "m0",
we can load it using {opt ddml estimate} with the {opt mname(m0)} and {opt replay} options and compare.{p_end}

{phang2}. {stata "ddml estimate, mname(m0) replay"}{p_end}
