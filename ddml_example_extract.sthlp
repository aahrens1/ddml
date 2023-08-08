{smcl}
{* *! version 31july2023}{...}
{smcl}
{pstd}{ul:ddml extract utility: Extracting stored information from ddml associative arrays}

{pstd}The examples below use the partially-linear model
and stacking regression using {helpb pystacked}.
We also request short-stacking.
The model name is the default name "m0".
For simplicity we use {helpb pystacked}'s default learners and settings.
{p_end}

{pstd}Preparation and estimation:{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "ddml E[Y|X]: pystacked net_tfa $X, type(reg)"}{p_end}
{phang2}. {stata "ddml E[D|X]: pystacked e401 $X, type(reg)"}{p_end}
{phang2}. {stata "ddml crossfit, shortstack"}{p_end}
{phang2}. {stata "ddml estimate, robust"}{p_end}

{pstd}{ul:{opt show(something)} option examples}{p_end}

{pstd}{opt show} option examples: report standard (pystacked) and short-stacked weights.
Standard stacking weights displayed here are mean weights across cross-fit folds.{p_end}

{phang2}. {stata "ddml extract, show(stweights)"}{p_end}
{phang2}. {stata "ddml extract, show(ssweights)"}{p_end}

{pstd}The {opt show} option leaves results in r(.) macros.{p_end}

{phang2}. {stata "mat list r(Y_net_tfa_ss)"}{p_end}
{phang2}. {stata "mat list r(D_e401_ss)"}{p_end}

{pstd}{opt show} option examples: examine the learner weights and MSEs by fold reported by {cmd:pystacked}.{p_end}

{phang2}. {stata "ddml extract, show(pystacked)"}{p_end}

{pstd}{ul:List keys examples}{p_end}

{pstd}List keys of associative arrays used in model m0.
Associative array m0.eqnAA is an "equation AA" and has one key,
which is is the name of the variable for which conditional expectations are estimated.
Associative array m0.estAA is an "estimation AA" and has two keys.
The objects stored on this AA are either estimation results,
AAs that have sets of estimation results, or objects with information about the estimations.{p_end}

{phang2}. {stata "ddml extract, keys"}{p_end}

{pstd}List keys relating to equation for D variable, e401.
Keys for two associative arrays are reported.
Associative array e401.lrnAA is a "learner AA" and has two keys; it stores e.g. an estimation specification.
Associative array e401.resAA is a "results AA" and has three keys; it stores e.g. estimation results.{p_end}

{phang2}. {stata "ddml extract, keys vname(e401)"}{p_end}

{pstd}{ul:Working with model estimation results}{p_end}

{pstd}Extract the estimated beta for the short-stack specification ("ss"), resample 2.
Provide the keys for the AA with the results for the specification and resampling,
and the subkeys for this AA to obtain the posted beta.{p_end}

{phang2}. {stata "ddml extract, key1(ss) key2(2) subkey1(b) subkey2(post)"}{p_end}

{pstd}As above, but store as a Mata object "bmat".
This is done by providing this name after "ddml extract".{p_end}

{phang2}. {stata "ddml extract bmat, key1(ss) key2(2) subkey1(b) subkey2(post)"}{p_end}
{phang2}. {stata "mata: bmat"}{p_end}

{pstd}By default, the object is saved as a Mata object.
To save as a Stata macro r(bmat), use the {opt Stata} option:{p_end}

{phang2}. {stata "ddml extract bmat, key1(ss) key2(2) subkey1(b) subkey2(post) stata"}{p_end}
{phang2}. {stata "mat list r(bmat)"}{p_end}

{pstd}More examples of the above, relating to specification ss and
various resamples or the mean/median across resamples.
(The list of available results was already displayed above by {stata "ddml extract, keys"}.){p_end}

{phang2}. {stata "ddml extract, key1(ss) key2(1) subkey1(D_e401_ss_mse) subkey2(scalar)"}{p_end}
{phang2}. {stata "ddml extract, key1(ss) key2(2) subkey1(D_e401_ss_mse_folds) subkey2(matrix)"}{p_end}
{phang2}. {stata "ddml extract, key1(ss) key2(mn) subkey1(V) subkey2(post)"}{p_end}
{phang2}. {stata "ddml extract, key1(ss) key2(md) subkey1(title) subkey2(local)"}{p_end}

{pstd}{ul:Working with equation estimation results}{p_end}

{pstd}Display information stored on learner AA e401.lrnAA
about the specification of conditional expectations for variable e401.{p_end}

{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(est_main)"}{p_end}
{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(stack_base_est)"}{p_end}

{pstd}Display information stored on results AA e401.resAA
about the estimation results for resamplings 1 and 2.{p_end}

{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(MSE_folds) key3(1)"}{p_end}
{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(MSE_folds) key3(2)"}{p_end}
{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(stack_weights) key3(1)"}{p_end}
{phang2}. {stata "ddml extract, vname(e401) key1(D1_pystacked) key2(stack_weights) key3(2)"}{p_end}

{pstd}{ul:Working directly with an equation associative array}{p_end}

{pstd}Extract the associative AA for the estimation of conditional expectations for variable e401.
Store it as a Mata object called AA_e401.
Note: the {cmd:crossfit} command returns an equation associative array,
so this step is unnecessary when using this command.{p_end}

{phang2}. {stata "ddml extract AA_e401, vname(e401)"}{p_end}
{phang2}. {stata "mata: AA_e401"}{p_end}

{pstd}Examples of working with this equation associative array.
Note that the {opt ename} option must be used.{p_end}

{phang2}. {stata "ddml extract, ename(AA_e401) key1(D1_pystacked) key2(MSE) key3(1)"}{p_end}
{phang2}. {stata "ddml extract, ename(AA_e401) key1(D1_pystacked) key2(MSE) key3(2)"}{p_end}

{pstd}{ul:Using Mata's associative array commands}{p_end}

{pstd}If preferred, Mata's associative array commands can be used directly.
Note that all keys are strings.{p_end}

{phang2}. {stata "mata: m0.estAA.keys()"}{p_end}
{phang2}. {stata `"mata: AA_e1_r2 = (m0.estAA).get(("ss","2"))"'}{p_end}
{phang2}. {stata "mata: AA_e1_r2.keys()"}{p_end}
{phang2}. {stata `"mata: AA_e1_r2.get(("b","post"))"'}{p_end}
