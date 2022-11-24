{smcl}
{* *! version 21nov2022}{...}
{hline}
{cmd:help stacking}{right: v0.5}
{hline}


{title:stacking - Stata program for stacking regression}

{pstd}
{opt stacking} implements stacking regression ({helpb stacking##Wolpert1992:Wolpert, 1992}),
a way of combining multiple supervised
machine learners (the "base" or "level-0" learners) into a meta learner.
Any base learner with standard Stata {opt reg y x} syntax and with
a postestimation {opt predict, xb} command is supported.

{pstd}
Full syntax:

{p 8 14 2}
{cmd:stacking} depvar,
{opt estring(string)}
[{opt kfolds(integer)}
{opt foldvar(varname)}
{opt norandom}
{opt cv:alid}
{opt replace}
{opt noprefix}
{opt prefix(name)}
{opt tab:le}
{opt holdout}[{cmd:(}{it:varname}{cmd:)}]
{opt graph}[{cmd:(}{it:options}{cmd:)}]
{opt lgraph}[{cmd:(}{it:options}{cmd:)}]
{opt holdout}[{cmd:(}{it:varname}{cmd:)}]
]

{pstd}
Postestimation syntax:

{p 8 14 2}
{cmd:stacking},
{opt tab:le}
{opt holdout}[{cmd:(}{it:varname}{cmd:)}]
{opt graph}[{cmd:(}{it:options}{cmd:)}]
{opt lgraph}[{cmd:(}{it:options}{cmd:)}]
{opt holdout}[{cmd:(}{it:varname}{cmd:)}]
]


{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:General}
{synoptline}
{synopt:{opt estring(string)}}
A set of estimation commands separated by "||", e.g. "reg y x1 x2 || rlasso y x1 x2".
These are the base learners that will be "stacked" to obtain the ensemble prediction.
{p_end}
{synopt:{opt kfolds(integer)}}
Number of randomly drawn folds; ignored if {opt foldvar(varlist)} is specified; default=5.
{p_end}
{synopt:{opt foldvar(varname)}}
Integer variable with user-specified cross-fitting folds.
{p_end}
{synopt:{opt norandom}}
Use observations in existing order instead of randomizing before splitting into folds.
{p_end}
{synopt:{opt cv:alid}}
Save cross-validated out-of-sample predictions of each learner as Stata variables.
{p_end}
{synopt:{opt replace}}
Overwrite existing variables if they exist when saving CV OOS predictions.
{p_end}
{synopt:{opt noprefix}}
Don't add a prefix to fitted base learner and CV OOS predictions.
{p_end}
{synopt:{opt prefix(name)}}
Add the prefix {it:name} (default="_stack_") to fitted base learner and CV OOS predictions.
{p_end}

{synoptset 20}{...}
{synopthdr:Tables and graphs}
{synoptline}
{synopt:{opt holdout}[{cmd:(}{it:varname}{cmd:)}]}
Use a holdout sample for tables/graphs (default=all observations not in the estimation sample).
{p_end}
{synopt:{opt table}}
Report an MSPE table for the estimation sample and (if specified) the holdout sample.
{p_end}
{synopt:{opt graph}[{cmd:(}{it:options}{cmd:)}]}
Report graphs of in-sample or holdout sample observations;
{it:options} controls graphing options for the combined graph.
{p_end}
{synopt:{opt lgraph}[{cmd:(}{it:options}{cmd:)}]}
{it:options} controls graphing options for the individual learner graphs.
{p_end}


{title:Prediction}

{pstd}
To get predicted values or residuals:

{p 8 14 2}
{cmd:predict}
{it:type} {it:newname} 
[{cmd:if} {it:exp}] [{cmd:in} {it:range}]
{bind:[{cmd:,}}
{opt xb}
{opt resid}
{opt cv:alid}
]

{pstd}
To get fitted values or residuals for each base learner:

{p 8 14 2}
{cmd:predict}
{it:type} {it:stub} 
[{cmd:if} {it:exp}] [{cmd:in} {it:range}]
{bind:[{cmd:,}}
{opt transf:orm}
{opt cv:alid}
]

{synoptset 20}{...}
{synopthdr:Option}
{synoptline}
{synopt:{opt xb}}
Predicted value (default).
{p_end}
{synopt:{opt resid}}
Residuals.
{p_end}
{synopt:{opt transf:orm}}
Predicted values or residuals for each base learner.
The variable names are the {it:stub} with the number of the base learner appended.
{p_end}
{synopt:{opt cv:alid}}
Use cross-validated OOS predictions (default = use base learners re-fitted on full estimation sample).
{p_end}
{synoptline}

{pstd}
{it:Note:} Predicted values (in- and out-of-sample)
are calculated using the base learners re-fit on the full estimation sample.
To obtain cross-validated OOS predictions, use the {opt cvalid} option with {opt stacking}.

{pstd}
{it:Note:} Predicted values of re-fitted base learners are automatically created by {opt stacking}
and their names are saved in the macro {opt e(base_yhat)}.
Postestimation support for (re-)creating these predicted values
is to support standard Stata {opt predict} syntax.


{marker summary}{...}
{title:Summary}

{pstd}
{opt stacking} implements stacking regression ({helpb stacking##Wolpert1992:Wolpert, 1992}),
an ensemble method that combines multiple supervised
machine learners (the "base" or "level-0" learners) into a meta learner.
Each learner is fit on K-1 folds and the cross-validated out-of-sample (OOS) predictions
are obtained for the Kth fold.
This process is repeated so that each fold serves once as the holdout fold
for which OOS predictions created.
After all learners have been fit on all splits,
a full set of OOS predictions is available for each learner.
The final "stacked" predictor is obtained as a weighted average of the predictions of each learner,
where in this case the predictions are obtained by re-fitting the learners on the full sample.
The weights are obtained by fitting the outcome variable on the OOS prediction using
non-negative least squares (NNLS) without an intercept and with the constraints
that weights sum to one and individual learner weights lie in the interval [0,1].
{p_end}

{pstd}
{opt stacking} reports the NNLS weights
and optionally saves the OOS predictions of the individual learners.
The ensemble prediction is obtained postestimation using {opt predict}.
{p_end}

{pstd}
Any base learner with standard Stata {opt reg y x} syntax and with
a postestimation {opt predict, xb} command is supported.
{p_end}

{pstd}
{opt stacking}'s default is to generate a single random split into folds.
This can be overridden by specifying user-defined fold variables,
or by the {opt norandom} option (indicating that the split use the data in the existing order).
{p_end}

{pstd}
After estimation, {opt stacking} can report a table of in-sample
and, optionally, out-of-sample (holdout sample) MSPE (mean squared prediction errors)
for both the stacking regression and the base learners.
The default holdout sample used for out-of-sample performance with the {opt holdout} option
is all observations not included in the estimation.
Alternatively, the user can specify the holdout sample explicitly
using the syntax {opt holdout(varname)}.
The table can be requested postestimation or as part of the {opt stacking} estimation command.

{pstd}
{opt stacking} can also report graphs of in-sample
and, optionally, out-of-sample (holdout sample) performance
for both the stacking regression and the base learners.
The graphs compare predicted vs actual values of {it:depvar}.
As with the {opt table} option, the default holdout sample used for out-of-sample performance
is all observations not included in the estimation,
but the user can instead specify the holdout sample explicitly.
The table can be requested postestimation or as part of the {opt stacking} estimation command.

{pstd}
The {opt graph} option on its own reports the graphs using {opt pystacked}'s default settings.
Because graphs are produced using Stata's {helpb twoway} command,
the user can control either the combined graph ({opt graph(options)})
or the individual learner graphs ({opt lgraph(options)}) appear by passing options to these commands.

{pstd}
{helpb pystacked} is a program by the same authors as {opt stacking}
that supports stacking regression via Python and
{browse "https://scikit-learn.org/stable/index.html":scikit-learn}'s 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingRegressor.html":sklearn.ensemble.StackingRegressor} and 
{browse "https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.StackingClassifier.html":sklearn.ensemble.StackingClassifier}. 
The main differences are that {opt stacking} will support any Stata estimator with standard syntax
({helpb pystacked} is limited to Python learners supported by sklearn)
but {helpb pystacked} is typically faster.

{pstd}
{helpb pystacked} is useful here because it can be used as a front-end for any single learner supported by sklearn.
This means {opt stacking} can be used to stack sklearn learners with learners from other packages,
as long as they conform to standard Stata syntax.
An example is provided below.


{marker examples}{...}
{title:Examples}

{pstd}Load housing data.{p_end}
{phang2}. {stata "insheet using https://statalasso.github.io/dta/housing.csv, clear"}

{pstd}
Stacking regression with regress and the "rigorous" (plug-in) lasso
available from the package {helpb lassopack}.
Default is 5-fold cross-validation.
Request cross-validated predictions with the {opt cvalid} option;
use {opt replace} if they already exist in memory.
{p_end}
{phang2}. {stata "stacking medv, estring(regress medv crim-lstat || rlasso medv crim-lstat) cv replace"}{p_end}

{pstd}
The weights determine how much each base learner contributes
to the final stacking prediction.{p_end}

{pstd}
Request the MSPE table:{p_end}
{phang2}. {stata "stacking, table"}{p_end}

{pstd}
Re-estimate using the first 400 observations, and request the MSPE table.
MSPEs for in-sample, cross-validated and the default holdout sample (all unused observations) are reported.
Use the {opt norandom} to force a split based on the existing order of observations:{p_end}
{phang2}. {stata "stacking medv if _n<=400, estring(regress medv crim-lstat || rlasso medv crim-lstat) cv replace norandom"}{p_end}
{phang2}. {stata "stacking, table holdout"}{p_end}

{pstd}
Graph predicted vs actual for the holdout sample:{p_end}
{phang2}. {stata "stacking, graph holdout"}{p_end}

{pstd}
Storing the predicted values:{p_end}
{phang2}. {stata "predict double yhat, xb"}{p_end}

{pstd}
We can also save the predicted values of each base learner:{p_end}
{phang2}. {stata "predict double yhat, transform"}{p_end}

{pstd}
{opt stacking} vs {opt pystacked}.
The first example uses {opt pystacked} to do both base learner estimation and stacking.
The second example uses separate calls to {opt pystacked} to estimate the base learners,
and {opt pystacked} does the stacking.
The stacking weights are essentially the same.

{phang2}. {stata "pystacked medv zn-rad if _n<=200, type(regress) methods(ols lassoic) folds(2)"}{p_end}

{phang2}. {stata "stacking medv if _n<=200, estring(pystacked medv zn-rad, type(regress) methods(ols) || pystacked medv zn-rad, type(regress) methods(lassoic)) kfolds(2) norandom"}{p_end}


{marker results}{title:Saved results}

{p}{opt stacking} saves the following results in {cmd:e()}:

Scalars
{col 4}{opt e(N)}{col 25}Number of observations.
{col 4}{opt e(mcount)}{col 25}Number of base learners.
{col 4}{opt e(cvalid)}{col 25}=1 if cross-validated predictions were created, =0 if not.

Macros
{col 4}{opt e(estring)}{col 25}Estimation string with base learner estimation commands.
{col 4}{opt e(estring_1)}{col 25}Estimation string for base learner 1.
{col 4}{opt e(estring_2)}{col 25}Estimation string for base learner 2.
{col 4}...
{col 4}{opt e(base_est)}{col 25}List of base learners, prefixed by Y1_, Y2_, ....
{col 4}{opt e(base_cv)}{col 25}Varlist of cross-validated predictions of base learners.
{col 4}{opt e(base_yhat)}{col 25}Varlist of re-estimated (full sample) predictions of base learners.
{col 4}{opt e(depvar)}{col 25}Dependent (outcome) variable.

Matrices
{col 4}{opt e(weights)}{col 25}Stacking weights by base learner.
{col 4}{opt e(rmspe)}{col 25}Root MSPEs by base learner - in-sample, CV and holdout (created by {opt table} option).
{col 4}{opt e(N_list)}{col 25}Sample size by base learner.
{col 4}{opt e(N_folds_list)}{col 25}Sample size by base learner and fold.
{col 4}{opt e(mse_list)}{col 25}Cross-validated MSPEs by base learner.
{col 4}{opt e(mse_folds_list)}{col 25}Cross-validated MSPEs by base learner and fold.


{marker references}{title:References}

{phang}
Ahrens, A., Hansen, C.B. and M.E. Schaffer. 2020.
lassopack: model selection and prediction with regularized regression in Stata.
{it:The Stata Journal}, 20(1):176-235.
{browse "https://journals.sagepub.com/doi/abs/10.1177/1536867X20909697"}.
Working paper version: {browse "https://arxiv.org/abs/1901.05397"}.{p_end}

{phang}
Ahrens, A., Hansen, C.B. and M.E. Schaffer. 2022.
Stacking generalization and machine learning in Stata.
Working paper version: {browse "https://arxiv.org/abs/2208.10896"}.{p_end}

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}

{title:Authors}

{pstd}
Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland  {break}
achim.ahrens@gess.ethz.ch

{pstd}
Christian B. Hansen, University of Chicago, USA {break}
Christian.Hansen@chicagobooth.edu

{pstd}
Mark E Schaffer, Heriot-Watt University, UK {break}
m.e.schaffer@hw.ac.uk   

{pstd}
Thomas Wiemann, University of Chicago, USA {break}
wiemann@uchicago.edu


{title:Also see (if installed)}

{pstd}
Help: {helpb ddml}, {helpb qddml}, {helpb lasso2}, {helpb cvlasso}.{p_end}
