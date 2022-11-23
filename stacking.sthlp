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
{opt cvoos}
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
{synopt:{opt cvoos}}
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
{opt cvoos}
]

{pstd}
To get fitted values or residuals for each base learner:

{p 8 14 2}
{cmd:predict}
{it:type} {it:stub} 
[{cmd:if} {it:exp}] [{cmd:in} {it:range}]
{bind:[{cmd:,}}
{opt transf:orm}
{opt cvoos}
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
{synopt:{opt cvoos}}
Use cross-validated OOS predictions (default = use base learners re-fitted on full estimation sample).
{p_end}
{synoptline}

{pstd}
{it:Note:} Predicted values (in- and out-of-sample)
are calculated using the base learners re-fit on the full estimation sample.
To obtain cross-validated OOS predictions, use the {opt cvoos} option with {opt stacking}.

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


{marker examples}{...}
{title:Examples}

{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}

{pstd}Note that the variable created is called yhat_1 because the number of resamples defaults to 1.{p_end}
{phang2}. {stata "crossfit, estring(reg earnings $X) gen(yhat) kfolds(3)"}{p_end}
{phang2}. {stata "sum earnings yhat_1"}{p_end}

{pstd}As above but using 5 resamples.{p_end}
{phang2}. {stata "crossfit, estring(reg earnings $X) gen(yhat) kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "sum earnings yhat*"}{p_end}

{pstd}A simple example of 3-fold cross-validation with 5 resamples using {opt crossfit}.
The example uses {opt lasso2} from {opt lassopack}; click on {stata "ssc install lassopack"} to install.
We estimate using the following values of the lambda parameter: 2000, 1000, 500, 250.
Each time we call {opt crossfit} to obtain the residuals (prediction errors).
These could be used after cross-fitting to calculate the MSPE (mean squared prediction error),
but the MSPE is one of the returned results of {opt crossfit} so we just report that.
The specification that minimizes the MSPE for all 5 resamples is lambda=250.
{p_end}
{phang2}. {stata "crossfit, estring(lasso2 earnings $X, lglmnet lambda(2000)) gen(ehat2000) resid kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "mat list r(mse_list)"}{p_end}
{phang2}. {stata "crossfit, estring(lasso2 earnings $X, lglmnet lambda(1000)) gen(ehat1000) resid kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "mat list r(mse_list)"}{p_end}
{phang2}. {stata "crossfit, estring(lasso2 earnings $X, lglmnet lambda(500)) gen(ehat500) resid kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "mat list r(mse_list)"}{p_end}
{phang2}. {stata "crossfit, estring(lasso2 earnings $X, lglmnet lambda(250)) gen(ehat250) resid kfolds(3) reps(5)"}{p_end}
{phang2}. {stata "mat list r(mse_list)"}{p_end}


{marker results}{title:Saved results}

{p}{opt crossfit} saves the following results in {cmd:r()}:

Scalars
{col 4}{opt r(N)}{col 25}Number of observations.
{col 4}{opt r(mse)}{col 25}Mean squared prediction error in the last resample.

Macros
{col 4}{opt r(cmd_list)}{col 25}Estimation command

Matrices
{col 4}{opt r(N_list)}{col 25}Sample size; rows are resamples.
{col 4}{opt r(mse_list)}{col 25}MSPE; rows are resamples.
{col 4}{opt r(N_folds_list)}{col 25}Sample size by fold; rows are resamples.
{col 4}{opt r(mse_folds_list)}{col 25}MSPE by fold; rows are resamples.


{marker references}{title:References}

{phang}
Ahrens, A., Hansen, C.B. and M.E. Schaffer. 2020.
lassopack: model selection and prediction with regularized regression in Stata.
{it:The Stata Journal}, 20(1):176-235.
{browse "https://journals.sagepub.com/doi/abs/10.1177/1536867X20909697"}.
Working paper version: {browse "https://arxiv.org/abs/1901.05397"}.{p_end}

{phang}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

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
