{smcl}
{* *! version 21jan2023}{...}
{hline}
{cmd:help crossfit}{right: v1.2}
{hline}


{title:crossfit - Stata program for cross-fitting}

{pstd}
{opt crossfit} fits a supervised machine learner on K-1 folds
and returns the out-of-sample predicted values (or residuals) for the holdout fold.
This is done iteratively to obtain out-of-sample ("cross-fitted") fitted values (or residuals)
for the whole sample.

{p 8 14 2}
{cmd:crossfit} , 
{opt estring(string)}
{opt g:enerate(stubname)}
[{opt kfolds(integer)}
{opt foldvar(varlist)}
{opt norandom}
{opt reps(integer)}
{opt vtype(string)}]

{pstd}
{opt crossfit} is an auxiliary program that is internally used by 
{helpb ddml} and {helpb qddml}, but can be used for other purposes.


{marker syntax}{...}
{title:Options}

{synoptset 20}{...}
{synopthdr:General}
{synoptline}
{synopt:{opt estring(string)}}
An estimation string, e.g. "reg y x1 x2", that will be 
repeatedly invoked. See note on compatible programs 
{helpb ddml##compatibility:here}.
{p_end}
{synopt:{opt g:enerate(stubname)}}
Name of the new variable to be created;
the resample number is appended to the end of the variable name.
Note that if the variable (including the resample number) already exists, it is overwritten.
{p_end}
{synopt:{opt resid}}
Create residuals (prediction errors); default is predicted values of outcome variable.
{p_end}
{synopt:{opt kfolds(integer)}}
Number of randomly drawn folds; ignored if {opt foldvar(varlist)} is specified; default=5.
{p_end}
{synopt:{opt foldvar(varlist)}}
Integer variable(s) with user-specified cross-fitting folds; one foldvar per resample.
{p_end}
{synopt:{opt norandom}}
Use observations in existing order instead of randomizing before splitting into folds;
if multiple resamples, applies to first resample only;
ignored if user-defined fold variables are provided in {opt foldvar(varlist)}.
{p_end}
{synopt:{opt reps(integer)}}
Number of resampling iterations, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds;
ignored if {opt foldvar(varlist)} is specified;
default=1.
{p_end}
{synopt:{opt vtype(string)}}
Variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank.
{p_end}


{marker summary}{...}
{title:Summary}

{pstd}
{opt crossfit} fits a supervised machine learner on K-1 folds
and returns the out-of-sample predicted values or residuals for the holdout fold.
This process is repeated so that each fold serves once as the holdout fold
for which predictions or residuals are created.
At the end of the cross-fitting, a full set of predictions or residuals is available
in the new variable specified by the {opt generate} option.
The "supervised machine learner" can be any Stata estimator
that supports standard postestimation prediction.
{p_end}

{pstd}
{opt crossfit}'s default is to generate a single random split into folds.
This can be overridden by specifying user-defined fold variables,
or by the {opt norandom} option (indicating that the split use the data in the existing order).
{p_end}

{pstd}
{opt crossfit} allows multiple resampling,
meaning that the procedure is applied repeatedly
using multiple fold variables that indicate different fold splits.
This can be done via the {opt reps} option,
or by providing multiple user-defined fold variables.
The resample number is appended to the generated predictions/residuals.
{p_end}

{pstd}
The output of {opt crossfit} can be seen as the intermediate step
of standard K-fold cross-validation.
In a typical cross-validation exercise, a search is conducted across a range of specifications (e.g. values for a tuning parameter).
The prediction errors (residuals) for the holdout folds are assembled for each specification,
and the specification with the best prediction performance (e.g. smallest mean squared prediction error) is chosen.
A simple example of how to use {opt crossfit} to do this is below.
{p_end}


{marker compatibility}{...}
{title:Compatible programs}

{pstd} 
See {help ddml##compatibility:here}.


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
