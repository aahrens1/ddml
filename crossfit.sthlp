{smcl}
{* *! version 15sep2022}{...}
{hline}
{cmd:help ddml}{right: v0.5}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{hi: crossfit} {hline 2}}Stata program for Double Debiased Machine Learning{p_end}
{p2colreset}{...}

{pstd}
{opt crossfit} fits a supervised machine learner on K-1
random folds and returns the out-of-sample predicted 
values (or residuals) for the hold-out fold. This is done iteratively
to obtain out-of-sample ("cross-fitted") fitted values (or residuals)
for the whole sample.

{p 8 14 2}
{cmd:crossfit} , 
{opt estring(string)}
{opt vtilde(vtilde)}
[{opt foldvar(varname)}
{opt kfolds(integer)}
{opt reps(integer)}
{opt predopt(string)}
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
an estimation string, e.g. "reg y x1 x2", that will be 
repeatedly invoked. See note on compatible programs 
{helpb ddml##compatibility:here}.
{p_end}
{synopt:{opt vtilde(vtilde)}}
name of the new variable to be created.
{p_end}
{synopt:{opt foldvar(varname)}}
integer variable with user-specified cross-fitting folds.
{p_end}
{synopt:{opt kfolds(integer)}}
number of randomly drawn folds (ignored if {opt foldvar(varname)} is specified).
{p_end}
{synopt:{opt reps(integer)}}
number of re-sampling iterations, i.e., how often the cross-fitting procedure is
repeated on randomly generated folds (ignored if {opt foldvar(varname)} is specified).
{p_end}
{synopt:{opt predopt(string)}}
{cmd:predict} option to be used to get predicted values. 
Typical values could be {opt xb} or {opt pr}. Default is 
blank.
{p_end}
{synopt:{opt vtype(string)}}
variable type of the variable to be created. Defaults to {it:double}. 
{it:none} can be used to leave the type field blank.
{p_end}

{marker compatibility}{...}
{title:Compatible programs}

{pstd} 
See {helpb ddml##compatibility:here}.

{marker examples}{...}
{title:Examples}

{phang2}. {stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta,clear"}{p_end}
{phang2}. {stata "global Y earnings"}{p_end}
{phang2}. {stata "global X sex age married black hispanic"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "crossfit, estring(reg $Y $X) vtilde(test) kfolds(3)"}{p_end}

{marker references}{title:References}

{pstd}
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
Help: {helpb ddml}, {helpb qddml}, {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}, {helpb pystacked}.{p_end}
