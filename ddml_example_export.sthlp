{smcl}
{* *! version 27july2023}{...}
{smcl}
{pstd}{ul:ddml export utility - Basic example with {help pystacked}}{p_end}

{pstd}Load the data, define global macros, set the seed and initialize the model.
Use 2-fold cross-fitting with two repetitions (resamples)
Use {help pystacked}'s default learners as the supervised learners.
We explicitly name the model to be estimated as m_sip1991;
this is the name of the Mata global containing the model details,
and will also prefix the sample and fold indicators created by {opt ddml}.
{p_end}

{phang2}. {stata "use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear"}{p_end}
{phang2}. {stata "global Y net_tfa"}{p_end}
{phang2}. {stata "global D e401"}{p_end}
{phang2}. {stata "global X tw age inc fsize educ db marr twoearn pira hown"}{p_end}
{phang2}. {stata "set seed 42"}{p_end}
{phang2}. {stata "ddml init partial, kfolds(2) reps(2) mname(m_sip1991)"}{p_end}
{phang2}. {stata "ddml E[Y|X], mname(m_sip1991): pystacked $Y $X"}{p_end}
{phang2}. {stata "ddml E[D|X], mname(m_sip1991): pystacked $D $X"}{p_end}
{phang2}. {stata "ddml crossfit, mname(m_sip1991)"}{p_end}
{phang2}. {stata "ddml estimate, mname(m_sip1991)"}{p_end}

{pstd}It will often be a good idea to include an ID variable that identifies the observation number.
This dataset doesn't include an ID variable, so we create one.{p_end}

{phang2}. {stata "gen long m_sip1991_id = _n"}{p_end}

{pstd}To include the ID variable with everything else, we use the {opt addvars(.)} option.
The data will be exported to a CSV file called "m_ddml_sip1991.csv".{p_end}

{phang2}. {stata "ddml export using m_ddml_sip1991.csv, mname(m_sip1991) replace addvars(m_sip1991_id)"}{p_end}
