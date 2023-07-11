{smcl}
{* *! version 11july2023}{...}
{hline}
{cmd:help ddml overlap}{right: v1.3}
{hline}

{title:ddml overlap commands for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables. 

{pstd}
{cmd:ddml overlap} reports overlap plots following estimation of the {opt ddml} { {it:interactive} and {it:interactiveiv} models.
One overlap (line) plot of propensity scores is reported for each treatment variable learner;
by default, propensity scores for all crossfit samples are plotted.
Overlap plots for the treatment variables are combined using {helpb graph combine}.

{marker syntax}{...}
{title:Syntax}

{p 8 14}{cmd:ddml overlap} [ {opt mname(name)} {opt replist(numlist)} {opt pslist(namelist)} {opt n(integer)} {opt kernel(name)}
{opt name(name [, replace])} {opt title(string)} {opt subtitle(string)} {opt lopt0(string)}
{opt lopt1(string)}{bind: ]}

{synoptset 20}{...}
{synopthdr:Options}
{synoptline}
{synopt:{opt mname(name)}}
name of the DDML model. Defaults to {it:m0}.
{p_end}
{synopt:{opt replist(numlist)}}
list of crossfitting resamples to plot. Defaults to all.
{p_end}
{synopt:{opt pslist(namelist)}}
varnames of propensity scores to plot (excluding the resample number). Defaults to all.
{p_end}
{synopt:{opt n(integer)}}
see {helpb teffects overlap}.
{p_end}
{synopt:{opt kernel(name)}}
see {helpb teffects overlap}.
{p_end}
{synopt:{opt name(name)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt title(string)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt subtitle(string)}}
see {helpb graph combine}.
{p_end}
{synopt:{opt lopt0(string)}}
options for line plot of untreated; default is solid/navy; see {helpb line}.
{p_end}
{synopt:{opt lopt0(string)}}
options for line plot of treated; default is short dash/dark orange; see {helpb line}.
{p_end}
{synoptline}

{p2colreset}{...}
{pstd}

{title:Examples}

For examples of usage see {help ddml##examples:help ddml}.


{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}

{marker installation}{title:Installation}

{pstd}
To get the latest stable version of {cmd:ddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.

{pstd}
To verify that {cmd:ddml} is correctly installed, 
click on or type {stata "whichpkg ddml"} 
(which requires {helpb whichpkg} 
to be installed; {stata "ssc install whichpkg"}).

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
Help: {helpb lasso2}, {helpb cvlasso}, {helpb rlasso}, {helpb ivlasso},
 {helpb pdslasso}, {helpb pystacked}.{p_end}
