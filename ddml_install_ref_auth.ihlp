{smcl}
{marker installation}{title:Installation}

{pstd}
To verify that {opt ddml} and {opt qddml} are correctly installed, 
click on or type {stata "whichpkg ddml"} 
(which requires {helpb whichpkg} 
to be installed; {stata "ssc install whichpkg"}).

{pstd}
To get the latest stable versions of {opt ddml} and {opt qddml} from our website, 
check the installation instructions at {browse "https://statalasso.github.io/installation/"}.
We update the stable website version more frequently than the SSC version.


{marker references}{title:References}

{pstd}
Chernozhukov, V., Chetverikov, D., Demirer, M., 
Duflo, E., Hansen, C., Newey, W. and Robins, J. (2018), 
Double/debiased machine learning for 
treatment and structural parameters. 
{it:The Econometrics Journal}, 21: C1-C68. {browse "https://doi.org/10.1111/ectj.12097"}

{marker Hastie2009}{...}
{pstd}
Hastie, T., Tibshirani, R., & Friedman, J. (2009). 
The elements of statistical learning: data mining, inference,
and prediction. Springer Science & Business Media.

{marker Wolpert1992}{...}
{pstd}
Wolpert, David H. Stacked generalization. {it:Neural networks} 5.2 (1992): 241-259.
{browse "https://doi.org/10.1016/S0893-6080(05)80023-1"}


{marker authors}{title:Authors}

{pstd}
Achim Ahrens, Public Policy Group, ETH Zurich, Switzerland {break}
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


{marker also_see}{title:Also see (if installed)}

{pstd}
Help:
{helpb pystacked},
{helpb lasso2},
{helpb cvlasso},
{helpb rlasso},
{helpb ivlasso},
{helpb pdslasso}.{p_end}
