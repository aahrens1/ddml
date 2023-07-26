{smcl}
{* *! {version 26july2023}{...}
{vieweralsosee "ddml main page" "ddml"}{...}
{vieweralsosee "Other" "ddml_stacking##also_see"}{...}
{hline}
{cmd:help ddml examples}{right: v1.4}
{hline}

{title:ddml examples for Double Debiased Machine Learning}

{pstd}
{opt ddml} implements algorithms for causal inference aided by supervised
machine learning as proposed in 
{it:Double/debiased machine learning for treatment and structural parameters}
(Econometrics Journal, 2018). Five different models are supported, allowing for 
binary or continuous treatment variables and endogeneity, high-dimensional 
controls and/or instrumental variables.

{pstd}
Below is a list of links to all {opt ddml} examples.
All the examples in the help files have clickable links and can be run by the user.

{pstd}{ul:Partially-linear model}:{p_end}
{pstd}{help ddml_example_partial_pystacked_basic:Basic example of the partially-linear model with pystacked}{p_end}
{pstd}{help ddml_example_partial_pystacked_detailed:Detailed example of the partially-linear model with pystacked}{p_end}
{pstd}{help ddml_example_partial_anylearner_detailed:Detailed general example of the partially-linear model with any learner(s)}{p_end}
{pstd}{help ddml_example_partial_pystacked_multitreat:Estimating the partially-linear model with multiple treatments}{p_end}

{pstd}{ul:Interactive model (ATE, ATET)}:{p_end}
{pstd}{help ddml_example_interactive_pystacked_basic:Basic example of the interactive model (ATE, ATET) with pystacked}{p_end}
{pstd}{help ddml_example_interactive_pystacked_detailed:Detailed example of the interactive model (ATE, ATET) with pystacked}{p_end}

{pstd}{ul:Partially-linear IV model}:{p_end}
{pstd}{help ddml_example_partialiv_pystacked_basic:Basic example of the partially-linear IV model with pystacked}{p_end}
{pstd}{help ddml_example_partialiv_anylearner_basic:Basic example of the partially-linear IV model with any learner(s)}{p_end}

{pstd}{ul:Flexible partially-linear IV model}:{p_end}
{pstd}{help ddml_example_flexiv_anylearner_basic:Basic example of the flexible partially-linear IV model with any learner(s)}{p_end}
{pstd}{help ddml_example_flexiv_anylearner_detailed:Detailed example of the flexible partially-linear IV model with any learner(s)}{p_end}

{pstd}{ul:Ineractive IV model (LATE)}:{p_end}
{pstd}{help ddml_example_interactiveiv_pystacked_basic:Basic example of the interactive IV model (LATE) with pystacked}{p_end}
{pstd}{help ddml_example_interactiveiv_pystacked_detailed:Detailed example of the interactive IV model (LATE) with pystacked}{p_end}

{pstd}{ul:Stacking regression with ddml}:{p_end}
{pstd}{help ddml_example_stacking:Detailed discussion of stacking with ddml plus examples}{p_end}

{pstd}{ul:ddml utilities}:{p_end}
{pstd}{help ddml_example_extract:Extracting stored information from ddml associative arrays}{p_end}
{pstd}{help ddml_example_describe:Describe the model setup and/or results}{p_end}
{pstd}{help ddml_example_export:Save estimated conditional expectations etc. to a csv file}{p_end}
{pstd}{help ddml_example_overlap:Overlap plots for interactive models}{p_end}

{pstd}{ul:Cluster cross-fitting with ddml}:{p_end}
{pstd}{help ddml_example_fcluster:Cluster sampling with cross-fit folds}{p_end}


{smcl}
INCLUDE help ddml_install_ref_auth
