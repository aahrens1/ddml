// wrapper to execute all ddml-related cert scripts

adopath ++ "/Users/kahrens/MyProjects/ddml"
cd "/Users/kahrens/MyProjects/ddml/cert"

// general
do ddml_cert

// detailed
do ddml_cert_crossfit
do ddml_cert_fiv
do ddml_cert_interactive
do ddml_cert_interactiveiv
do ddml_cert_misc
do ddml_cert_partial
do ddml_cert_partial_iv

// other
do ddml_cert_helpfiles
do qddml_cert
