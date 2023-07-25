// do file to run all clickable examples in ddml help files
// this file requires the Stata program loghelp
// install from SSC if available or from github:
// net install loghelp, ///
//   from(https://raw.githubusercontent.com/markeschaffer/stata-utilities/master) ///
//   replace

clear all

cap log close
log using "ddml_helpfiles_cert", replace smcl
which ddml
mata: whichddml()
which qddml
which crossfit
which pystacked
log close

foreach sthlpfile in									///
	ddml_example_partial_pystacked_basic.sthlp			///
	ddml_example_partial_pystacked_detailed.sthlp		///
	ddml_example_partial_anylearner_detailed.sthlp		///
	ddml_example_partial_pystacked_multitreat.sthlp		///
	ddml_example_interactive_pystacked_basic.sthlp		///
	ddml_example_interactive_pystacked_detailed.sthlp	///
	ddml_example_partialiv_pystacked_basic.sthlp		///
	ddml_example_partialiv_anylearner_basic.sthlp		///
	ddml_example_flexiv_anylearner_basic.sthlp			///
	ddml_example_flexiv_anylearner_detailed.sthlp		///
	ddml_example_interactiveiv_pystacked_basic.sthlp	///
	ddml_example_interactiveiv_pystacked_detailed.sthlp	///
	ddml_example_extract.sthlp							///
	ddml_example_stacking.sthlp							///
	ddml_example_describe.sthlp							///
	ddml_example_export.sthlp							///
	ddml_example_overlap.sthlp							///
	ddml_example_fcluster.sthlp							///
	qddml.sthlp											///
	crossfit.sthlp										///
	{

	log using "ddml_helpfiles_cert", append smcl
	di
	di "Help file: `sthlpfile'"
	di
	log close

	findfile `sthlpfile'
	loghelp using ddml_helpfiles_cert,					///
		inputname(`r(fn)') append smcl
	
	}
