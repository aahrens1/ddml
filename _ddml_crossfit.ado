*! ddml v1.5.0
*! last edited: 18dec2025
*! authors: aa/ms

*** ddml cross-fitting
program _ddml_crossfit, eclass sortpreserve
	version 16
	syntax [anything] ,								/// 
							[						///
							debug					/// 
							mname(name)				///
							shortstack				///
							poolstack				///
							NOSTDstack				/// no standard stacking (pystacked only); use voting to get learners
							NOIsily					///
							Verbose 				///
							crossfitother			/// force use of general (not-pystacked-specific) code
							NOENFORCElie			/// do not enforce LIE in fiv model
							*						/// other options
							]

	local debugflag		= "`debug'"~=""
	if "`noisily'"==""	local qui qui
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens((`mname'.nameD)))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("reps",strofreal(`mname'.nreps))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'

	// reps = total number of reps; crossfitted = reps done so far (=0 if none)
	mata: st_local("reps", strofreal(`mname'.nreps))
	mata: st_local("crossfitted", strofreal(`mname'.crossfitted))

	// clear equation and estimation results from model struct if starting over (reps=crossfitted)
	if `reps'==`crossfitted' {
		mata: clear_model_results(`mname')
		local crossfitted = 0
	}
	else {
		// just clear any preexisting estimation results from the model struct
		mata: clear_model_estimation(`mname')
	}
	
	// standard stacking for model estimation (all eqns) possible only if pystacked used for every eqn
	mata: st_local("stdflag", strofreal(`mname'.stdflag))
	
	// set flags
	local ssflag		= "`shortstack'"~=""
	local psflag		= "`poolstack'"~=""
	// update stdflag if nostdstack specified
	local stdflag		= "`nostdstack'"=="" & `stdflag'
	// fiv model only - update on the model struct
	if "`model'"=="fiv" & "`noenforcelie'"=="" {
		mata: `mname'.lieflag = 1
	}
	else {
		mata: `mname'.lieflag = 0
	}

	mata: st_local("prefixflag",strofreal(`mname'.prefixflag))
	if `prefixflag'		local prefix `mname'_
	local touse `mname'_sample
	
	local firstrep	= `crossfitted'+1
	local lastrep	= `reps'
	// fold IDs
	forvalues m=`firstrep'/`lastrep' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	`qui' di as text "Fold IDs: `fidlist'"
	
	// check that required learners have been added
	if "`nameY'"=="" {
		di as err "error - no Y learner defined"
		exit 198
	}
	if `numeqnD'==0 {
		di as err "error - no D learner defined"
		exit 198
	}
	if `numeqnZ'==0 & ("`model'"=="interactiveiv" | "`model'"=="iv") {
		di as err "error - no Z learner defined"
		exit 198
	}
	
	// check interactive model
	if "`model'"=="interactive" {
		if `numeqnD' > 1 {
			di as err "error - interactive model allows only one D variable"
			exit 198
		}
		qui count if `nameD'~=1 & `nameD'~=0 & `touse'
		if r(N) > 0 {
			di as err "error - interactive model supported only for D=0 or D=1"
			exit 198
		}
	}
	
	// check interactive IV (LATE) model
	if "`model'"=="interactiveiv" {
		if `numeqnD' > 1 {
			di as err "error - interactiveiv model allows only one D variable"
			exit 198
		}
		if `numeqnZ' > 1 {
			di as err "error - interactiveiv model allows only one Z variable"
			exit 198
		}
		qui count if `nameZ'~=1 & `nameZ'~=0 & `touse'
		if r(N) > 0 {
			di as err "error - interactiveiv model supported only for Z=0 or Z=1"
			exit 198
		}
		qui count if `nameD'==1 & `nameZ'==0 & `touse'
		if r(N)>0 {
			di as text "note: treatment (`nameD') = 1 in `r(N)' cases when assignment (`nameZ') = 0"
		}
	}
	
	// equations and learners
	
	// will be set to zero if any eqn doesn't use pystacked,
	// or if any equation has pystacked with a single learner
	
	// will always be a Y eqn
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	`qui' di as text "Y eqn learners: `vtlistY'"
	local vtlist `vtlistY'
	
	// will always be a D eqn
	`qui' di as text "D equations (`numeqnD'): `nameD'"
	foreach var of varlist `nameD' {
		`qui' di as text _col(5) "D equation `var':"
		mata: `eqn' = (`mname'.eqnAA).get("`var'")
		mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
		`qui' di as text _col(10) "learners: `vtlistD'"
		local vtlist `vtlist' `vtlistD'
	}
	
	// Z eqn exists for late, iv models
	if `numeqnZ' {
		`qui' di as text "Z equations (`numeqnZ'): `nameZ'"
		foreach var of varlist `nameZ' {
			`qui' di as text _col(5) "Z equation `var':"
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			`qui' di as text _col(10) "learners: `vtlistZ'"
			local vtlist `vtlist' `vtlistZ'
		}
	}

	// update model struct flags for stacking
	if `stdflag'	mata: `mname'.stdflag	= 1
	else 			mata: `mname'.stdflag	= 0
	if `ssflag'		mata: `mname'.ssflag	= 1
	else 			mata: `mname'.ssflag	= 0
	if `psflag' 	mata: `mname'.psflag	= 1
	else			mata: `mname'.psflag	= 0

	
	************************** Y equation **************************
	// will always be a Y eqn
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	// used to track minimum number of learners in an equation; must be >1 for short/pool stacking
	mata: st_local("numlnrY", strofreal(`eqn'.nlearners))
	local minlearners = `numlnrY'
	// flag for integrated pystacked code
	mata: st_local("pystacked_Y", strofreal(`eqn'.pystackedmulti))
	
	// shortstack and poolstack variable names
	if `ssflag'		mata: `eqn'.shortstack = "`prefix'Y_`nameY'"
	else			mata: `eqn'.shortstack = ""
	if `psflag'		mata: `eqn'.poolstack = "`prefix'Y_`nameY'"
	else			mata: `eqn'.poolstack = ""

	// set/clear treatvar macro
	if "`model'"=="interactive" {
		local treatvar	`nameD'
	}
	else if ("`model'"=="interactiveiv") {
		local treatvar	`nameZ'
	}
	else if ("`model'"=="fiv") {
		local treatvar
	} 
	else {
		local treatvar
	}
	
	if ("`model'"=="interactive")	di as text "Cross-fitting E[y|X,D] equation: `nameY'"
	else							di as text "Cross-fitting E[y|X] equation: `nameY'"

	if `ssflag' & `numlnrY'==1 & `pystacked_Y'<=1 {
		di as text "Note: short-stacked fitted values = fitted values from single learner"
	}
	
	if `psflag' & (`pystacked_Y'==0) {
		di as text "Note: poolstack option requires pystacked as the single learner; option ignored"
		mata: `eqn'.poolstack = ""
	}
	else if `psflag' & `numlnrY'==1 & `pystacked_Y'<=1 {
		di as text "Note: pool-stacked fitted values = fitted values from single learner"
	}
	
	crossfit if `touse',						///
		ename(`eqn') noreplace					///
		pystackedmulti(`pystacked_Y')			///
		foldvar(`fidlist')						///
		firstrep(`firstrep')					///
		treatvar(`treatvar')					///
		model(`model')							///
		`nostdstack'							///
		`crossfitother'							///
		`noisily'								///
		`options'
	// resinsert into model struct AA with equations
	mata: (`mname'.eqnAA).put("`nameY'",`eqn')

	************************** D equation **************************
	// will always be a D eqn
	if ("`model'"=="interactiveiv") {
		local treatvar	`nameZ'
	}
	else {
		// clear local
		local treatvar
	}
	foreach var of varlist `nameD' {
		mata: `eqn' = (`mname'.eqnAA).get("`var'")
		mata: st_local("numlnrD",strofreal(cols(`eqn'.vtlist)))
		// flag for integrated pystacked code
		mata: st_local("pystacked_D", strofreal(`eqn'.pystackedmulti))

		// shortstack and poolstack variable names
		if `ssflag'		mata: `eqn'.shortstack = "`prefix'D_`var'"				
		else			mata: `eqn'.shortstack = ""
		if `psflag'		mata: `eqn'.poolstack = "`prefix'D_`var'"				
		else			mata: `eqn'.poolstack = ""
		if ("`model'"=="partial")						di as text "Cross-fitting E[D|X] equation: `var'"
		if ("`model'"=="interactiveiv")					di as text "Cross-fitting E[D|X,Z] equation: `var'"
		if ("`model'"=="fiv")							di as text "Cross-fitting E[D|X,Z] and E[D|X] equation: `var'"
		if ("`model'"=="interactive"|"`model'"=="iv")	di as text "Cross-fitting E[D|X] equation: `var'"
		
		if `ssflag' & `numlnrD'==1 & `pystacked_D'==1 {
			di as text "Note: short-stacked fitted values = fitted values from single learner"
		}
		if `psflag' & (`pystacked_D'==0) {
			di as text "Note: poolstack option requires pystacked as the single learner; option ignored"
			mata: `eqn'.poolstack = ""
		}
		else if `psflag' & `numlnrD'==1 & `pystacked_D'<=1 {
			di as text "Note: pool-stacked fitted values = fitted values from single learner"
		}

		// All learners for each D eqn
		crossfit if `touse',						///
			ename(`eqn') noreplace					///
			pystackedmulti(`pystacked_D')			///
			foldvar(`fidlist')						///
			firstrep(`firstrep')					///
			treatvar(`treatvar')					///
			model(`model')							///
			`noenforcelie'							/// applicable only to D eqns
			`nostdstack'							///
			`crossfitother'							///
			`noisily'								///
			`options'
		if ("`model'"=="interactiveiv") {
			// LATE special case - perfect assignment to treatment
			if r(perfectflag)==1 {
				mata: `mname'.perfectflag=1
			}
		}
		if ("`model'"=="fiv" & "`r(fiv_warn)'"~="") {
			// LIE + shortstacking + pystacked integration special case
			mata: `mname'.lieflag_ss = 0
		}
		mata: (`mname'.eqnAA).put("`var'",`eqn')
	}

	************************** Z equation **************************
	// iv and interactiveiv only
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("numlnrZ",strofreal(cols(`eqn'.vtlist)))
			// flag for integrated pystacked code
			mata: st_local("pystacked_Z", strofreal(`eqn'.pystackedmulti))
			// shortstack and poolstack variable names
			if `ssflag'		mata: `eqn'.shortstack = "`prefix'Z_`var'"				
			else			mata: `eqn'.shortstack = ""
			if `psflag'		mata: `eqn'.poolstack = "`prefix'Z_`var'"				
			else			mata: `eqn'.poolstack = ""
			di as text "Cross-fitting E[Z|X]: `var'"
			if `ssflag' & `numlnrZ'==1 & `pystacked_Z'<=1 {
				di as text "Note: short-stacked fitted values = fitted values from single learner"
			}
			if `psflag' & (`pystacked_Z'==0) {
				di as text "Note: poolstack option requires pystacked as the single learner; option ignored"
				mata: `eqn'.poolstack = ""
			}
			else if `psflag' & `numlnrZ'==1 & `pystacked_Z'<=1 {
				di as text "Note: pool-stacked fitted values = fitted values from single learner"
			}
			
			// All learners for each Z eqn
			crossfit if `touse',						///
				ename(`eqn') noreplace					///
				pystackedmulti(`pystacked_Z')			///
				foldvar(`fidlist')						///
				firstrep(`firstrep')					///
				model(`model')							///
				`nostdstack'							///
				`crossfitother'							///
				`noisily'								///
				`options'
			mata: (`mname'.eqnAA).put("`var'",`eqn')
		}
	}
	
	// create sample dfn variable by resample
	create_sample_indicators, mname(`mname')
	
	// set flag on model struct
	// reps=resamplings to be done, crossfitted=resamplings done so far
	// if appending, crossfitted updated to complete set of reps
	mata: `mname'.crossfitted = `mname'.nreps

	// report results by equation type with resamplings grouped together
	if ("`verbose'"!="") {
		di
		di as res "Reporting crossfitting results:"
		_ddml_describe `mname'
	}
	
	// drop temp eqn struct	
	mata: mata drop `eqn'
	
end

// creates sample indicators and leaves them in memory
// same naming as main sample variable but with a _m resample subscript
program create_sample_indicators
	version 16
	syntax [anything], mname(name)
	
	*** extract details of estimation
	mata: model_chars(`mname')
	local nreps		= r(nreps)
	local numeqnD	= r(numeqnD)
	local numeqnZ	= r(numeqnZ)

	// collect names of Y variables
	local vlist `r(Y)' `r(Y_L)'

	// collect names of D variables
	forvalues i=1/`numeqnD' {
		local vlist `vlist' `r(D`i')' `r(D`i'_L)' `r(D`i'_h)'
	}

	// collect names of Z variables
	forvalues i=1/`numeqnZ' {
		local vlist `vlist' `r(Z`i')' `r(Z`i'_L)'
	}
	
	forvalues m=1/`nreps' {
		cap drop `mname'_sample_`m' 
		qui gen byte `mname'_sample_`m' = `mname'_sample
		label var `mname'_sample_`m' "Sample indicator for rep `m'"
		foreach vt in `vtlist' {
			// check that variable exists (may not if e.g. no std stacking with pystacked)
			cap confirm variable `vt', exact
			if _rc==0	qui replace `mname'_sample_`m' = 0 if `vt'_`m'==.
		}
	}
	
	// create mean/median sample indicators
	if `nreps'>1 {
		cap drop `mname'_sample_mn
		qui gen byte `mname'_sample_mn = `mname'_sample_1
		forvalues m=2/`nreps' {
			qui replace `mname'_sample_mn = `mname'_sample_`m' if `mname'_sample_mn==0 & `mname'_sample_`m'==1
		}
		cap drop `mname'_sample_md
		qui gen byte `mname'_sample_md = `mname'_sample_mn
	}
	
end

program report_debugging
	version 16
	syntax name(name=mname), [ fidlist(varlist) ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()

	// locals used below
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'

	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	
	local var `nameY'
	foreach vtilde in `vtlist' {
		mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
		mata: st_local("est_main", return_learner_item(`eqn',"`vtilde'","est_main"))
		mata: st_local("est_options", return_learner_item(`eqn',"`vtilde'","est_options"))
		di
		di as res "Estimating equation `i', `var'/`vtilde':"
		di as res "est cmd: `estring'"
		di as res "(full sample, for debugging; no crossfit)"
		// estimate
		`est_main' if `mname'_sample, `est_options'
		foreach fid of varlist `fidlist' {
			di
			di as res "By resample (foldvar=`fid')
			forvalues k=1/`kfolds' {
				di
				di as res "By fold (`fid'=`k'):
				di as res "est cmd: `est_main' if `mname'_sample & `fid'==`k', `est_options'"
				`est_main' if `mname'_sample & `fid'==`k', `est_options'
			}
		}
	}

	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlist",invtokens(`eqn'.vtlist))
			
			foreach vtilde in `vtlist' {
				mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
				mata: st_local("est_main", return_learner_item(`eqn',"`vtilde'","est_main"))
				mata: st_local("est_options", return_learner_item(`eqn',"`vtilde'","est_options"))
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				di as res "est cmd: `estring'"
				di as res "(full sample, for debugging; no crossfit)"
				// estimate
				`est_main' if `mname'_sample, `est_options'
				foreach fid of varlist `fidlist' {
					di
					di as res "By resample (foldvar=`fid')
					forvalues k=1/`kfolds' {
						di
						di as res "By fold (`fid'=`k'):
						di as res "est cmd: `est_main' if `mname'_sample & `fid'==`k', `est_options'"
						`est_main' if `mname'_sample & `fid'==`k', `est_options'
					}
				}
			}
		}
	}


	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlist",invtokens(`eqn'.vtlist))
			foreach vtilde in `vtlist' {
				mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
				mata: st_local("est_main", return_learner_item(`eqn',"`vtilde'","est_main"))
				mata: st_local("est_options", return_learner_item(`eqn',"`vtilde'","est_options"))
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				di as res "est cmd: `estring'"
				di as res "(full sample, for debugging; no crossfit)"
				// estimate
				`est_main' if `mname'_sample, `est_options'
				foreach fid of varlist `fidlist' {
					di
					di as res "By resample (foldvar=`fid')
					forvalues k=1/`kfolds' {
						di
						di as res "By fold (`fid'=`k'):
						di as res "est cmd: `est_main' if `mname'_sample & `fid'==`k', `est_options'"
						`est_main' if `mname'_sample & `fid'==`k', `est_options'
					}
				}
			}
		}
	}


end
