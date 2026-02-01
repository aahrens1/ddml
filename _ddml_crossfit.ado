*! ddml v1.5.0
*! last edited: 1feb2026
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
							CVCbootnum(integer 500)	///
							yeqn					/// crsossfit Y eqn only
							deqn					/// crsossfit D eqn(s) only
							zeqn					/// crsossfit Z eqn(s) only
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

	// reps = total number of crossfit reps; crossfitted = reps done so far (=0 if none)
	mata: st_local("reps", strofreal(`mname'.nreps))
	mata: st_local("crossfitted", strofreal(`mname'.crossfitted))

	if `reps'==`crossfitted' & "`yeqn'`deqn'`zeqn'"=="" {
		// clear all estimation and crossfit (equation) results from model struct if starting over from scratch
		// reps = #prev crossfits and not crossfitting indiv eqns
		mata: clear_model_results(`mname')
		local crossfitted = 0
	}
	else {
		// just clear any preexisting estimation results from the model struct; leave crossfit results
		mata: clear_model_estimation(`mname')
	}
	// collect number of crossfit reps by eqn
	// always just one Y
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("ycrossfitted", strofreal(`eqn'.crossfitted))
	// always at least one D
	local nameD_1 : word 1 of `nameD'
	mata: `eqn' = (`mname'.eqnAA).get("`nameD_1'")
	mata: st_local("dcrossfitted", strofreal(`eqn'.crossfitted))
	// Z eqn exists for late, iv models
	if `numeqnZ' {
		local nameZ_1 : word 1 of `nameZ'
		mata: `eqn' = (`mname'.eqnAA).get("`nameZ_1'")
		mata: st_local("zcrossfitted", strofreal(`eqn'.crossfitted))
	}
	else {
		local zcrossfitted=0
	}
	// first rep may be >1 if appending; last rep is always nreps
	local lastrep	= `reps'
	
	// set flags
	// standard stacking for model estimation (all eqns) possible only if pystacked used for every eqn
	mata: st_local("stdflag", strofreal(`mname'.stdflag))
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
	// default is to crossfit all eqns
	local yeqnflag		= "`yeqn'"~="" | "`yeqn'`deqn'`zeqn'"==""
	local deqnflag		= "`deqn'"~="" | "`yeqn'`deqn'`zeqn'"==""
	local zeqnflag		= "`zeqn'"~="" | "`yeqn'`deqn'`zeqn'"==""
	// will always be y and d eqns, but z eqn for only iv and interactiveiv
	if `numeqnZ'==0		local zeqnflag=0

	mata: st_local("prefixflag",strofreal(`mname'.prefixflag))
	if `prefixflag'		local prefix `mname'_
	local touse `mname'_sample
	
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
	if `yeqnflag' {
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)

		// if all crossfits already complete, start from the beginning
		if `ycrossfitted'==`reps'		local firstrep = 1
		// otherwise start from the first rep to do
		else							local firstrep	= `ycrossfitted'+1
		
		// fold IDs; first clear list, then assemble
		local fidlist
		forvalues m=`firstrep'/`lastrep' {
			local fidlist `fidlist' `mname'_fid_`m'
		}
		`qui' di as text "Y Fold IDs: `fidlist'"
		
		// number of learners must be >1 for short/pool stacking
		mata: st_local("numlnrY", strofreal(`eqn'.nlearners))
	
		// flag for integrated pystacked code
		mata: st_local("pystacked_Y", strofreal(`eqn'.pystackedmulti))
		
		// shortstack and poolstack variable names
		if `ssflag'		mata: `eqn'.shortstack = "`prefix'Y_`nameY'"
		else			mata: `eqn'.shortstack = ""
		if `psflag'		mata: `eqn'.poolstack = "`prefix'Y_`nameY'"
		else			mata: `eqn'.poolstack = ""
	
		// set/clear treatvar macro
		if "`model'"=="interactive"				local treatvar	`nameD'
		else if ("`model'"=="interactiveiv")	local treatvar	`nameZ'
		else 									local treatvar
		
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
			cvcbootnum(`cvcbootnum')				///
			`nostdstack'							///
			`crossfitother'							///
			`noisily'								///
			`options'
		// resinsert into model struct AA with equations
		mata: (`mname'.eqnAA).put("`nameY'",`eqn')
		
		// update
		local dcrossfitted = `lastrep'
	}
	
	************************** D equation **************************
	if `deqnflag' {

		// set/clear treatvar macro
		if ("`model'"=="interactiveiv") 	local treatvar	`nameZ'
		else								local treatvar

		// if all crossfits already complete, start from the beginning
		if `dcrossfitted'==`reps'		local firstrep = 1
		// otherwise start from the first rep to do
		else							local firstrep	= `dcrossfitted'+1
		
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			
			// fold IDs; first clear list, then assemble
			local fidlist
			forvalues m=`firstrep'/`lastrep' {
				local fidlist `fidlist' `mname'_fid_`m'
			}
			`qui' di as text "D Fold IDs: `fidlist'"
		
			// number of learners must be >1 for short/pool stacking
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
				cvcbootnum(`cvcbootnum')				///
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
			// resinsert into model struct AA with equations
			mata: (`mname'.eqnAA).put("`var'",`eqn')
		}
		
		// update
		local dcrossfitted = `lastrep'
	}

	************************** Z equation **************************
	// iv and interactiveiv only
	if `zeqnflag' {

		// if all crossfits already complete, start from the beginning
		if `zcrossfitted'==`reps'		local firstrep = 1
		// otherwise start from the first rep to do
		else							local firstrep	= `zcrossfitted'+1
		
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			
			// fold IDs; first clear list, then assemble
			local fidlist
			forvalues m=`firstrep'/`lastrep' {
				local fidlist `fidlist' `mname'_fid_`m'
			}
			`qui' di as text "Z Fold IDs: `fidlist'"
		
			// number of learners must be >1 for short/pool stacking
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
				cvcbootnum(`cvcbootnum')				///
				`nostdstack'							///
				`crossfitother'							///
				`noisily'								///
				`options'
			// resinsert into model struct AA with equations
			mata: (`mname'.eqnAA).put("`var'",`eqn')
		}
		
		// update
		local zcrossfitted = `lastrep'
	}
	
	// create sample dfn variable by resample
	create_sample_indicators, mname(`mname')
	
	// set number of fully completed crossfit reps 
	// reps=resamplings to be done, crossfitted=resamplings done so far
	// if appending, crossfitted updated to complete set of reps
	if `numeqnZ'	mata: `mname'.crossfitted = min((`ycrossfitted',`dcrossfitted',`zcrossfitted'))
	else			mata: `mname'.crossfitted = min((`ycrossfitted',`dcrossfitted'))

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
