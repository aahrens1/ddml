*! ddml v1.2
*! last edited: 20 july 2023
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
							*						///
							]

	local debugflag		= "`debug'"~=""
	if "`noisily'"==""	local qui qui
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
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
	
	// set flags
	local ssflag	= "`shortstack'"~=""
	local psflag	= "`poolstack'"~=""
	local stdflag	= "`nostdstack'"==""
	
	*** extract details of estimation
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens((`mname'.nameD)))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("reps",strofreal(`mname'.nreps))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
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
	
	// equations and learners
	
	// will be set to zero if any eqn doesn't use pystacked, or if any
	// or if any equation has pystacked with a single learner
	local allpystackedmulti=1
	
	// will always be a Y eqn
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	`qui' di as text "Y eqn learners: `vtlistY'"
	local vtlist `vtlistY'
	// used to track minimum number of learners in an equation; must be >1 for short/pool stacking
	mata: st_local("minlearners", strofreal(`eqn'.nlearners))
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	local allpsm		= 1
	if `pystackedmulti'	local minlearners=`pystackedmulti'
	else				local allpsm=0
	
	// if pystackedmulti=0, not pystacked so every eqn goes to general crossfit code
	// if pystackedmulti=1, pystacked with 1 learner so every eqn goes to general crossfit code
	if `pystackedmulti'<=1	local allpystackedmulti=0

	// will always be a D eqn
	`qui' di as text "D equations (`numeqnD'): `nameD'"
	foreach var of varlist `nameD' {
		`qui' di as text _col(5) "D equation `var':"
		mata: `eqn' = (`mname'.eqnAA).get("`var'")
		mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
		`qui' di as text _col(10) "learners: `vtlistD'"
		local vtlist `vtlist' `vtlistD'
		// update minlearners
		mata: st_local("numlnrD", strofreal(`eqn'.nlearners))
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti'			local numlnrD=`pystackedmulti'
		else						local allpsm=0
		if `numlnrD'<`minlearners'	local minlearners=`numlnrD'
		if `pystackedmulti'<=1		local allpystackedmulti=0
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
			// update minlearners
			mata: st_local("numlnrZ", strofreal(`eqn'.nlearners))
			mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
			if `pystackedmulti'			local numlnrZ=`pystackedmulti'
			else						local allpsm=0
			if `numlnrZ'<`minlearners'	local minlearners `numlnrZ'
			if `pystackedmulti'<=1		local allpystackedmulti=0
		}
	}

	// check minlearners
	// short-stacking and pooled-stacking require multiple learners in all equations
	if (`ssflag' | `psflag') & `minlearners'==1 {
		di as text "`shortstack' `poolstack' requested but must have multiple learners in all equations; option ignored"
		mata: `mname'.ssflag = 0
		local ssflag=0
		local shortstack
		mata: `mname'.psflag = 0
		local psflag=0
		local poolstack
	}
	// pooled-stacking requires pystacked multilearner in all equations
	if `psflag' & (~`allpsm' | ~`stdflag') {
		di as text "poolstack requires a single use of pystacked with multiple learners in all equations; option ignored"
		mata: `mname'.psflag = 0
		local psflag=0
		local poolstack
	}
	// fiv/lie model special case - pystacked multilearner does not support short-stacking
	if `ssflag' & "`model'"=="fiv" & `allpystackedmulti' {
		di as text "shortstack requested but fiv model does not support pystacked integration"
		di as text "to short-stack with fiv model, must have multiple learners in all equations; option ignored"
		mata: `mname'.ssflag = 0
		local ssflag=0
		local shortstack
	}
	// fiv/lie model does not support pooled-stacking
	if `psflag' & "`model'"=="fiv" {
		di as text "poolstack requested but not supported for fiv model; option ignored"
		mata: `mname'.psflag = 0
		local psflag=0
		local poolstack
	}
	// update allpystackedmulti
	if `allpystackedmulti' & "`crossfitother'"=="" {
		mata: `mname'.allpystackedmulti = 1
		`qui' di as text "crossfitting all eqns will use pystacked-specific code"
	}
	else {
		local allpystackedmulti = 0
		mata: `mname'.allpystackedmulti = 0
		// general (non-pystacked) code supports short-stacking but not standard or pooled stacking
		local stdflag = 0
		local psflag = 0
		`qui' di as text "crossfitting all eqns will use general (not pystacked-specific) code"
	}

	// update model struct flags for stacking
	if `stdflag'	mata: `mname'.stdflag	= 1
	else 			mata: `mname'.stdflag	= 0
	if `ssflag'		mata: `mname'.ssflag	= 1
	else 			mata: `mname'.ssflag	= 0
	if `psflag' 	mata: `mname'.psflag	= 1
	else			mata: `mname'.psflag	= 0

	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname', fidlist(`fidlist')
	}
	else {
		// Y equation, all learners
		// will always be a Y eqn
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		
		// shortstack and poolstack variable names
		if `ssflag'		mata: `eqn'.shortstack = "Y_`nameY'"
		else			mata: `eqn'.shortstack = ""
		if `psflag'		mata: `eqn'.poolstack = "Y_`nameY'"
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
		
		if ("`model'"=="partial") di as text "Cross-fitting E[y|X] equation: `nameY'"
		if ("`model'"=="fiv"|"`model'"=="interactiveiv") di as text "Cross-fitting E[y|X,Z] equation: `nameY'"
		if ("`model'"=="interactive") di as text "Cross-fitting E[y|X,D] equation: `nameY'"
		if ("`model'"=="iv") di as text "Cross-fitting E[y|X] equation: `nameY'"

		crossfit if `touse',						///
			ename(`eqn') noreplace					///
			pystackedmulti(`allpystackedmulti')		///
			foldvar(`fidlist')						///
			firstrep(`firstrep')					///
			treatvar(`treatvar')					///
			model(`model')							///
			`options' `nostdstack' `noisily'
		// resinsert into model struct AA with equations
		mata: (`mname'.eqnAA).put("`nameY'",`eqn')

		// D equations
		if `numeqnD' {
			if ("`model'"=="interactiveiv") {
				local treatvar	`nameZ'
				local allowallzero allowallzero // for the case where D is always zero when Z=0
			}
			else {
				// clear local
				local treatvar
			}
			foreach var of varlist `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("numlnrD",strofreal(cols(`eqn'.vtlist)))
				// shortstack and poolstack variable names
				if `ssflag'		mata: `eqn'.shortstack = "D_`var'"				
				else			mata: `eqn'.shortstack = ""
				if `psflag'		mata: `eqn'.poolstack = "D_`var'"				
				else			mata: `eqn'.poolstack = ""
				if ("`model'"=="partial") di as text "Cross-fitting E[D|X] equation: `var'"
				if ("`model'"=="interactiveiv") di as text "Cross-fitting E[D|X,Z] equation: `var'"
				if ("`model'"=="fiv") di as text "Cross-fitting E[D|X,Z] and E[D|X] equation: `var'"
				if ("`model'"=="interactive"|"`model'"=="iv") di as text "Cross-fitting E[D|X] equation: `var'"

				// All learners for each D eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					pystackedmulti(`allpystackedmulti')		///
					foldvar(`fidlist')						///
					firstrep(`firstrep')					///
					treatvar(`treatvar')					///
					model(`model')							///
					`options' `nostdstack' `noisily'		///
					`allowallzero'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
		// Z equations
		if `numeqnZ' {
			foreach var of varlist `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("numlnrZ",strofreal(cols(`eqn'.vtlist)))
				// shortstack and poolstack variable names
				if `ssflag'		mata: `eqn'.shortstack = "Z_`var'"				
				else			mata: `eqn'.shortstack = ""
				if `psflag'		mata: `eqn'.poolstack = "Z_`var'"				
				else			mata: `eqn'.poolstack = ""
				di as text "Cross-fitting E[Z|X]: `var'"
				// All learners for each Z eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					pystackedmulti(`allpystackedmulti')		///
					foldvar(`fidlist')						///
					firstrep(`firstrep')					///
					model(`model')							///
					`options' `nostdstack' `noisily'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
		
		// create sample dfn variable by resample
		create_sample_indicators, mname(`mname') stdflag(`stdflag') ssflag(`ssflag') psflag(`psflag')
		
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
	}
	
	// drop temp eqn struct	
	mata: mata drop `eqn'
	
end

// creates sample indicators and leaves them in memory
// same naming as main sample variable but with a _m resample subscript
program create_sample_indicators
	version 16
	syntax [anything], mname(name) [ stdflag(integer 1) ssflag(integer 1) psflag(integer 1) ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	mata: st_local("model",`mname'.model)
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",invtokens(`mname'.nameY))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'
	
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("shortstack", `eqn'.shortstack)
	mata: st_local("poolstack", `eqn'.poolstack)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	if `ssflag'		local vtlistY `vtlistY' `shortstack'_ss
	if `psflag'		local vtlistY `vtlistY' `poolstack'_ps
	
	local lieflag = 0
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: st_local("poolstack", `eqn'.poolstack)
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			if `ssflag'		local vtlistD `vtlistD' `shortstack'_ss
			if `psflag'		local vtlistD `vtlistD' `poolstack'_ps
			local Dt_list `Dt_list' `vtlistD'
		}
	}
	
	if `numeqnD' & `lieflag' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: st_local("poolstack", `eqn'.poolstack)
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			foreach vn in `vtlistD' {
				local vtlistD_h `vtlistD_h' `vn'_h
			}
			if `ssflag'	local vtlistD_h `vtlistD_h' `shortstack'_h_ss
			if `psflag'	local vtlistD_h `vtlistD_h' `poolstack'_h_ps
			local DHt_list `DHt_list' `vtlistD_h'
		}
	}
	
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: st_local("poolstack", `eqn'.poolstack)
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			if `ssflag'		local vtlistZ `vtlistZ' `shortstack'_ss
			if `psflag'		local vtlistZ `vtlistZ' `poolstack'_ps
			local Zt_list `Zt_list' `vtlistZ'
		}
	}

	forvalues m=1/`nreps' {
		
		// reset local
		local vtlist
		
		cap drop `mname'_sample_`m' 
		qui gen byte `mname'_sample_`m' = `mname'_sample
		label var `mname'_sample_`m' "Sample indicator for rep `m'"
				
		// Y
		if "`model'"=="interactive" | "`model'"=="interactiveiv" {
			foreach vt in `vtlistY' {
				local vtlist `vtlist' `vt'0
				local vtlist `vtlist' `vt'1
			}
		}
		else {
			local vtlist `vtlistY'
		}
		// D and Z
		if "`model'"=="interactiveiv" {
			foreach vt in `vtlistD' {
				local vtlist `vtlist' `vt'0
				local vtlist `vtlist' `vt'1
			}
			local vtlist `vtlist' `Zt_list'
		}
		else {
			local vtlist `vtlist' `Dt_list' `DHt_list' `Zt_list'
		}
		
		foreach vt in `vtlist' {
			qui replace `mname'_sample_`m' = 0 if `vt'_`m'==.
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
	
	// no longer needed
	mata: mata drop `eqn'
	
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
