*** ddml cross-fitting

* notes:
* why is vtype field needed?
* crossfitting code in separate program
* reporting code in separate subroutine below
* debug reporting in separate subroutine below
* number of resamples set here in reps(.) option
* in eqn struct, (*p).idVtilde is a problem ... but is it needed? currently commented out.
* add noisily option

program _ddml_crossfit, eclass sortpreserve

	syntax [anything] ,								/// 
							[						///
							NOIsily					///
							debug					/// 
							Robust					///
							mname(name)				///
							shortstack				///
							Verbose 				///
							]

	// no checks included yet
	
	local debugflag		= "`debug'"~=""
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// clear any preexisting equation results from the model struct
	mata: clear_model_results(`mname')
	
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
	
	// set shortstack indictor on model struct
	if "`shortstack'" ~= "" {
		if ~(`reps' > 1) {
			di as err "warning - shortstack option ignored - reps must be > 1"
			// clear macro
			local shortstack
		}
		else {
			mata: `mname'.ssflag = 1
		}
	}
	
	di as text "Model: `model'"
	// fold IDs
	forvalues m=1/`reps' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	di as text "Fold IDs: `fidlist'"
	
	// equations and learners
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	di as text "Y eqn learners (`numlnrY'): `vtlistY'"
	local vtlist `vtlistY'
	if `numeqnD' {
		di as text "D equations (`numeqnD'): `nameD'"
		foreach var of varlist `nameD' {
			di as text _col(5) "D equation `var':"
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			di as text _col(10) "learners: `vtlistD'"
			local vtlist `vtlist' `vtlistD'
		}
	}
	if `numeqnZ' {
		di as text "Z equations (`numeqnZ'): `nameZ'"
		foreach var of varlist `nameZ' {
			di as text _col(5) "Z equation `var':"
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			di as text _col(10) "learners: `vtlistZ'"
			local vtlist `vtlist' `vtlistZ'
		}
	}
	//di as text "All learners: `vtlist'" // AA: don't think this output is required
	
	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname', fidlist(`fidlist')
	}
	else {
		// Y equation, all learners
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		if "`model'"=="interactive" {
			local treatvar	`nameD'
			local resid
			local residy
		}
		else if ("`model'"=="late") {
			local treatvar	`nameZ'
			local resid
			local residy 
		}
		else if ("`model'"=="ivhd") {
			local treatvar
			local resid
			local residy resid
		} 
		else {
			local treatvar
			local resid resid
			local residy resid
		}
		if "`shortstack'"~="" {
			local ssname `nameY'_ss
			mata: `eqn'.shortstack = "`ssname'"
		}
		else {
			// clear local
			local ssname		
		}
		if ("`model'"=="partial") di as text "Cross-fitting E[Y|X] equation: `nameY'"
		if ("`model'"=="ivhd"|"`model'"=="late") di as text "Cross-fitting E[Y|X,Z] equation: `nameY'"
		if ("`model'"=="interactive") di as text "Cross-fitting E[Y|X,D] equation: `nameY'"
		if ("`model'"=="iv") di as text "Cross-fitting E[Y|X] equation: `nameY'"
		crossfit if `touse',						///
			ename(`eqn') noreplace					///
			foldvar(`fidlist')						///
			shortstack(`ssname')					///
			treatvar(`treatvar')					///
			`residy' `noisily'
		// resinsert into model struct AA with equations
		mata: (`mname'.eqnAA).put("`nameY'",`eqn')

		// D equations
		if `numeqnD' {
			if ("`model'"=="late") {
				local treatvar	`nameZ'
			}
			else {
				// clear local
				local treatvar
			}
			foreach var of varlist `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				if "`shortstack'"~="" {
					local ssname `var'_ss
					mata: `eqn'.shortstack = "`ssname'"
				}
				else {
					// clear local
					local ssname		
				}
				if ("`model'"=="partial") di as text "Cross-fitting E[D|X] equation: `var'"
				if ("`model'"=="late") di as text "Cross-fitting E[D|X,Z] equation: `var'"
				if ("`model'"=="ivhd") di as text "Cross-fitting E[D|X,Z] and E[D|X] equation: `var'"
				if ("`model'"=="interactive"|"`model'"=="iv") di as text "Cross-fitting E[D|X] equation: `var'"
				// All learners for each D eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					treatvar(`treatvar')					///
					`resid' `noisily'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
		// Z equations
		if `numeqnZ' {
			foreach var of varlist `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				if "`shortstack'"~="" {
					local ssname `var'_ss
					mata: `eqn'.shortstack = "`ssname'"
				}
				else {
					// clear local
					local ssname		
				}
				di as text "Cross-fitting E[Z|X]: `var'"
				// All learners for each Z eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					`resid' `noisily'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
		
		// create sample dfn variable by resample
		create_sample_indicators, mname(`mname')
		
		// set flag on model struct
		mata: `mname'.crossfitted = 1
	
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

	syntax [anything], mname(name)
	
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
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			local Dt_list `Dt_list' `vtlistD'
		}
	}
	
	if `numeqnD' & `lieflag' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			foreach vn in `vtlistD' {
				local vtlistD_h `vtlistD_h' `vn'_h
			}
			local DHt_list `DHt_list' `vtlistD_h'
		}
	}
	
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			local Zt_list `Zt_list' `vtlistZ'
		}
	}

	forvalues m=1/`nreps' {
		
		// reset local
		local vtlist
		
		cap drop `mname'_sample_`m' 
		qui gen byte `mname'_sample_`m' = `mname'_sample
				
		// Y
		if "`model'"=="interactive" | "`model'"=="late" {
			foreach vt in `vtlistY' {
				local vtlist `vtlist' `vt'0
				local vtlist `vtlist' `vt'1
			}
		}
		else {
			local vtlist `vtlistY'
		}
		// D and Z
		if "`model'"=="late" {
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
		qui gen `mname'_sample_md = `mname'_sample_mn
	}
	
	
end

program report_debugging

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
