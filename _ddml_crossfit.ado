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
	di as text "All learners: `vtlist'"
	
	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname', fidlist(`fidlist')
	}
	else {
		// Y equation, all learners
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		if "`model'"=="interactive" {
			local treatvar	`nameD'
		}
		else if ("`model'"=="late") {
			local treatvar	`nameZ'
		}
		else {
			// clear local
			local treatvar
		}
		mata: st_local("ssname",`eqn'.shortstack)
		di as text "Cross-fitting Y equation: `nameY'"
		crossfit if `touse',						///
			ename(`eqn') noreplace					///
			foldvar(`fidlist')						///
			shortstack(`ssname')					///
			treatvar(`treatvar')					///
			resid `noisily'
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
				mata: st_local("ssname",`eqn'.shortstack)
				di as text "Cross-fitting D equation: `var'"
				// All learners for each D eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					treatvar(`treatvar')					///
					resid `noisily'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
		// Z equations
		if `numeqnZ' {
			foreach var of varlist `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("ssname",`eqn'.shortstack)
				di as text "Cross-fitting Z equation: `var'"
				// All learners for each Z eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					resid `noisily'
				mata: (`mname'.eqnAA).put("`var'",`eqn')
			}
		}
	
		// set flag on model struct
		mata: `mname'.crossfitted = 1
	
		// report results by equation type with resamplings grouped together
		di
		di as res "Reporting crossfitting results:"
		_ddml_describe `mname'
	}
	
	// drop temp eqn struct	
	mata: mata drop `eqn'
	
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
