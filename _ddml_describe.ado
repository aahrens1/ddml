* notes:


program define _ddml_describe

	syntax name(name=mname), [NOCMD all byfold]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	local showcmd	= ("`nocmd'"=="")
	local showall	= ("`all'"~="")

	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	// fold IDs
	forvalues m=1/`nreps' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	
	di as res "Model: `model'"
	di as res "ID: `mname'_id"
	di as res "Fold ID(s): `fidlist'"
	di as res "Sample indicator: `mname'_sample" _c
	qui count if `mname'_sample
	di as res " (N=`r(N)')"
	di as res "Number of folds     =" %3.0f `kfolds'
	di as res "Number of resamples =" %3.0f `nreps'
	
	// equations and learners
	di as res "Dependent variable (Y): `nameY'"
	// mata: `eqn' = (*(`mname'.peqnAA)).get(`mname'.nameY)
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	di as res _col(2) "`nameY' learners:" as res _col(25) "`vtlistY'"
	if `numeqnD' {
		di as res "D equations (`numeqnD'): `nameD'"
		foreach var of varlist `nameD' {
			// mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			di as res _col(2) "`var' learners:" _col(25) "`vtlistD'"
		}
	}
	if `numeqnZ' {
		di as res "Z equations (`numeqnZ'): `nameZ'"
		foreach var of varlist `nameZ' {
			// mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			di as res _col(2) "`var' learners:" _col(25) "`vtlistZ'"
		}
	}
	// learners in detail
	di as res "Y learners (detail):"
	desc_learners `mname', vname(`nameY') etype(yeq)
	if `numeqnD' {
		di as res "D learners (detail):"
		foreach var of varlist `nameD' {
			desc_learners `mname', vname(`var') etype(deq)
		}
	}
	if `numeqnZ' {
		di as res "Z learners (detail):"
		foreach var of varlist `nameZ' {
			desc_learners `mname', vname(`var') etype(zeq)
		}
	}


	/*
	if `showall' {
		di
		di as res "Other:"
		di
		di as res "liststruct(.):"
		mata: liststruct(`mname')
		di
		di as res "Equation pointers:"
		mata: `mname'.eqnlist
		di as res "Corresponding tilde names:"
		mata: `mname'.eqnlistNames
	}
	
	if `crossfitted' {
		_ddml_crossfit_report `mname', `byfold'
	}
	*/
	
	// clear this global from Mata
	mata: mata drop `eqn'
	
	
end

prog define desc_learners

	syntax name(name=mname), vname(string) etype(string)	// etype is yeq, deq or zeq (not dheq)

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("nreps",strofreal(`mname'.nreps))
	
	// used below to indicate set of crossfitting results to report
	local pairs		= 0
	local heqn		= 0
	if ("`etype'"=="yeq") & ("`model'"=="interactive" | "`model'"=="late") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="late") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="ivhd") {
		// includes both deq and dheq
		local heqn	= 1
	}
	
	// mata: `eqn' = (*(`mname'.peqnAA)).get("`vname'")
	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	foreach vtilde in `vtlist' {
		di as res _col(2) "Learner:" _col(15) "`vtilde'"
		mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
		if `heqn' {
			di as res _col(15) "est cmd (D): `estring'"
		}
		else {
			di as res _col(15) "est cmd: `estring'"
		}
		if `crossfitted' {
			di as res _col(15) "<full crossfitting results will go here>"
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`vtilde'","MSE","`m'")))
					di _col(15) "MSE, resample `m': " `mse'
				}
			}
			else {
				forvalues m=1/`nreps' {
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`vtilde'","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`vtilde'","MSE1","`m'")))
					di _col(15) "MSE, resample `m' (0/1): " `mse0' "     " `mse1'
				}
			}
		}
		if `heqn' {
			mata: st_local("estring_h", return_learner_item(`eqn',"`vtilde'","estring_h"))
			di as res _col(15) "est cmd (H): `estring_h'"
			if `crossfitted' {
				di as res _col(15) "<full crossfitting results will go here>"
				forvalues m=1/`nreps' {
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`vtilde'","MSE_h","`m'")))
					di _col(15) "MSE, resample `m': " `mse'
				}
			}
		}
	}

	// clear this global from Mata
	mata: mata drop `eqn'
	
end
