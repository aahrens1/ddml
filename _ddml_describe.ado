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
	mata: st_local("ncombos",strofreal(`mname'.ncombos))
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
		local rslist `rslist' `mname'_sample_`m'
	}
	
	di as text "Model: `model'"
	di as text "ID: `mname'_id"
	di as text "Sample indicator: `mname'_sample" _c
	qui count if `mname'_sample
	di as text " (N=`r(N)')"
	di as text "Number of resamples =" %3.0f `nreps'
	di as text "Number of folds     =" %3.0f `kfolds'
	di as text "Fold ID:" _col(16) _c
	forvalues m=1/`nreps' {
		local fid : word `m' of `fidlist'
		di as text %~12s "`fid'" _c
	}
	di
	di as text "Sample indic.:" _col(16) _c
	forvalues m=1/`nreps' {
		local rs : word `m' of `rslist'
		di as text %~12s "`rs'" _c
	}
	di
	di as text "Estimation N:" _col(16) _c
	forvalues m=1/`nreps' {
		if `ncombos' {
			local rs : word `m' of `rslist'
			qui count if `rs'
			local N "`: di %2.0f r(N)'"
		}
		else {
			local N "(n.a.)"
		}
		di as text %~12s "`N'" _c
	}
	di
	di
	
	// equations and learners
	di as res "Dependent variable (Y): `nameY'"
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	di as res _col(2) "`nameY' learners:" as res _col(25) "`vtlistY'"
	if `numeqnD' {
		di as res "D equations (`numeqnD'): `nameD'"
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			di as res _col(2) "`var' learners:" _col(25) "`vtlistD'"
		}
	}
	if `numeqnZ' {
		di as res "Z equations (`numeqnZ'): `nameZ'"
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			di as res _col(2) "`var' learners:" _col(25) "`vtlistZ'"
		}
	}
	
	// learners in detail
	di
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
	// results in detail
	di
	di "Crossfit results (detail):"
	desc_learners `mname', vname(`nameY') etype(yeq) results header
	if `numeqnD' {
		foreach var of varlist `nameD' {
			desc_learners `mname', vname(`var') etype(deq) results
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

	syntax name(name=mname), vname(string) etype(string) [ results header ]	// etype is yeq, deq or zeq (not dheq)
	
	local showcmd = "`results'"==""
	local showheader = "`header'"~=""
	local vnabbrev = abbrev("`vname'",10)

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
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
	
	if ~`showcmd' & `showheader' {
		di as res _col(38) "All" _c
		di as res _col(45) "By fold:" _c
		di
		di as res "Cond. exp." _c
		di as res _col(13) "Learner" _c
		di as res _col(26) "rep" _c
		if `pairs' {
			di as res _col(31) "tv" _c
		}
		di as res _col(38) "MSE" _c
		forvalues k=1/`kfolds' {
			di as res "        `k' " _c
		}
		di
	}
	
	
	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("ssvname",invtokens(`eqn'.shortstack))
	local firstrow = 1
	foreach vtilde in `vtlist' {
		if `showcmd' {
			di as res _col(2) "Learner:" _col(15) "`vtilde'"
			mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
			if `heqn' {
				di as res _col(15) "est cmd (D): `estring'"
			}
			else {
				di as res _col(15) "est cmd: `estring'"
			}
		}
		else if `crossfitted' {
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse_folds
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`vtilde'","MSE","`m'")))
					mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`vtilde'","MSE_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					local lrnabbrev = abbrev("`vtilde'",10)
					di as res _col(12) "`lrnabbrev'" _c
					di _col(26) %2.0f `m' _c
					di _col(34) %8.2f `mse' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse_folds',1,`k') _c
					}
					di
				}
			}
			else {
				forvalues m=1/`nreps' {
					tempname mse0_folds mse1_folds
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`vtilde'","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`vtilde'","MSE1","`m'")))
					mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`vtilde'","MSE0_folds","`m'"))
					mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`vtilde'","MSE1_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					forvalues i=0/1 {
						local lrnabbrev = abbrev("`vtilde'",10)
						di as res _col(12) "`lrnabbrev'" _c
						di _col(26) %2.0f `m' _c
						di _col(31) %2.0f `i' _c
						di _col(34) %8.2f `mse`i'' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse`i'_folds',1,`k') _c
						}
						di
					}
				}
			}
		}
		if `heqn' & `showcmd' {
			mata: st_local("estring_h", return_learner_item(`eqn',"`vtilde'","estring_h"))
			di as res _col(15) "est cmd (H): `estring_h'"
		}
		else if `heqn' & `crossfitted' {
			forvalues m=1/`nreps' {
				tempname mse_h_folds
				mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`vtilde'","MSE_h","`m'")))
				mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`vtilde'","MSE_h_folds","`m'"))
				local lrnabbrev = abbrev("`vtilde'_h",10)
				di as res _col(12) "`lrnabbrev'" _c
				di _col(26) %2.0f `m' _c
				di _col(34) %8.2f `mse_h' _c
				forvalues k=1/`kfolds' {
					di "  " %8.2f el(`mse_h_folds',1,`k') _c
				}
				di
			}
		}
	}
	if "`ssvname'"~="" & `crossfitted' & ~`showcmd' {
		if `pairs'==0 {
			forvalues m=1/`nreps' {
				tempname mse_folds
				mata: st_local("mse", strofreal(return_result_item(`eqn',"`ssvname'","MSE","`m'")))
				mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`ssvname'","MSE_folds","`m'"))
				di as res _col(12) "shortstack" _c
				di _col(26) %2.0f `m' _c
				di _col(34) %8.2f `mse' _c
				forvalues k=1/`kfolds' {
					di "  " %8.2f el(`mse_folds',1,`k') _c
				}
				di
			}
			if `heqn' {
				forvalues m=1/`nreps' {
					tempname mse_h_folds
					mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`ssvname'","MSE_h","`m'")))
					mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`ssvname'","MSE_h_folds","`m'"))
					di as res _col(12) "shortstack_h" _c
					di _col(26) %2.0f `m' _c
					di _col(34) %8.2f `mse_h' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse_h_folds',1,`k') _c
					}
					di
				}
			}
		}
		else {
			forvalues m=1/`nreps' {
				tempname mse0_h_folds mse1_h_folds
				mata: st_local("mse0", strofreal(return_result_item(`eqn',"`ssvname'","MSE0","`m'")))
				mata: st_local("mse1", strofreal(return_result_item(`eqn',"`ssvname'","MSE1","`m'")))
				mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`ssvname'","MSE0_folds","`m'"))
				mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`ssvname'","MSE1_folds","`m'"))
				forvalues i=0/1 {
					local lrnabbrev = abbrev("`vtilde'",10)
					di as res _col(12) "shortstack" _c
					di _col(26) %2.0f `m' _c
					di _col(31) %2.0f `i' _c
					di _col(34) %8.2f `mse`i'' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse`i'_folds',1,`k') _c
					}
					di
				}
			}
		}
	}
	// clear this global from Mata
	mata: mata drop `eqn'
	
end
