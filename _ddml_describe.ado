*! ddml v1.2
*! last edited: 21 jan 2023
*! authors: aa/ms

program define _ddml_describe

	syntax name(name=mname), [LEARNers CROSSfit ESTimates SAMple all *]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	local all		= "`all'"~=""
	local lflag		= "`learners'"~=""	| `all'
	local cflag		= "`crossfit'"~=""	| `all'
	local eflag		= "`estimates'"~=""	| `all'
	local sflag		= "`sample'"~=""	| `all'
	
	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))	// flag for crossfitting results available
	mata: st_local("ncombos",strofreal(`mname'.ncombos))			// flag for estimation results available
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

	// basic info about equations and learners - always displayed
	di
	di as text "Model:" _col(25) as res "`model', crossfit folds k=" `kfolds' ", resamples r=" `nreps'
	mata: st_local("fclustvar",`mname'.fclustvar)
	if "`fclustvar'"~="" {
		di as res _col(25) "Folds respect clustering by `fclustvar'"
	}
	if "`nameY'"~="" {
		di as text "Dependent variable (Y):" _col(25) as res "`nameY'"
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
		local numlnrY : word count `vtlistY'
		di as text _col(2) "`nameY' learners:" as text _col(25) as res "`vtlistY'"
	}
	if `numeqnD' {
		di as text "D equations (`numeqnD'):" _col(25) as res "`nameD'"
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			local numlnrD : word count `vtlistD'
			local comboD `comboD' * `numlnrD'
			di as text _col(2) "`var' learners:" _col(25) as res "`vtlistD'"
		}
	}
	if `numeqnZ' {
		di as text "Z equations (`numeqnZ'):" _col(25) as res "`nameZ'"
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			local numlnrZ : word count `vtlistZ'
			local comboZ `comboZ' * `numlnrZ'
			di as text _col(2) "`var' learners:" _col(25) as res "`vtlistZ'"
		}
	}
	// report number of specifications
	if "`model'"=="interactive" | "`model'"=="late" {
		local comboY `numlnrY' * `numlnrY'
	}
	else {
		local comboY `numlnrY'
	}
	if "`model'"=="late" | "`model'"=="fiv" {
		local comboD `comboD' `comboD'
	}
	di as text "Specifications:" _col(25) as res `comboY' `comboD' `comboZ' " possible specs" _c
	if `nreps' > 1 {
		di as res " * " `nreps' " crossfit splits = " `comboY' `comboD' `comboZ' * `nreps'
	}
	else {
		di
	}
	
	// sample and folds in detail
	if `sflag' {
		di
		di as text "ID:" _col(25) as res "`mname'_id"
		di as text "Full sample indic.:" _col(25) as res "`mname'_sample" _c
		qui count if `mname'_sample
		di as res " (N=`r(N)')"
		if "`fclustvar'"~="" {
			qui tab `fclustvar' if `mname'_sample
			di as text "Cluster variable:" _col(25) as res "`fclustvar' (N_clust=`r(r)')"
		}
		// di as text "Number of resamples =" _col(25) as res %3.0f `nreps'
		// di as text "Number of folds     =" _col(25) as res %3.0f `kfolds'
		di as text "Fold ID:" _col(25) as res _c
		forvalues m=1/`nreps' {
			local fid : word `m' of `fidlist'
			di as res %~12s "`fid'" _c
		}
		di
		di as text "Fold sample indic.:" _col(25) as res _c
		forvalues m=1/`nreps' {
			local rs : word `m' of `rslist'
			di as res %~12s "`rs'" _c
		}
		di
		di as text "Estimation N:" _col(25) as res _c
		forvalues m=1/`nreps' {
			if `ncombos' {
				local rs : word `m' of `rslist'
				qui count if `rs'
				local N "`: di %2.0f r(N)'"
			}
			else {
				local N "(n.a.)"
			}
			di as res %~12s "`N'" _c
		}
		di
	}
	
	// learners in detail
	if `lflag' {
		di
		di as text "Y learners (detail):"
		desc_learners `mname', vname(`nameY') etype(yeq)
		if `numeqnD' {
			di as text "D learners (detail):"
			foreach var of varlist `nameD' {
				desc_learners `mname', vname(`var') etype(deq)
			}
		}
		if `numeqnZ' {
			di as text "Z learners (detail):"
			foreach var of varlist `nameZ' {
				desc_learners `mname', vname(`var') etype(zeq)
			}
		}
	}
		
	// crossfit results in detail
	if `cflag' & `crossfitted' {
		di
		di as text "Crossfit results (detail):"
		desc_learners `mname', vname(`nameY') etype(yeq) results header
		if `numeqnD' {
			foreach var of varlist `nameD' {
				desc_learners `mname', vname(`var') etype(deq) results
			}
		}
	}
	else if `cflag' {
		di
		di as text "No crossfitting results to display."
	}
	
	// estimate results in detail
	if `eflag' & ("`model'"=="interactive" | "`model'"=="late") & `ncombos' {
		di
		_ddml_estimate_ate_late `mname', `options' results
	}
	else if `eflag' & `ncombos' {
		di
		_ddml_estimate_linear `mname', `options' results
	}
	else if `eflag' {
		di
		di as text "No estimation results to display."
	}
	
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
	if ("`etype'"=="deq") & ("`model'"=="fiv") {
		// includes both deq and dheq
		local heqn	= 1
	}
	
	if ~`showcmd' & `showheader' {
		di as text _col(38) "All" _c
		di as text _col(45) "By fold:" _c
		di
		di as text "Cond. exp." _c
		di as text _col(13) "Learner" _c
		di as text _col(26) "rep" _c
		if `pairs' {
			di as text _col(31) "tv" _c
		}
		di as text _col(38) "MSE" _c
		forvalues k=1/`kfolds' {
			di as text "        `k' " _c
		}
		di
	}
	
	
	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for standard stacking
	mata: st_local("nostdstack", strofreal(`eqn'.nostdstack))
	local stdflag = 1-`nostdstack'

	local firstrow = 1
	foreach vtilde in `vtlist' {
		if `showcmd' {
			di as res _col(2) "Learner:" _col(15) "`vtilde'"
			mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
			// remove tabs and extraneous spaces
			local estring = subinstr("`estring'","	"," ",.)
			local estring = strtrim(stritrim("`estring'"))
			if `heqn' {
				di as res _col(15) "est cmd (D): `estring'"
			}
			else {
				di as res _col(15) "est cmd: `estring'"
			}
		}
		else if `crossfitted' & `stdflag' {
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
		else if `heqn' & `crossfitted' & `stdflag' {
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
	if `crossfitted' & ~`showcmd' {
		if `pairs'==0 {
			forvalues m=1/`nreps' {
				tempname mse_folds
				mata: st_local("mse", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE","`m'")))
				mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE_folds","`m'"))
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
					mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE_h","`m'")))
					mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE_h_folds","`m'"))
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
				tempname mse0_folds mse1_folds
				mata: st_local("mse0", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE0","`m'")))
				mata: st_local("mse1", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE1","`m'")))
				mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE0_folds","`m'"))
				mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE1_folds","`m'"))
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
