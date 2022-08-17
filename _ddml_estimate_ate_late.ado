*** ddml estimation: interactive model
* notes:
* add check somewhere that only a single D is allowed.

program _ddml_estimate_ate_late, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ATET 				///
								ROBust				/// has no effect - vcv always robust or cluster-robust
								CLUster(varname)	///
								vce(string)			///
								ALLest				/// show all regression outputs
								NOTable				/// suppress summary table
								FULLtable			/// show full summary table
								clear				/// deletes all tilde-variables (to be implemented)
								spec(string)		/// specification to post/display
								REP(string)			/// resampling iteration to post/display
								replay				/// model has been estimated, just display results
								avplot				///
								trim(real 0.01)		///
								debug				///
								tnumrows(int 10)	/// for debugging use only
								* ]
	
	if "`fulltable'"~="" {
		// display all rows
		local tnumrows	=.
	}
	if "`debug'"==""	local qui qui
	
	marksample touse
	
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	// model needs to be estimated
	local estflag = "`replay'"==""
	// display summary table
	local tableflag = "`notable'"==""
	// display all regression outpus
	local allflag = "`allest'"~=""
	
	** standard errors
	// local vce is the argument to the Stata option vce(.)
	// SEs are always either robust or cluster-robust
	if "`cluster'"~=""	local vce cluster `cluster'
	else				local vce robust
	
	if ~`crossfitted' {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}

	if "`spec'"=="" {
		local spec "mse"
	}

	// allowable forms
	if "`spec'"=="shortstack"	local spec ss
	if "`spec'"=="minmse"		local spec mse
	if "`rep'"=="mean"			local rep mn
	if "`rep'"=="median"		local rep md

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// locals used below
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	mata: st_local("nreps",strofreal(`mname'.nreps))

	// if rep not specified, default is rep=1 when nreps==1; md if nreps>1
	if "`rep'"=="" & `nreps'>1 {
		local rep md
	}
	else if "`rep'"=="" & `nreps'==1 {
		local rep 1
	}
	if real("`rep'")==. {
		// rep is an integer or mn/md
		if "`rep'"~="mn" & "`rep'"~="md" {
			di as err "error - illegal rep(.) option `rep'"
			exit 198
		}
	}
	else {
		if (real("`rep'")<1) | (real("`rep'")~=int(real("`rep'"))) {
			di as err "error - illegal rep(.) option `rep'"
			exit 198
		}
	}
	// check that rep, if integer, isn't larger than nreps
	if real("`rep'")!=. {
		if `rep'>`nreps' {
			di as err "rep() cannot be larger than `nreps'"
			exit 198
		}
	}

	// check whether shortstack is available for all equations
	if `ssflag' {
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("numlnr",strofreal(cols(`eqn'.vtlist)))
		if `numlnr'==1 {
			local ssflag = 0
		}
		foreach var in `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("numlnr",strofreal(cols(`eqn'.vtlist)))
			if `numlnr'==1 {
				local ssflag = 0
			}
		}
		foreach var in `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("numlnr",strofreal(cols(`eqn'.vtlist)))
			if `numlnr'==1 {
				local ssflag = 0
			}
		}
		if `ssflag' == 0 {
			di as err "warning - shortstack not available for all equations; option ignored"
		}
	}
	
	local ateflag=("`model'"=="interactive")

	// should only ever be 1 D or 1 Z eqn
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	if `numeqnD'>1 {
		di as err "error - model `model' supports only a single treatment variable"
		exit 198
	}
	if `ateflag' & (`numeqnZ'>1) {
		di as err "error - model `model' supports only a single instrument"
		exit 198
	}
	
	************* ESTIMATE ************
	
	// estimate or not
	if ~`estflag' {
		// model has already been estimated, so just recover ncombos macro
		mata: st_local("ncombos", strofreal(`mname'.ncombos))
		if `ncombos'==0 {
			di as err "internal ddml error - model has not be estimated, no results to display"
			exit 198
		}
		// recover matrices
		tempname nmat bmat semat
		mata: `nmat' = (`mname'.estAA).get(("nmat","all"))
		mata: `bmat' = (`mname'.estAA).get(("bmat","all"))
		mata: `semat' = (`mname'.estAA).get(("semat","all"))
		// recover min MSE specs
		forvalues m=1/`nreps' {
			mata: st_local("optspec",(`mname'.estAA).get(("optspec","`m'")))
			local optspec`m' = `optspec'
		}
	}
	else {
		// ATE and LATE both use Y0/Y1
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		mata: st_local("Ytilde",invtokens(`eqn'.vtlist))
		// ATE and LATE both have a D eqn
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameD)
		mata: st_local("Dtilde",invtokens(`eqn'.vtlist))
		if `ateflag' {
			// ATE has a single Dtilde
			_ddml_allcombos `Ytilde'- `Ytilde' - `Dtilde',		///
				addprefix("")									///
				`debug'  
			local Dlist `r(colstr3)'
		}
		else {
			// LATE has a D0/D1 and a single Z
			mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameZ)
			mata: st_local("Ztilde",invtokens(`eqn'.vtlist))
			_ddml_allcombos `Ytilde' - `Ytilde' - `Dtilde' - `Dtilde' - `Ztilde' ,	///
				`debug'																///
				addprefix("")
			local d0list `r(colstr3)'
			local d1list `r(colstr4)'
			local Zlist `r(colstr5)' 
		}
		local y0list `r(colstr1)'
		local y1list `r(colstr2)'
		local ncombos = r(ncombos)
		local tokenlen = `ncombos'*2 -1
			
		tempname nmat bmat semat
		mata: `nmat' = J(`ncombos',5,"")
		mata: `bmat' = J(`ncombos'*`nreps',`numeqnD',.)
		mata: `semat' = J(`ncombos'*`nreps',`numeqnD',.)
		
		// simplest if put into a Mata string matrix
		tokenize `y0list' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',1] = strtrim("``idx''")
		}
		tokenize `y1list' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',2] = strtrim("``idx''")
		}
		if `ateflag' {
			// ATE has a single D
			tokenize `Dlist' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',3] = strtrim("``idx''")
			}
		}
		else {
			// LATE has D0/D1 and a single Z
			tokenize `d0list' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',3] = strtrim("``idx''")
			}
			tokenize `d1list' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',4] = strtrim("``idx''")
			}
			tokenize `Zlist' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',5] = strtrim("``idx''")
			}
		}
		
		*** shortstack names
		if `ssflag' {
			// code works for both ATE and LATE
			local Y0ss	`nameY'_ss
			local Y1ss	`nameY'_ss
			local Dss	`nameD'_ss
			local D0ss	`nameD'_ss
			local D1ss	`nameD'_ss
			local Zss	`nameZ'_ss
		}
		
		forvalues m=1/`nreps' {
			
			// reset locals
			local Y0opt
			local Y1opt
			local Dopt
			local D0opt
			local D1opt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Y0opt",return_learner_item(`eqn',"opt0","`m'"))
			mata: st_local("Y1opt",return_learner_item(`eqn',"opt1","`m'"))
			mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
			if `ateflag' {
				mata: st_local("Dopt",return_learner_item(`eqn',"opt","`m'"))
			}
			else {
				mata: st_local("D0opt",return_learner_item(`eqn',"opt0","`m'"))
				mata: st_local("D1opt",return_learner_item(`eqn',"opt1","`m'"))
				mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
				mata: st_local("Zopt",return_learner_item(`eqn',"opt","`m'"))
			}
			
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			forvalues i = 1/`ncombos' {
				mata: st_local("y0",`nmat'[`i',1])
				mata: st_local("y1",`nmat'[`i',2])
				if `ateflag' {
					mata: st_local("d",`nmat'[`i',3])
				}
				else {
					mata: st_local("d0",`nmat'[`i',3])
					mata: st_local("d1",`nmat'[`i',4])
					mata: st_local("z",`nmat'[`i',5])
				}
				// check if opt for this resample
				// code works for both ATE and LATE
				local isopt
				local isY0opt	: list Y0opt == y0
				local isY1opt	: list Y1opt == y1
				local isDopt	: list Dopt == d
				local isD0opt	: list D0opt == d0
				local isD1opt	: list D1opt == d1
				local isZopt	: list Zopt == z
				local title "DDML model, specification `i'`stext'"
				if `isY0opt' & `isY1opt' & `isDopt' & `isD0opt' & `isD1opt' & `isZopt' {
					local optspec`m' = `i'
					local isopt *
					local title Min MSE `title'
					// save in AA
					mata: (`mname'.estAA).put(("optspec","`m'"),"`i'")
				}
				// code works for both ATE and LATE
				`qui' _ddml_ate_late if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')			///
					y0tilde(`y0') y1tilde(`y1')							///
					dtilde(`d') d0tilde(`d0') d1tilde(`d1')				///
					ztilde(`z')											///
					spec(`i') rep(`m')									///
					mname(`mname')										///
					title(`title')										///
					trim(`trim')										///
					vce(`vce') `atet'
				
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
				
			}
			
			if `ssflag' {
				
				// code works for both ATE and LATE
				local title "Shortstack DDML model`stext'"
				`qui' _ddml_ate_late if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')			///
					y0tilde(`Y0ss') y1tilde(`Y1ss')						///
					dtilde(`Dss') d0tilde(`D0ss') d1tilde(`D1ss')		///
					ztilde(`Zss')										///
					spec(ss) rep(`m')									///
					mname(`mname')										///
					title(`title')										///
					trim(`trim')										///
					vce(`vce') `atet'
			
			}
		}

		// we make a copy of the MSE-optimal model for each m
		forvalues m=1/`nreps' {
			tempname Bopt
			mata: st_local("optspec",(`mname'.estAA).get(("optspec","`m'")))
			mata: `Bopt' = (`mname'.estAA).get(("`optspec'","`m'"))
			mata: (`mname'.estAA).put(("mse","`m'"),`Bopt')
			mata: mata drop `Bopt'
		}
		
		// aggregate across resamplings
		if `nreps' > 1 {
 			`qui' _ddml_ate_late, mname(`mname') spec(mse) medmean(mn) title("Mean over min-mse specifications") vce(`vce') `atet' // min-mse specification
 			`qui' _ddml_ate_late, mname(`mname') spec(mse) medmean(md) title("Median over min-mse specifications") vce(`vce') `atet' // min-mse specification
			// numbered specifications
			forvalues i = 1/`ncombos' {
				local title "DDML model, specification `i' (mean)"
				`qui' _ddml_ate_late, mname(`mname') spec(`i') medmean(mn) title(`title') vce(`vce') `atet'
				local title "DDML model, specification `i' (median)"
				`qui' _ddml_ate_late, mname(`mname') spec(`i') medmean(md) title(`title') vce(`vce') `atet'
			}
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean)"
				`qui' _ddml_ate_late, mname(`mname') spec(ss) medmean(mn) title(`title') vce(`vce') `atet'
				local title "Shortstack DDML model (median)"
				`qui' _ddml_ate_late, mname(`mname') spec(ss) medmean(md) title(`title') vce(`vce') `atet'
			}
		}
		
		// estimation complete
		mata: `mname'.ncombos = `ncombos'
		mata: (`mname'.estAA).put(("nmat","all"),`nmat')
		mata: (`mname'.estAA).put(("bmat","all"),`bmat')
		mata: (`mname'.estAA).put(("semat","all"),`semat')
		
	}
	
	************** REPORT RESULTS **************

	if `allflag' {
		forvalues m=1/`nreps' {
			// all combos including min MSE model
			forvalues i=1/`ncombos' {
				di
				_ddml_ate_late, mname(`mname') spec(`i') rep(`m') replay
				di
			}
		}
	}
	if `ssflag' & `allflag' {
		forvalues m=1/`nreps' {
			di
			_ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
			di
		}
	}
	if (`nreps' > 1) & `allflag' {
		// numbered specifications
		forvalues i = 1/`ncombos' {
			di
			_ddml_ate_late, mname(`mname') spec(`i') rep(mn) replay
			di
			_ddml_ate_late, mname(`mname') spec(`i') rep(md) replay
			di
		}
		// shortstack
		if `ssflag' {
			_ddml_ate_late, mname(`mname') spec(ss) rep(mn) replay
			di
			_ddml_ate_late, mname(`mname') spec(ss) rep(md) replay
			di
		}
	}
	
	*** Results ***
	// counter for number of rows in summary table
	local rowcount 0
	// optional table of all results
	if `tableflag' {
		di
		di as text "DDML estimation results:"
		di as text "spec  r" %14s "Y0 learner" _c
		di as text           %14s "Y1 learner" _c
		if `ateflag' {
			di as text           %14s "D learner" _c
		}
		else {
			di as text           %14s "D0 learner" _c
			di as text           %14s "D1 learner" _c
		}
		di as text %10s "b" %10s "SE" _c
		if ~`ateflag' {
			di as text           %14s "Z learner" _c
		}
		di
		forvalues m=1/`nreps' {
			forvalues i=1/`ncombos' {
				local ++rowcount
				if `rowcount' <= `tnumrows' {
					mata: st_local("yt0",abbrev(`nmat'[`i',1],13))
					mata: st_local("yt1",abbrev(`nmat'[`i',2],13))
					if "`optspec`m''"=="`i'" {
						di "*" _c
					}
					else {
						di " " _c
					}
					local specrep `: di %3.0f `i' %3.0f `m''
					// pad out to 6 spaces
					local specrep = (6-length("`specrep'"))*" " + "`specrep'"
					local rcmd stata ddml estimate `mname', spec(`i') rep(`m') replay notable
					di %6s "{`rcmd':`specrep'}" _c
					di as res %14s "`yt0'" _c
					di as res %14s "`yt1'" _c
					if `ateflag' {
						mata: st_local("dt",`nmat'[`i',3])
						di as res %14s "`dt'" _c
					}
					else {
						mata: st_local("dt0",abbrev(`nmat'[`i',3],13))
						mata: st_local("dt1",abbrev(`nmat'[`i',4],13))
						di as res %14s "`dt0'" _c
						di as res %14s "`dt1'" _c
					}
					mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
					mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
					di as res %10.3f `b' _c
					local pse (`: di %6.3f `se'')
					di as res %10s "`pse'" _c
					if ~`ateflag' {
						mata: st_local("zt",abbrev(`nmat'[`i',5],13))
						di as res %14s "`zt'" _c
					}
					di

					if `ssflag' & (`rowcount' <= `tnumrows') {
						local ++rowcount
						qui _ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
						tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
						mat `btemp' = e(b)
						mat `Vtemp' = e(V)
						local specrep `: di "ss" %3.0f `m''
						// pad out to 6 spaces
						local specrep = "  " + "`specrep'"
						local rcmd stata ddml estimate `mname', spec(ss) rep(`m') replay notable
						di %6s "{`rcmd':`specrep'}" _c
						di as res %14s "[shortstack]" _c
						di as res %14s "[ss]" _c
						di as res %14s "[ss]" _c
						if ~`ateflag' {
							di as res %14s "[ss]" _c
						}
						di as res %10.3f el(`btemp',1,1) _c
						local pse (`: di %6.3f sqrt(el(`Vtemp',1,1))')
						di as res %10s "`pse'" _c
						if ~`ateflag' {
							di as res %14s "[ss]" _c
						}
						di
					}
				}
			}
		}
		if `rowcount' > `tnumrows' {
			local rcmd stata ddml estimate `mname', replay fulltable
			di %6s "{`rcmd':   ...  }" _c
			di as text "<-click or type " as res "ddml estimate, replay full" as text " to display full summary"
		}
		di as res "*" as text " = minimum MSE specification for that resample."
	}
		
	if `nreps' > 1 {
		di
		di as text "Mean/med.  Y0 learner" _c
		di as text           %14s "Y1 learner" _c
		if `ateflag' {
			di as text           %14s "D learner" _c
		}
		else {
			di as text           %14s "D0 learner" _c
			di as text           %14s "D1 learner" _c
		}
		di as text %10s "b" %10s "SE" _c
		if ~`ateflag' {
			di as text           %14s "Z learner" _c
		}
		di
		foreach medmean in mn md {
			qui _ddml_ate_late, mname(`mname') spec(mse) rep(`medmean') replay
			tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
			mat `btemp' = e(b)
			mat `Vtemp' = e(V)
			local specrep `: di "mse" %3s "`medmean'"' //'
			// pad out to 6 spaces
			local specrep = " " + "`specrep'"
			local rcmd stata ddml estimate `mname', spec(mse) rep(`medmean') replay notable
			di %6s "{`rcmd':`specrep'}" _c
			di as res %14s "[min-mse]" _c
			di as res %14s "[mse]" _c
			di as res %14s "[mse]" _c
			if ~`ateflag' {
				di as res %14s "[mse]" _c
			}
			di as res %10.3f el(`btemp',1,1) _c
			local pse (`: di %6.3f sqrt(el(`Vtemp',1,1))')
			di as res %10s "`pse'" _c
			if ~`ateflag' {
				di as res %14s "[mse]" _c
			}
			di
			if `ssflag' {
				qui _ddml_ate_late, mname(`mname') spec(ss) rep(`medmean') replay
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep `: di "ss" %3s "`medmean'"' //'
				// pad out to 6 spaces
				local specrep = "  " + "`specrep'"
				local rcmd stata ddml estimate `mname', spec(ss) rep(`medmean') replay notable
				di as res %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				di as res %14s "[ss]" _c
				di as res %14s "[ss]" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				di as res %10.3f el(`btemp',1,1) _c
				local pse (`: di %6.3f sqrt(el(`Vtemp',1,1))')
				di as res %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				di
			}
		}
	}

	// post selected estimates; rep is the resample number (default=1)
	di
	_ddml_ate_late, mname(`mname') spec(`spec') rep(`rep') replay
	di
	
	if `nreps' > 1 & ("`rep'"=="mn" | "`rep'"=="md") {
		tempvar bhat
		svmat e(b_resamples), names(`bhat')
		// variables in Stata will look like _000000A1, _000000A2, etc. and will disappear as temps after exit
		local dnames : colnames e(b)
		
		di as text "Summary over " `nreps' " resamples:"
		di as text %12s "D eqn" %10s "mean" %10s "min" %10s "p25" %10s "p50" %10s "p75" %10s "max"
		local i 1
		foreach vn in `dnames' {
			di %12s as text "`vn'" _col(15) _c
			qui sum `bhat'`i', detail
			di as res %10.4f r(mean) _c
			di as res %10.4f r(min) _c
			di as res %10.4f r(p25) _c
			di as res %10.4f r(p50) _c
			di as res %10.4f r(p75) _c
			di as res %10.4f r(max)
			local ++i
		}
	}
	
	// temp Mata objects no longer needed
	foreach obj in `eqn' `nmat' `bmat' `semat' {
		cap mata: mata drop `obj'
	}


end

// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	syntax [anything] , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end

