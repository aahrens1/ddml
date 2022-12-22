*! ddml v1.0
*! last edited: 13 dec 2022
*! authors: aa/ms

program _ddml_estimate_ate_late, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ATET 				///
								ATEU				///
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
	// teffect is either ATE, ATET, ATEU or LATE
	local check_opt : word count `atet' `ateu'
	if `check_opt' > 1 {
		di as err "error - may request only one of atet or ateu"
		exit 198
	}
	else if `check_opt'==1 & ~`ateflag' {
		di as err "error - `atet'`ateu' not available with LATE (interactiveiv)"
		exit 198
	}
	if ~`ateflag' {
		// it's LATE
		local teffect LATE
	}
	else if "`atet'`ateu'"=="" {
		// default
		local teffect ATE
	}
	else if "`atet'"~="" {
		local teffect ATET
	}
	else {
		local teffect ATEU
	}
	
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
				local title "DDML model, specification `i'`stext' (`teffect')"
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
					vce(`vce')											///
					te(`teffect')
				
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
				
			}
			
			if `ssflag' {
				
				// code works for both ATE and LATE
				local title "Shortstack DDML model`stext' (`teffect')"
				`qui' _ddml_ate_late if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')			///
					y0tilde(`Y0ss') y1tilde(`Y1ss')						///
					dtilde(`Dss') d0tilde(`D0ss') d1tilde(`D1ss')		///
					ztilde(`Zss')										///
					spec(ss) rep(`m')									///
					mname(`mname')										///
					title(`title')										///
					trim(`trim')										///
					vce(`vce')											///
					te(`teffect')
			
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
 			`qui' _ddml_ate_late, mname(`mname') spec(mse) medmean(mn) title("Mean over min-mse specifications (`teffect')") vce(`vce') te(`teffect') // min-mse specification
 			`qui' _ddml_ate_late, mname(`mname') spec(mse) medmean(md) title("Median over min-mse specifications (`teffect')") vce(`vce') te(`teffect') // min-mse specification
			// numbered specifications
			forvalues i = 1/`ncombos' {
				local title "DDML model, specification `i' (mean) (`teffect')"
				`qui' _ddml_ate_late, mname(`mname') spec(`i') medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "DDML model, specification `i' (median) (`teffect')"
				`qui' _ddml_ate_late, mname(`mname') spec(`i') medmean(md) title(`title') vce(`vce') te(`teffect')
			}
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean) (`teffect')"
				`qui' _ddml_ate_late, mname(`mname') spec(ss) medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Shortstack DDML model (median) (`teffect')"
				`qui' _ddml_ate_late, mname(`mname') spec(ss) medmean(md) title(`title') vce(`vce') te(`teffect')
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
		di as text "DDML estimation results (`teffect'):"
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
					local rcmd stata ddml estimate, mname(`mname') spec(`i') rep(`m') replay notable
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
				}
			}

			if `ssflag' & (`rowcount' <= `tnumrows') {
				local ++rowcount
				qui _ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep `: di "ss" %3.0f `m''
				// pad out to 6 spaces
				local specrep = "  " + "`specrep'"
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`m') replay notable
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
		if `rowcount' > `tnumrows' {
			local rcmd stata ddml estimate, mname(`mname') replay fulltable
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
			local rcmd stata ddml estimate, mname(`mname') spec(mse) rep(`medmean') replay notable
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
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`medmean') replay notable
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


* code below currently supports only a single treatment variable, but coded for multiple variables in places
program _ddml_ate_late, eclass

	version 13

	syntax [anything] [if] [in] ,				///
						[						///
							yvar(varname)		///
							dvar(varname)		///
							zvar(varname)		///
							y0tilde(name)		///
							y1tilde(name)		///
							dtilde(name)		///
							d0tilde(name)		///
							d1tilde(name)		///
							ztilde(name)		///
							mname(name)			///
							spec(string)		///
							rep(string)			///
							replay				///
							title(string)		///
							medmean(string)		///
							vce(string) 		///
							TEffect(name)		/// either ATE, ATET or ATEU
							trim(real 0)	    /// value should be provided by calling program
						]
		
	mata: st_local("model",`mname'.model)
	
	marksample touse
	
	if "`replay'"=="" & "`medmean'"=="" {	// estimate from scratch
		
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		// 0/1 etc
		local y0		`y0tilde'
		local y1		`y1tilde'
		local d			`dtilde'
		local d0		`d0tilde'
		local d1		`d1tilde'
		local z			`ztilde'
		// add suffixes
		local y0_m		`y0tilde'0_`rep'
		local y1_m		`y1tilde'1_`rep'
		local d_m		`dtilde'_`rep'
		local d0_m		`d0tilde'0_`rep'
		local d1_m		`d1tilde'1_`rep'
		local z_m		`ztilde'_`rep'

		// estimation samples may differ across conditional expectations
		if "`model'"=="interactive" {
			markout `touse' `y0_m' `y1_m' `d_m'
		}
		else {
			markout `touse' `y0_m' `y1_m' `d0_m' `d1_m' `z_m'
		}
		
		local vce1: word 1 of `vce'
		if "`vce1'"=="cluster" {
			local clustvar : word 2 of `vce'
		}
		
		if "`model'"=="interactive" {
			mata: ATE("`teffect'","`yvar'","`dvar'","`y0_m'", "`y1_m'", "`d_m'","`touse'","`b'","`V'","`clustvar'","`mname'_fid_`rep'",`trim')
		}
		else {
			mata: LATE("`yvar'","`dvar'","`zvar'","`y0_m'", "`y1_m'", "`d0_m'","`d1_m'","`z_m'","`touse'","`b'","`V'","`clustvar'",`trim')
		}
		if "`clustvar'"=="" {
			// e(.) for basic robust
			local vce		robust
			local vcetype	Robust
		}
		else {
			// e(.) for cluster-robust; clustvar already defined
			local vce		cluster
			local vcetype	Robust
			local N_clust	=r(N_clust)
		}
		
		// store post objects
		mata: `A'.put(("N","post"),`r(N)')
		mata: `A'.put(("b","post"),st_matrix("r(b)"))
		mata: `A'.put(("V","post"),st_matrix("r(V)"))
		mata: `A'.put(("depvar","post"),"`yvar'")
		
		// for calling program
		ereturn clear
		mata: st_matrix("e(bmat)",st_matrix("r(b)"))
		mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("r(V)"))'))
		
		// store locals
		local list_local title yvar dvar y0 y0_m y1 y1_m vce vcetype teffect
		if "`model'"=="interactive" {
			local list_local `list_local' d d_m
		}
		else {
			local list_local `list_local' d0 d0_m d1 d1_m z z_m
		}
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// store scalars
		mata: `A'.put(("lltrim","scalar"),`r(lltrim)')
		mata: `A'.put(("ultrim","scalar"),`r(ultrim)')
		mata: `A'.put(("trim","scalar"),`trim')
		if "`clustvar'"~="" {
			mata: `A'.put(("N_clust","scalar"),`N_clust')
		}
		// additional estimation results
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yvar'")
		mata: st_local("shortstack",`eqn'.shortstack)
		// MSE
		mata: `A'.put(("`y0tilde'_mse","scalar"),return_result_item(`eqn',"`y0tilde'","MSE0","`rep'"))
		mata: `A'.put(("`y1tilde'_mse","scalar"),return_result_item(`eqn',"`y1tilde'","MSE1","`rep'"))
		// MSE folds
		mata: `A'.put(("`y0tilde'_mse_folds","matrix"),return_result_item(`eqn',"`y0tilde'","MSE0_folds","`rep'"))
		mata: `A'.put(("`y1tilde'_mse_folds","matrix"),return_result_item(`eqn',"`y1tilde'","MSE1_folds","`rep'"))
		// shortstack weights
		if "`spec'"=="ss" {
			mata: `A'.put(("`yvar'_ssw0","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights0","`rep'"))
			mata: `A'.put(("`yvar'_ssw1","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights1","`rep'"))
		}
		if "`model'"=="interactive" {
			// D eqn results
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			mata: st_local("shortstack",`eqn'.shortstack)
			// MSE
			mata: `A'.put(("`dtilde'_mse","scalar"),return_result_item(`eqn',"`dtilde'","MSE","`rep'"))
			// MSE folds
			mata: `A'.put(("`dtilde'_mse_folds","matrix"),return_result_item(`eqn',"`dtilde'","MSE_folds","`rep'"))
			// shortstack weights
			if "`spec'"=="ss" {
				mata: `A'.put(("`dvar'_ssw","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights","`rep'"))
			}
		}
		else {
			// D
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			mata: st_local("shortstack",`eqn'.shortstack)
			// MSE, D
			mata: `A'.put(("`d0tilde'_mse","scalar"),return_result_item(`eqn',"`d0tilde'","MSE0","`rep'"))
			mata: `A'.put(("`d1tilde'_mse","scalar"),return_result_item(`eqn',"`d1tilde'","MSE1","`rep'"))
			// MSE folds, D
			mata: `A'.put(("`d0tilde'_mse_folds","matrix"),return_result_item(`eqn',"`d0tilde'","MSE0_folds","`rep'"))
			mata: `A'.put(("`d1tilde'_mse_folds","matrix"),return_result_item(`eqn',"`d1tilde'","MSE1_folds","`rep'"))
			// shortstack weights, D
			if "`spec'"=="ss" {
				mata: `A'.put(("`dvar'_ssw0","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights0","`rep'"))
				mata: `A'.put(("`dvar'_ssw1","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights1","`rep'"))
			}
			// Z
			mata: `eqn' = (`mname'.eqnAA).get("`zvar'")
			mata: st_local("shortstack",`eqn'.shortstack)
			// MSE, Z
			mata: `A'.put(("`ztilde'_mse","scalar"),return_result_item(`eqn',"`ztilde'","MSE","`rep'"))
			// MSE folds, Z
			mata: `A'.put(("`ztilde'_mse_folds","matrix"),return_result_item(`eqn',"`ztilde'","MSE_folds","`rep'"))
			// shortstack weights, Z
			if "`spec'"=="ss" {
				mata: `A'.put(("`zvar'_ssw","matrix"),return_result_item(`eqn',"`shortstack'","ss_weights","`rep'"))
			}
		}
	
		mata: (`mname'.estAA).put(("`spec'","`rep'"),`A')
		
		// no longer needed
		foreach obj in `A' `eqn' {
			cap mata: mata drop `obj'
		}
	
	}
	else if "`replay'"=="" & "`medmean'"~="" {	// aggregate over resamples
		
		tempname b V bagg Vagg Vi
		tempname bvec sbvec bmed Vvec sVvec Vmed
		tempname nlltrim nultrim trimval
		tempvar esample
		tempname B
		
		// initialize
		mata: st_local("nameD",invtokens(`mname'.nameD))
		local K : word count `nameD'
		mata: st_local("nreps",strofreal(`mname'.nreps))
		mata: `B' = AssociativeArray()
		local isodd = mod(`nreps',2)
		local medrow = ceil(`nreps'/2)
		local N = 0
		
		// bvec a misnomer - usually a vector, but can be a matrix if multiple D variables
		mata: `bvec' = J(`nreps',`K',0)
		mata: `bagg' = J(1,`K',0)
		forvalues m=1/`nreps' {
			mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
			mata: `bvec'[`m',.] = `B'.get(("b","post"))
			// row/colnames etc. - need to do this only once
			if `m'==1 {
				mata: st_local("depvar",`B'.get(("depvar","post")))
				// retrieve locals; if empty, will be ""
				local list_local y0 y0_m y1 y1_m d d_m d0 d0_m d1 d1_m z z_m yvar dvar vce vcetype clustvar
				foreach obj in `list_local' {
					mata: st_local("`obj'",`B'.get(("`obj'","local")))
				}
				// retrieve scalars (as locals)
				local list_scalar
				if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
				foreach obj in `list_scalar' {
					mata: st_local("`obj'",strofreal(`B'.get(("`obj'","scalar"))))
				}
			}
			// possible that different estimation samples have different #obs
			qui count if `mname'_sample_`m'==1
			local N = `N' + r(N)
		}
		local N = round(`N'/`nreps')
		if "`medmean'"=="mn" {
			// mean beta
			mata: `bagg' = mean(`bvec')
			mata: st_matrix("`bagg'",`bagg')
		}
		else if "`medmean'"=="md" {
			// median beta
			forvalues k=1/`K' {
				// leave order of bvec unchanged
				mata: `sbvec' = sort(`bvec',`k')
				// mata: _sort(`bvec',`k')
				if `isodd' {
					mata: `bagg'[1,`k'] = `sbvec'[`medrow',`k']
				}
				else {
					mata: `bagg'[1,`k'] = (`sbvec'[`medrow',`k'] + `sbvec'[`medrow'+1,`k'])/2
				}
			}
			mata: st_matrix("`bagg'",`bagg')
		}
		else {
			di as err "_ddml_ate_late error - unrecognized option `medmean'"
			exit 198
		}
		
		mata: `Vagg' = J(`K',`K',0)
		mata: `Vvec' = J(`nreps',1,0)
		if "`medmean'"=="mn" {
			// harmonic mean
			// inefficient - does off-diagonals twice
			forvalues m=1/`nreps' {
				mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
				mata: `Vi' = `B'.get(("V","post"))
				forvalues j=1/`K' {
					forvalues k=1/`K' {
						// abs(.) needed?
						mata: `Vi'[`j',`k'] = `Vi'[`j',`k'] + abs((`bvec'[`m',`j'] - `bagg'[1,`j'])*(`bvec'[`m',`k'] - `bagg'[1,`k']))
					}
				}
				mata: `Vagg' = `Vagg' + 1:/`Vi'
			}
			mata: `Vagg' = `nreps' :/ `Vagg'
			mata: st_matrix("`Vagg'",`Vagg')
		}
		else if "`medmean'"=="md" {
			// median VCV
			// inefficient - does off-diagonals twice
			forvalues j=1/`K' {
				forvalues k=1/`K' {
					forvalues m=1/`nreps' {
						mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
						mata: `Vi' = `B'.get(("V","post"))
						mata: `Vvec'[`m'] = `Vi'[`j',`k']
					}
					// adjustment as per
					// https://docs.doubleml.org/stable/guide/resampling.html#repeated-cross-fitting-with-k-folds-and-m-repetition
					// (generalized to multiple D variables)
					mata: `Vvec' = `Vvec' + abs((`bvec'[.,`j'] :- `bagg'[1,`j']):*(`bvec'[.,`k'] :- `bagg'[1,`k']))
					// leave order of Vvec unchanged
					mata: `sVvec' = sort(`Vvec',1)
					// mata: _sort(`Vvec',1)
					if `isodd' {
						mata: `Vagg'[`j',`k'] = `sVvec'[`medrow',1]
					}
					else {
						mata: `Vagg'[`j',`k'] = (`sVvec'[`medrow',1] + `sVvec'[`medrow'+1,1])/2
					}
				}
			}
			mata: st_matrix("`Vagg'",`Vagg')
		}
		else {
			di as err "_ddml_ate_late error - unrecognized option `medmean'"
			exit 198
		}
		
		// count trim instances
		mata: `nlltrim' = 0
		mata: `nultrim' = 0
		forvalues m=1/`nreps' {
			mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
			mata: `nlltrim' = `nlltrim' + (`B'.get(("lltrim","scalar"))>0)
			mata: `nultrim' = `nultrim' + (`B'.get(("ultrim","scalar"))>0)
		}
		// retrieve from last rep (will be stored in all)
		mata: `trimval' = `B'.get(("trim","scalar"))
	
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),`bagg')
		mata: `A'.put(("V","post"),`Vagg')
		mata: `A'.put(("depvar","post"),"`depvar'")
		mata: `A'.put(("b_resamples","matrix"),`bvec')
		
		// store locals
		local list_local title yvar dvar y0 y1 vce vcetype teffect
		if "`model'"=="interactive" {
			local list_local `list_local' d
		}
		else {
			local list_local `list_local' d0 d1 z
		}
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// special case - "_m" subscript doesn't apply to mean/median over resamplings
		// so store without resample subscript
		foreach obj in title y0 y1 d d0 d1 z {
			mata: `A'.put(("`obj'_m","local"),"``obj''")
		}
		
		// store scalars
		local trim `trimval'	// hack, to fix
		local list_scalar nlltrim nultrim trim
		if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
		foreach obj in `list_scalar' {
			mata: `A'.put(("`obj'","scalar"),``obj'')
		}
		
		// store AA with median/mean results
		mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
		
		// no longer needed
		foreach obj in `A' `B' `bagg' `bvec' `sbvec' `Vagg' `Vvec' `sVvec' `Vi' `nlltrim' `nultrim' `trimval' {
			cap mata: mata drop `obj'
		}
			
	}
	else {
		// replay
				
		tempname B keys isscalar islocal ismatrix

		mata: `B' = AssociativeArray()
		mata: `B' = (`mname'.estAA).get(("`spec'","`rep'"))
		mata: `keys' = `B'.keys()
		mata: st_local("nentries",strofreal(rows(`keys')))
		mata: `isscalar'	= (`keys'[.,2] :== "scalar")
		mata: `islocal'		= (`keys'[.,2] :== "local")
		mata: `ismatrix'	= (`keys'[.,2] :== "matrix")
		
		tempname b V
		mata: st_matrix("`b'",`B'.get(("b","post")))
		mata: st_matrix("`V'",`B'.get(("V","post")))
		mata: st_local("N",strofreal(`B'.get(("N","post"))))
		mata: st_local("depvar",`B'.get(("depvar","post")))
		
		mata: st_local("yvar",`B'.get(("yvar","local")))
		mata: st_local("dvar",`B'.get(("dvar","local")))
		
		matrix rownames `b' = `depvar'
		matrix colnames `b' = `dvar'
	 	matrix colnames `V' = `dvar'
		matrix rownames `V' = `dvar'
		
		tempvar esample
		cap gen `esample' = `mname'_sample_`rep'
		if _rc>0 {
			// sample variable doesn't exist; ignore
			local esample
		}
		
		ereturn clear
		ereturn post `b' `V', depname(`depvar') obs(`N') esample(`esample')
		
		ereturn local cmd ddml
		ereturn local model `model'
		ereturn local rep `rep'
		ereturn local spec `spec'
		ereturn local tmname `mname'
		ereturn local teffect `teffect'
		
		// extract and post scalars, locals, matrices
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`isscalar'[`i']))
			if `topost' {
				mata: st_local("sname",substr(`keys'[`i',1],1,32))
				mata: st_numscalar("e(`sname')",`B'.get(`keys'[`i',.]))
			}
		}
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`islocal'[`i']))
			if `topost' {
				mata: st_local("lname",substr(`keys'[`i',1],1,32))
				mata: st_global("e(`lname')",`B'.get(`keys'[`i',.]))
			}
		}
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`ismatrix'[`i']))
			if `topost' {
				mata: st_local("tmname",substr(`keys'[`i',1],1,32))
				mata: st_matrix("e(`tmname')",`B'.get(`keys'[`i',.]))
			}
		}
		
		// no longer needed
		foreach obj in `B' `keys' `isscalar' `islocal' `ismatrix' {
			cap mata: mata drop `obj'
		}
		
		// display results
		di as text "`e(title)'"
		di as text "E[y|X,D=0]" _col(14) "= " as res "`e(y0_m)'" _c
		di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
		di as text "E[y|X,D=1]" _col(14) "= " as res "`e(y1_m)'"
		if "`e(model)'"=="interactive" {
			di as text "E[D|X]" _col(14)  "= " as res "`e(d_m)'"
		}
		else {
			di as text "E[D|X,Z=0]" _col(14)  "= " as res "`e(d0_m)'"
			di as text "E[D|X,Z=1]" _col(14)  "= " as res "`e(d1_m)'"
			di as text "E[Z|X]" _col(14)  "= " as res "`e(z_m)'"
		}
		ereturn display
		
		// report warning if clustered SEs requested but doesn't match clustered crossfitting
		mata: st_local("fclustvar",`mname'.fclustvar)
		if "`e(clustvar)'"~="" {
			if "`fclustvar'"=="" {
				di as res "Warning" as text ": crossfit folds do not respect cluster structure used for VCE."
			}
			else if "`fclustvar'"~="`e(clustvar)'" {
				di as res "Warning" as text ": cluster variable for VCE does not match cluster variable for crossfit folds."
			}
		}
	}
	
	// warn if any values trimmed
	if e(lltrim)>0 & e(lltrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(lltrim) " propensity scores trimmed to lower limit " e(trim) "."
	}
	if e(ultrim) & e(ultrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(ultrim) " propensity scores trimmed to upper limit " 1-e(trim) "."
	}
	// for mean/median over resamples
	if e(nlltrim)>0 & e(nlltrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(nlltrim) " resamples had propensity scores trimmed to lower limit " e(trim) "."
	}
	if e(ultrim) & e(ultrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(nultrim) " resamples had propensity scores trimmed to upper limit " 1-e(trim) "."
	}
	
end


********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

void ATE(   
			string scalar teffect,    // ATE, ATET or ATEU
			string scalar yvar,       // Y
			string scalar dvar,       // D
			string scalar y0tilde,    // E[Y|X,D=0]
			string scalar y1tilde,    // E[Y|X,D=1]
			string scalar dtilde,     // E[D|X]
			string scalar sample,     // sample
			string scalar outate,     // output: name of matrix to store b
			string scalar outatese,   // output: name of matrix to store V
			string scalar clustvar,   //
			string scalar foldvar, 	  //
			real scalar trim          // trim the propensity score
			)
{
	st_view(my_d0x,.,y0tilde,sample)
	st_view(my_d1x,.,y1tilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)
	st_view(fid,.,foldvar,sample)
	// copy since we may trim it
	md_x=st_data(.,dtilde,sample)
	if (clustvar!="") {
		st_view(clustid,.,clustvar,sample)
		clustid_uni=uniqrows(clustid)
		nclust = rows(clustid_uni)
	}

	n = rows(y)
	
	// first a vector
	lltrim = md_x :< trim
	ultrim = md_x :> (1-trim)
	// trim
	md_x = md_x :* (1:-lltrim) + trim * lltrim
	md_x = md_x :* (1:-ultrim) + (1-trim) * ultrim
	// now a scalar
	lltrim = sum(lltrim)
	ultrim = sum(ultrim)

	// psi = psi_b + psi_a*theta, e.g. equation 5.62
	if (teffect=="ATE") {
		psi_b  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x 
		psi_a  = J(n,1,-1) 
	}
	else {
		// ATET and ATEU
		// calculate mean of d by fold
		fid_uni = uniqrows(fid)
		folds = rows(fid_uni)
		p_hat = J(n,1,.)
		for (j=1;j<=folds;j++) {
			k=fid_uni[j,1]
			sel = selectindex(fid:==k)
			meank = mean(d[sel])
			p_hat[sel] = J(length(sel), 1, meank)
		}
		if (teffect=="ATET") {
			psi_b = (d :* (y :- my_d0x) :/ p_hat) :-  md_x :* (1 :- d) :* (y :- my_d0x) :/ (p_hat :*(1 :- md_x)) 
			psi_a = -d :/ p_hat
		}
		else if (teffect=="ATEU") {
			psi_b = ((1:-d) :* (y :- my_d1x) :/ (1:-p_hat)) :- (1:-md_x) :* d :* (y :- my_d1x) :/ ((1:-p_hat) :* md_x) 
			psi_a = (1:-d) :/ (1:-p_hat)
		}
		else {
			errprintf("internal ddml error - teffect estimator %s not defined\n",teffect)
			exit(499)
		}
	}
	theta = -mean(psi_b) / mean(psi_a)
	psi = psi_a :* theta :+ psi_b

	if (clustvar=="") {
		V =  mean(psi:^2) / (mean(psi_a):^2) / n
	}
	else {
		gamma = 0
		jhat = 0
		for (i=1;i<=nclust;i++) {
			psi_c = select(psi,clustid:==clustid_uni[i,1])
			psi_a_c = select(psi_a,clustid:==clustid_uni[i,1])
			gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
			jhat = jhat :+  1/nclust :* sum(psi_a_c)
		}
		V = gamma / jhat:^2 / nclust
		st_numscalar("r(N_clust)",nclust)
	}

	st_numscalar("r(N)",n)
	st_matrix("r(b)",theta)
	st_matrix("r(V)",V)
	st_numscalar("r(lltrim)",lltrim)
	st_numscalar("r(ultrim)",ultrim)
}

void LATE(  string scalar yvar,      // Y
            string scalar dvar,      // D
            string scalar zvar,      // Z
            string scalar y0tilde,   // E[Y|X,Z=0]
            string scalar y1tilde,   // E[Y|X,Z=1]
            string scalar d0tilde,   // E[D|X,Z=0]
            string scalar d1tilde,   // E[D|X,Z=1]
            string scalar ztilde,    // E[Z|X]
            string scalar sample,    // sample
            string scalar outlate,   // output: name of matrix to store b
            string scalar outlatese,  // output: name of matrix to store V
            string scalar clustvar,
 			real scalar trim          // trim the propensity score
           )
{
    st_view(my_z0x,.,y0tilde,sample)
    st_view(my_z1x,.,y1tilde,sample)
    st_view(md_z0x,.,d0tilde,sample)
    st_view(md_z1x,.,d1tilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)
    st_view(z,.,zvar,sample)
	// copy since we may trim it
    mz_x=st_data(.,ztilde,sample)
    if (clustvar!="") {
		st_view(clustid,.,clustvar,sample)
		clustid_uni=uniqrows(clustid)
		nclust = rows(clustid_uni)
	}

    n = rows(y)

	// first a vector
	lltrim = mz_x :< trim
	ultrim = mz_x :> (1-trim)
	// trim
	mz_x = mz_x :* (1:-lltrim) + trim * lltrim
	mz_x = mz_x :* (1:-ultrim) + (1-trim) * ultrim
	// now a scalar
	lltrim = sum(lltrim)
	ultrim = sum(ultrim)
	
    psi_b =  z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x 
    psi_a =  -(z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x)

	theta = -mean(psi_b) / mean(psi_a)
	psi = psi_a :* theta :+ psi_b

	if (clustvar=="") {
		V =  mean(psi:^2) :/ mean(psi_a):^2 :/ n
	}
	else {
		gamma = 0
		jhat = 0
		for (i=1;i<=nclust;i++) {
			psi_c = select(psi,clustid:==clustid_uni[i,1])
			psi_a_c = select(psi_a,clustid:==clustid_uni[i,1])
			gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
			jhat = jhat :+  1/nclust :* sum(psi_a_c)
		}
		V = gamma / jhat:^2 / nclust
		st_numscalar("r(N_clust)",nclust)
	}

    st_numscalar("r(N)",n)
    st_matrix("r(b)",theta)
    st_matrix("r(V)",V)
	st_numscalar("r(lltrim)",lltrim)
	st_numscalar("r(ultrim)",ultrim)
}

end
