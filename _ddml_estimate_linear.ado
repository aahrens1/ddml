*! ddml v1.1
*! last edited: 6 jan 2023
*! authors: aa/ms

program _ddml_estimate_linear, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ROBust				///
								CLUster(varname)	///
								NOConstant			/// suppress constant in estimation
								SHOWConstant		/// display constant in summary table
								vce(string)			///
								allcombos			/// estimate and show all combinations (dfn changed below)
								NOTable				/// suppress summary table
								clear				/// deletes all tilde-variables (to be implemented)
								spec(string)		/// specification to post/display
								REP(string)			/// resampling iteration to post/display or mean/median
								replay				/// model has been estimated, just display results
								debug				///
								* ]
	
	if "`debug'"==""	local qui qui
	
	marksample touse
	
	// consflag
	local consflag = ("`noconstant'"=="")
	local showconsflag = ("`showconstant'"~="" & `consflag')	// ignored if nocons
	// replay existing results
	local replayflag = "`replay'"~=""
	// display summary table
	local tableflag = "`notable'"==""
	// request estimation/reporting of all combinations
	local doallcombos = "`allcombos'"~=""
	// remaining flags
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	// reestimation necessary unless replay specified
	if `replayflag' {
		// estimated macro =0/1 indicating estimation results exist
		mata: st_local("estimated", strofreal(`mname'.estimated))
		// initial ncombos; will be 0 if all combos not (yet) estimated
		mata: st_local("ncombos", strofreal(`mname'.ncombos))
		// error checks
		if `estimated'==0 {
			di as err "error - replay specified but model not yet estimated"
			exit 198
		}
		if `ncombos'==0 & "`spec'"~="" & real("`spec'")<. {
			di as err "error - spec(`spec') not available; add 'allcombos' to estimate all combinations"
			di as err "add 'replay' to retrieve one of the available estimations stored in memory"
			exit 198
		}
	}
	else {
		// error checks
		if `doallcombos'==0 & "`spec'"~="" & real("`spec'")<. {
			di as err "error - spec(`spec') not available; add 'allcombos' to estimate all combinations"
			exit 198
		}
		mata: clear_model_estimation(`mname')
		local estimated = 0
		local ncombos = 0
	}
	
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

	** standard errors
	// local vce is the argument to the Stata option vce(.)
	if "`robust'"!=""	local vce robust
	if "`cluster'"~=""	local vce cluster `cluster'
	
	if ~`crossfitted' {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}

	// default spec is ss if available, otherwise mse	
	if "`spec'"=="" & `ssflag' {
		local spec "ss"
	}
	else if "`spec'"=="" {
		local spec "mse"
	}
	// allowed forms of spec and rep
	if "`spec'"=="shortstack"	local spec ss
	if "`spec'"=="minmse"		local spec mse
	if "`rep'"=="mean"			local rep mn
	if "`rep'"=="median"		local rep md

	// if rep not specified, default is rep=1 when nreps==1; md if nreps>1
	if "`rep'"=="" & `nreps'>1 {
		local rep md
	}
	else if "`rep'"=="" & `nreps'==1 {
		local rep 1
	}
	
	// checks
	if "`spec'"~="ss" & "`spec'"~="mse" & real("`spec'")==. {
		di as err "error - invalid spec(`spec')"
		exit 198
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
			di as err "rep() cannot be larger than `nreps' in current model specification"
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
			di as err "error - shortstack not available for all equations"
			exit 198
		}
	}

	// shortstack names
	if `ssflag' {
		local Yss `nameY'_ss
		foreach var in `nameD' {
			local Dss `Dss' `var'_ss
			local DHss `DHss' `var'_ss_h
		}
		foreach var in `nameZ' {
			local Zss `Zss' `var'_ss
		}
	}
	
	************* ESTIMATE ************
	
	if `estimated'==0 {
		// enter if no estimates exist
	
		// Loop over resamples and estimate/save the min mse and ss model for each
		forvalues m=1/`nreps' {
			
			// reset locals
			local Yopt
			local Dopt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
			foreach var in `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
				local Dopt `Dopt' `oneDopt'
				// DHopt is stored in list Zopt
				if "`model'"=="fiv" {
					mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
					local Zopt `Zopt' `oneDHopt'
				}
			}
			// nameZ is empty for fiv model
			foreach var in `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneZopt",return_learner_item(`eqn',"opt","`m'"))
				local Zopt `Zopt' `oneZopt'
			}
			
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			if "`model'"=="fiv" {
				local dlist
				local zlist
				forvalues j = 1/`numeqnD' {
					tempvar zvar`j' 
					tempvar dx`j'
					local dh : word `j' of `z'
					local dt : word `j' of `d'
					local dd : word `j' of `nameD'
					qui gen double `zvar`j'' = `dt'_`m'-`dh'_`m' // E[D|ZX]-E[D|X] = instrument
					qui gen double `dx`j'' = `dd'-`dh'_`m' // D-E[D|X] = endogenous regressor
					local dlist `dlist' `dx`j''
					local zlist `zlist' `zvar`j''
				}
				local dvtnames `d'
				local zvtnames `z'
				local d `dlist'
				local z `zlist'
				local norep norep
			}
			local title Min MSE DDML model`stext'
			`qui' _ddml_reg if `mname'_sample_`m' & `touse',					///
					`noconstant' vce(`vce')										///
					y(`Yopt') yname(`nameY')									///
					d(`Dopt') dnames(`nameD') dvtnames(`dvtnames')		 		///
					z(`Zopt') znames(`nameZ') zvtnames(`zvtnames')				///
					mname(`mname') spec(mse) rep(`m') title(`title') `norep'
			
			// estimate shortstack for this rep
			if `ssflag' {
				if "`model'"=="fiv" {
					local dlist
					local zlist
					forvalues j = 1/`numeqnD' {
						tempvar zvar`j' 
						tempvar dx`j'
						local dh : word `j' of `DHss'
						local dt : word `j' of `Dss'
						local dd : word `j' of `nameD'
						qui gen double `zvar`j'' = `dt'_`m'-`dh'_`m' // E[D|ZX]-E[D|X] = instrument
						qui gen double `dx`j'' = `dd'-`dh'_`m' // D-E[D|X] = endogenous regressor
						local dlist `dlist' `dx`j''
						local zlist `zlist' `zvar`j''
					}
					local dvtnames `Dss'
					local zvtnames `DHss'
					local d `dlist'
					local z `zlist'
					local norep norep
				}
				else {
					local d `Dss'
					local z `Zss'
				}
				
				local title "Shortstack DDML model`stext'"
				`qui' _ddml_reg if `mname'_sample_`m' & `touse',					///
						`noconstant' vce(`vce')										///
						y(`Yss') yname(`nameY')										///
						d(`d') dnames(`nameD') dvtnames(`dvtnames') 				///
						z(`z') znames(`nameZ') zvtnames(`zvtnames')					///
						mname(`mname') spec(ss) rep(`m') title(`title') `norep'
			}
		}
		
		// have looped over reps to get each optimal model and shortstack per rep
		// now aggregate over reps to get mean/median
		if `nreps' > 1 {
			`qui' _ddml_reg, mname(`mname') spec(mse) medmean(mn) title("Mean over min-mse specifications") `noconstant'
			`qui' _ddml_reg, mname(`mname') spec(mse) medmean(md) title("Median over min-mse specifications") `noconstant'
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean over `nreps' resamples)"
				`qui' _ddml_reg, mname(`mname') spec(ss) medmean(mn) title(`title') `noconstant'
				local title "Shortstack DDML model (median over `nreps' resamples)"
				`qui' _ddml_reg, mname(`mname') spec(ss) medmean(md) title(`title') `noconstant'
			}
		}
		
		// Estimation of med/med/shortstack complete, mark as estimated.
		mata: `mname'.estimated = 1
		// (re-)set estimated
		mata: st_local("estimated", strofreal(`mname'.estimated))
	}
	
	******************************************
	
	if `ncombos' {
		// all combos have already been estimated, so just recover matrices
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
	else if `doallcombos' {
		// allcombos need to be estimated
		// get varlists
		_ddml_make_varlists, mname(`mname')
		local yvars `r(yvars)'
		local dvars `r(dvars)'
		local zvars `r(zvars)'
		local zpos	`r(zpos_start)'
	
		// obtain all combinations
		_ddml_allcombos `yvars' - `dvars' - `zvars' ,	///
			`debug'										///
			zpos(`zpos')		 						///
			addprefix("")
		
		local ncombos = r(ncombos)
		local tokenlen = `ncombos'*2
		local ylist `r(ystr)'
		local Dlist `r(dstr)'
		local Zlist `r(zstr)' 
		
		tempname nmat bmat semat
		mata: `nmat' = J(`ncombos',3,"")
		mata: `bmat' = J(`ncombos'*`nreps',`numeqnD'+`consflag',.)
		mata: `semat' = J(`ncombos'*`nreps',`numeqnD'+`consflag',.)
		
		// simplest if put into a Mata string matrix
		tokenize `ylist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',1] = strtrim("``idx''")
		}
		tokenize `Dlist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',2] = strtrim("``idx''")
		}
		tokenize `Zlist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',3] = strtrim("``idx''")
		}
		
		forvalues m=1/`nreps' {
			
			// reset locals
			local Yopt
			local Dopt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
			foreach var in `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
				local Dopt `Dopt' `oneDopt'
				// DHopt is stored in list Zopt
				if "`model'"=="fiv" {
					mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
					local Zopt `Zopt' `oneDHopt'
				}
			}
			// nameZ is empty for fiv model
			foreach var in `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneZopt",return_learner_item(`eqn',"opt","`m'"))
				local Zopt `Zopt' `oneZopt'
			}
			
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			forvalues i = 1/`ncombos' {
				mata: st_local("y",`nmat'[`i',1])
				mata: st_local("d",`nmat'[`i',2])
				mata: st_local("z",`nmat'[`i',3])
				// check if opt for this resample; note === for D so order doesn't matter
				local isopt
				local isYopt : list Yopt == y
				local isDopt : list Dopt === d
				local isZopt : list Zopt === z
				local title "DDML model, specification `i'`stext'"
				if `isYopt' & `isDopt' & `isZopt' {
					local optspec`m' = `i'
					local isopt *
					local title Min MSE `title'
					// save in AA
					mata: (`mname'.estAA).put(("optspec","`m'"),"`i'")
				}
				if "`model'"=="fiv" {
					local dlist
					local zlist
					forvalues j = 1/`numeqnD' {
						tempvar zvar`j' 
						tempvar dx`j'
						local dh : word `j' of `z'
						local dt : word `j' of `d'
						local dd : word `j' of `nameD'
						qui gen double `zvar`j'' = `dt'_`m'-`dh'_`m' // E[D|ZX]-E[D|X] = instrument
						qui gen double `dx`j'' = `dd'-`dh'_`m' // D-E[D|X] = endogenous regressor
						local dlist `dlist' `dx`j''
						local zlist `zlist' `zvar`j''
					}
					local dvtnames `d'
					local zvtnames `z'
					local d `dlist'
					local z `zlist'
					local norep norep
				}
				`qui' _ddml_reg if `mname'_sample_`m' & `touse',				///
						`noconstant' vce(`vce')									///
						y(`y') yname(`nameY')									///
						d(`d') dnames(`nameD') dvtnames(`dvtnames')		 		///
						z(`z') znames(`nameZ') zvtnames(`zvtnames')				///
						mname(`mname') spec(`i') rep(`m') title(`title') `norep'
				
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
			}
			
			if `ssflag' {
				if "`model'"=="fiv" {
					local dlist
					local zlist
					forvalues j = 1/`numeqnD' {
						tempvar zvar`j' 
						tempvar dx`j'
						local dh : word `j' of `DHss'
						local dt : word `j' of `Dss'
						local dd : word `j' of `nameD'
						qui gen double `zvar`j'' = `dt'_`m'-`dh'_`m' // E[D|ZX]-E[D|X] = instrument
						qui gen double `dx`j'' = `dd'-`dh'_`m' // D-E[D|X] = endogenous regressor
						local dlist `dlist' `dx`j''
						local zlist `zlist' `zvar`j''
					}
					local dvtnames `Dss'
					local zvtnames `DHss'
					local d `dlist'
					local z `zlist'
					local norep norep
				}
				else {
					local d `Dss'
					local z `Zss'
				}
				local title "Shortstack DDML model`stext'"
				`qui' _ddml_reg if `mname'_sample_`m' & `touse',					///
						`noconstant' vce(`vce')										///
						y(`Yss') yname(`nameY')										///
						d(`d') dnames(`nameD') dvtnames(`dvtnames') 				///
						z(`z') znames(`nameZ') zvtnames(`zvtnames')					///
						mname(`mname') spec(ss) rep(`m') title(`title') `norep'
			}
	
		}
	
		// estimation of all combos complete; ncombos > 0 indicates all combos estimated
		mata: `mname'.ncombos = `ncombos'
		mata: (`mname'.estAA).put(("nmat","all"),`nmat')
		mata: (`mname'.estAA).put(("bmat","all"),`bmat')
		mata: (`mname'.estAA).put(("semat","all"),`semat')
	}

	************** REPORT RESULTS **************
	
	*** Results ***
	// optional table of results
	if `tableflag' {
		di
		di as text "DDML estimation results:"
		di as text "spec  r" %14s "Y learner" _c
		forvalues j=1/`numeqnD' {
			di as text %14s "D learner" %10s "b" %10s "SE" _c
		}
		if `showconsflag' {
			di as text %10s "_cons" %10s "SE" _c
		}
		if "`model'"=="fiv" {
			forvalues j=1/`numeqnD' {
				di as text %14s "DH learner" _c
			}
		}
		forvalues j=1/`numeqnZ' {
			di as text %14s "Z learner" _c
		}
		di
		forvalues m=1/`nreps' {
			if `doallcombos' {
				// all combos available, so loop through
				forvalues i=1/`ncombos' {
					mata: st_local("yt",abbrev(`nmat'[`i',1],13))
					mata: st_local("dtlist",invtokens(abbrev(tokens(`nmat'[`i',2]),13)))
					mata: st_local("ztlist",`nmat'[`i',3])
					if "`ztlist'"~="" {
						mata: st_local("ztlist",invtokens(abbrev(tokens(`nmat'[`i',3]),13)))
					}
					if "`optspec`m''"=="`i'" {
						di "*" _c
					}
					else {
						di " " _c
					}
					local specrep `: di %3.0f `i' %3.0f `m''
					// pad out to 6 spaces
					local specrep = (6-length("`specrep'"))*" " + "`specrep'"
					local rcmd stata ddml estimate, mname(`mname') spec(`i') rep(`m') replay notable `noconstant'
					di %6s "{`rcmd':`specrep'}" _c
					di as res %14s "`yt'" _c
					forvalues j=1/`numeqnD' {
						local vt : word `j' of `dtlist'
						mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
						mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
						di as res %14s "`vt'" _c
						di as res %10.3f `b' _c
						local pse (`: di %6.3f `se'')
						di as res %10s "`pse'" _c
					}
					if `showconsflag' {
						// set j by hand; cons is in last column
						local j = `numeqnD' + 1
						mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
						mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
						di as res %10.3f `b' _c
						local pse (`: di %6.3f `se'')
						di as res %10s "`pse'" _c
					}
					forvalues j=1/`numeqnZ' {
						local vt : word `j' of `ztlist'
						di as res %14s "`vt'" _c
					}
					if "`model'"=="fiv" {
						forvalues j=1/`numeqnD' {
							local vt : word `j' of `ztlist'
							di as res %14s "`vt'" _c
						}
					}
					di
				}
			}
			else {
				// only mse/ss specs available/reported
				`qui' _ddml_reg, mname(`mname') spec(mse) rep(`m') replay `noconstant'
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep `: di "opt" %3.0f `m''
				// pad out to 6 spaces
				local specrep = " " + "`specrep'"
				local rcmd stata ddml estimate, mname(`mname') spec(mse) rep(`m') replay notable `noconstant'
				di %6s "{`rcmd':`specrep'}" _c
				mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
				mata: st_local("yt",return_learner_item(`eqn',"opt","`m'"))
				di as res %14s "`yt'" _c
				forvalues j=1/`numeqnD' {
					local dd : word `j' of `nameD'
					mata: `eqn' = (`mname'.eqnAA).get("`dd'")
					mata: st_local("vt",return_learner_item(`eqn',"opt","`m'"))
					di as res %14s "`vt'" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if `showconsflag' {
					// set j by hand; cons is in last column
					local j = `numeqnD' + 1
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				forvalues j=1/`numeqnZ' {
					local vt : word `j' of `ztlist'
					di as res %14s "`vt'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						local vt : word `j' of `ztlist'
						di as res %14s "`vt'" _c
					}
				}
				di
			}
			if `ssflag' {
				`qui' _ddml_reg, mname(`mname') spec(ss) rep(`m') replay `noconstant'
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep `: di "ss" %3.0f `m''
				// pad out to 6 spaces
				local specrep = "  " + "`specrep'"
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`m') replay notable `noconstant'
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ss]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if `showconsflag' {
					// set j by hand; cons is in last column
					local j = `numeqnD' + 1
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ss]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
					di as res %14s "[ss]" _c
				}
				di
			}
		}
		if `doallcombos' {
			di as res "*" _c
		}
		else {
			di as text "opt" _c
		}
		di as text " = minimum MSE specification for that resample."
	}

	if `nreps' > 1 & `tableflag' {
		di
		di as text "Mean/med.   Y learner" _c
		forvalues j=1/`numeqnD' {
			di as text %14s "D learner" %10s "b" %10s "SE" _c
		}
		if "`model'"=="fiv" {
			forvalues j=1/`numeqnD' {
				di as text %14s "DH learner" _c
			}
		}
		forvalues j=1/`numeqnZ' {
			di as text %14s "Z learner" _c
		}
		di
		foreach medmean in mn md {
			** mean and median over mse
			// force noconstant with mean/median
			`qui' _ddml_reg, mname(`mname') spec(mse) rep(`medmean') replay noconstant
			tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
			mat `btemp' = e(b)
			mat `Vtemp' = e(V)
			local specrep `: di "mse" %3s "`medmean'"'
			// pad out to 6 spaces
			local specrep = " " + "`specrep'"
			// force noconstant with mean/median
			local rcmd stata ddml estimate, mname(`mname') spec(mse) rep(`medmean') replay notable `noconstant'
			di %6s "{`rcmd':`specrep'}" _c
			di as res %14s "[min-mse]" _c
			forvalues j=1/`numeqnD' {
				di as res %14s "[mse]" _c
				di as res %10.3f el(`btemp',1,`j') _c
				local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
				di as res %10s "`pse'" _c
			}
			if "`model'"=="fiv" {
				forvalues j=1/`numeqnD' {
					di as res %14s "[mse]" _c
				}
			}
			forvalues j=1/`numeqnZ' {
				di as res %14s "[mse]" _c
			}
			di
			if `ssflag' {
				// force noconstant with mean/median
				`qui' _ddml_reg, mname(`mname') spec(ss) rep(`medmean') replay noconstant
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep `: di "ss" %4s "`medmean'"'
				// pad out to 6 spaces
				local specrep = " " + "`specrep'"
				// force noconstant with mean/median
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`medmean') replay notable noconstant
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ss]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ss]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
					di as res %14s "[ss]" _c
				}
				di
			}
		}
	}
	
	di
	if ("`rep'"=="mn" | "`rep'"=="md") {
		// force noconstant with mean/median
		_ddml_reg, mname(`mname') spec(`spec') rep(`rep') replay noconstant
	}
	else {
		_ddml_reg, mname(`mname') spec(`spec') rep(`rep') replay `noconstant'
	}
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
	
	// temp Mata object no longer needed
	foreach obj in `eqn' `nmat' `bmat' `semat' {
		cap mata: mata drop `obj'
	}


end

// adds model name prefixes to list of varnames
program define add_prefix, sclass
	syntax [anything] , prefix(name)

	// anything is a list of to-be-varnames that need prefix added to them
	foreach vn in `anything' {
		local vnames `vnames' `prefix'`vn' 
	}
	
	sreturn local vnames `vnames'
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

// does OLS/IV and reports with substitute yname and dnames
program define _ddml_reg, eclass
	syntax [anything] [if] [in] , [								///
				y(name) yname(name)								///
				d(namelist) dnames(namelist) dvtnames(namelist)	///
				z(namelist) znames(namelist) zvtnames(namelist)	///
				mname(name)										///
				spec(string) rep(string)						///
				vce(string)										///
				title(string)									///
				medmean(string)									///
				NOREP											///
				NOConstant										///
				replay											///
				*												///
				]

	// consflag
	local consflag = ("`noconstant'"=="")
	mata: st_local("model",`mname'.model)
	local fivflag	= "`model'"=="fiv"

	if "`replay'"=="" & "`medmean'"=="" {	// estimate from scratch
		
		marksample touse
		
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		// default vtilde names
		if "`dvtnames'"=="" {
			local dvtnames `d'
		}
		if "`zvtnames'"=="" {
			local zvtnames `z'
		}

		// add resample suffixes
		// y always gets a suffix
		local y_m `y'_`rep'
		if "`norep'"=="" {	
			add_suffix `d', suffix("_`rep'")
			local d_m `s(vnames)'
			add_suffix `z', suffix("_`rep'")
			local z_m `s(vnames)'
		}
		else {
			local d_m `d'
			local z_m `z'
		}
		
		// estimate
		if "`z_m'"=="" {
			qui reg `y_m' `d_m'         if `touse', vce(`vce') `noconstant' `options'
		}
		else {
			// old-style regress syntax: put IVs in parentheses
			qui reg `y_m' `d_m' (`z_m') if `touse', vce(`vce') `noconstant' `options'
		}
		tempname b V
		mat `b' = e(b)
		mat `V' = e(V)
		local N = e(N)

		local vce		`e(vce)'
		local vcetype	`e(vcetype)'
		local clustvar	`e(clustvar)'
		local N_clust	=e(N_clust)
		if `fivflag' {
			local d		`dvtnames'
			add_suffix	`dvtnames', suffix("_`rep'")
			local d_m	`s(vnames)'
			local dh	`zvtnames'
			add_suffix	`zvtnames', suffix("_`rep'")
			local dh_m	`s(vnames)'
			// clear locals (IV info is stored as dh not z)
			local z
			local z_m
		}
		
		// store post objects
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),st_matrix("`b'"))
		mata: `A'.put(("V","post"),st_matrix("`V'"))
		mata: `A'.put(("depvar","post"),"`yname'")
		
		// for calling program
		ereturn clear
		mata: st_matrix("e(bmat)",st_matrix("`b'"))
		mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("`V'"))'))
		
		// store locals
		local list_local title y y_m d d_m yname dnames vce vcetype
		if "`model'"=="iv" {
			local list_local `list_local' z z_m
		}
		else if "`model'"=="fiv" {
			local list_local `list_local' z z_m dh dh_m
		}
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// store scalars
		local list_scalar
		if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
		foreach obj in `list_scalar' {
			mata: `A'.put(("`obj'","scalar"),``obj'')
		}
		
		// additional estimation results
		ereturn scalar resample = `rep'
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yname'")
		// MSE
		mata: `A'.put(("`y'_mse","scalar"),return_result_item(`eqn',"`y'","MSE","`rep'"))
		// MSE folds
		mata: `A'.put(("`y'_mse_folds","matrix"),return_result_item(`eqn',"`y'","MSE_folds","`rep'"))
		// ss weights
		if "`spec'"=="ss" {
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: `A'.put(("`yname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'","ss_weights","`rep'"))
		}

		// D eqn results - uses vtnames
		local numeqnD	: word count `dnames'
		forvalues i=1/`numeqnD' {
			local dname : word `i' of `dnames'
			local vtilde : word `i' of `dvtnames'
			local vtilde_h : word `i' of `zvtnames'
			// remove the trailing "_h" so that the AA lookup uses the learner name
			local vtilde_h = substr("`vtilde_h'",1,strlen("`vtilde_h'")-2)
			mata: `eqn' = (`mname'.eqnAA).get("`dname'")
			// MSE
			mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
			// MSE folds
			mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
			// ss weights
			if "`spec'"=="ss" {
				mata: st_local("shortstack", `eqn'.shortstack)
				mata: `A'.put(("`dname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'","ss_weights","`rep'"))
			}
			if `fivflag' {
				// MSE
				mata: `A'.put(("`vtilde_h'_mse","scalar"),return_result_item(`eqn',"`vtilde_h'","MSE_h","`rep'"))
				// MSE folds
				mata: `A'.put(("`vtilde_h'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde_h'","MSE_h_folds","`rep'"))
				// ss weights
				if "`spec'"=="ss" {
						mata: st_local("shortstack", `eqn'.shortstack)
						mata: `A'.put(("`dname'_h_ssw","matrix"), return_result_item(`eqn',"`shortstack'","ss_weights_h","`rep'"))
				}
			}
		}
		if `fivflag'==0 {
			// Z eqn results; fiv won't enter
			local numeqnZ	: word count `znames'
			forvalues i=1/`numeqnZ' {
				local zname : word `i' of `znames'
				local vtilde : word `i' of `z'
				mata: `eqn' = (`mname'.eqnAA).get("`zname'")
				// MSE
				mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
				// MSE folds
				mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
				// ss weights
				if "`spec'"=="ss" {
					mata: st_local("shortstack", `eqn'.shortstack)
					mata: `A'.put(("`zname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'","ss_weights","`rep'"))
				}
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
		tempname bvec brow sbvec bmed Vvec sVvec Vmed
		tempvar esample
		tempname B
		
		// initialize
		mata: st_local("nameD",invtokens(`mname'.nameD))
		// don't aggregate the constant
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
			mata: `brow' = `B'.get(("b","post"))
			// don't aggregate the constant
			mata: `bvec'[`m',.] = `brow'[1,(1..(cols(`brow')-`consflag'))]
			// row/colnames etc. - need to do this only once
			if `m'==1 {
				mata: st_local("depvar",`B'.get(("depvar","post")))
				// retrieve locals; if empty, will be ""
				local list_local y d dh z yname dnames vce vcetype clustvar
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
				mata: `sbvec' = sort(`bvec',`k')
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
			di as err "_ddml_reg error - unrecognized option `medmean'"
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
				// don't aggregate the constant
				mata: `Vi' = `Vi'[(1::(cols(`Vi')-`consflag')),(1..(cols(`Vi')-`consflag'))]
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
					mata: `sVvec' = sort(`Vvec',1)
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
			di as err "_ddml_reg error - unrecognized option `medmean'"
			exit 198
		}
	
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),`bagg')
		mata: `A'.put(("V","post"),`Vagg')
		mata: `A'.put(("depvar","post"),"`depvar'")
		
		// store locals
		local list_local title y d yname dnames vce vcetype
		if "`model'"=="iv" {
			local list_local `list_local' z
		}
		else if "`model'"=="fiv" {
			local list_local `list_local' z dh
		}
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// store scalars
		local list_scalar
		if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
		foreach obj in `list_scalar' {
			mata: `A'.put(("`obj'","scalar"),``obj'')
		}
		// special case - "_m" subscript doesn't apply to mean/median over resamplings
		// so store without resample subscript
		foreach obj in y d dh z {
			mata: `A'.put(("`obj'_m","local"),"``obj''")
		}
		// special case - vector of betas available only for mean/median
		mata: `A'.put(("b_resamples","matrix"),`bvec')
		
		// store AA with median/mean results
		mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
		
		// no longer needed
		foreach obj in `A' `B' `bagg' `bvec' `brow' `sbvec' `Vagg' `Vvec' `sVvec' `Vi' {
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
		
		mata: st_local("yname",`B'.get(("yname","local")))
		mata: st_local("dnames",`B'.get(("dnames","local")))
		
		if `consflag' {
			// will be empty if no constant
			local consname "_cons"
		}
		
		matrix rownames `b' = `depvar'
		matrix colnames `b' = `dnames' `consname'
	 	matrix colnames `V' = `dnames' `consname'
		matrix rownames `V' = `dnames' `consname'
		
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
		ereturn local tmname `mname'	// temporary mname
		
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
		di as text "y-E[y|X]" _col(11) "= " as res "`e(y_m)'" _c
		di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
		if "`e(model)'"~="fiv" {
			di as text "D-" _c
		}
		di as text "E[D|X,Z]" _col(11)  "= " as res "`e(d_m)'"
		if "`e(model)'" == "iv" {
			di as text "Z-E[Z|X]" _col(11) "= " as res "`e(z_m)'"
		}
		else if "`e(model)'" == "fiv" {
			di as text "E[D|X]" _col(11) "= " as res "`e(dh_m)'"
		}
		if "`e(model)'" == "fiv" {
			di as text "Orthogonalised D = D - E[D|X]; optimal IV = E[D|X,Z] - E[D|X]."
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

end

program define _ddml_make_varlists, rclass

	syntax [anything], mname(name)
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()

	// locals used below
	mata: st_local("model",`mname'.model)
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'
	
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			tempname Dt_list
			local `Dt_list' `vtlistD'
			local Dorthog_lists `Dorthog_lists' `Dt_list'
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
			tempname DHt_list
			local `DHt_list' `vtlistD_h'
			local DHorthog_lists `DHorthog_lists' `DHt_list'
		}
	}
	
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			tempname Zt_list
			local `Zt_list' `vtlistZ'
			local Zorthog_lists `Zorthog_lists' `Zt_list'
		}
	}

	local dash
	foreach vl in `Dorthog_lists' {
		local Dtilde `Dtilde' `dash' ``vl''
		local dash -
	}
	local dash
	foreach vl in `Zorthog_lists' {
		local Ztilde `Ztilde' `dash' ``vl''
		local dash -
	}
	local dash
	foreach vl in `DHorthog_lists' {
		local DHtilde `DHtilde' `dash' ``vl''
		local dash -
	}

	// clear from Mata
	mata: mata drop `eqn'
	
	return scalar dpos_end = `numeqnD' + 1
	return scalar dpos_start = 2
	if (`numeqnZ'>0) {
		return scalar zpos_start = `numeqnD' +2
		return scalar zpos_end = `numeqnD' + `numeqnZ' + 1
	}
	else if ("`model'"=="fiv") {		
		return scalar zpos_start = `numeqnD' +2
		// return scalar zpos_end = `numeqnD' + `numeqnDH' + 1
		// not sure this will work
		return scalar zpos_end = 2*`numeqnD' + 1
	}
	else {
		return scalar zpos_start = 0
		return scalar zpos_end = 0
	}
	return scalar numD = `numeqnD'
	return scalar numZ = `numeqnZ'
	return local yvars `vtlistY'
	return local dvars `Dtilde'
	return local zvars `Ztilde' `DHtilde'

end
