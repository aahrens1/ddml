*** ddml estimation: linear models
program _ddml_estimate_linear, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	///
								clear			/// deletes all tilde-variables (to be implemented)
								post(string)	/// specification to post/display
								REP(string)		/// resampling iteration to post/display or mean/median
								RESULTSonly		/// model has been estimated, just display results
								* ]
	
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	// model needs to be estimated
	local estflag = "`resultsonly'"==""
	
	// if post not specified, post optimal model
	if "`post'"=="" {
		local post "opt"
	}
	if "`post'"=="shortstack" & ~`ssflag' {
		di as err "error - no shortstack crossfit estimates available to post"
		exit 198
	}
	// if rep not specified, default is rep=1
	if "`rep'"=="" {
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
	if ~`crossfitted' {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}
	// opt + mean/median not currently supported
	if "`post'"=="opt" & real("`rep'")==. {
		di as err "error - post(opt) model not yet avaiable with mean/median"
		exit 198
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
		// recover optimal specs
		forvalues m=1/`nreps' {
			mata: st_local("spec",(`mname'.estAA).get(("optspec","`m'")))
			local optspec`m' = `spec'
		}
	}
	else {
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
		mata: `bmat' = J(`ncombos'*`nreps',`numeqnD',.)
		mata: `semat' = J(`ncombos'*`nreps',`numeqnD',.)
		
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
			
		*** shortstack names
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
				if "`model'"=="ivhd" {
					mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
					local Zopt `Zopt' `oneDHopt'
				}
			}
			// nameZ is empty for ivhd model
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
					local title Optimal `title'
					// save in AA
					mata: (`mname'.estAA).put(("optspec","`m'"),"`i'")
				}
				if "`model'"=="ivhd" {
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
				qui _ddml_reg if `mname'_sample_`m',							///
						nocons `robust'											///
						y(`y') yname(`nameY')									///
						d(`d') dnames(`nameD') dvtnames(`dvtnames')		 		///
						z(`z') znames(`nameZ') zvtnames(`zvtnames')				///
						mname(`mname') spec(`i') rep(`m') title(`title') `norep'
				
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
			}
			
			if `ssflag' {
				if "`model'"=="ivhd" {
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
				qui _ddml_reg if `mname'_sample_`m',							///
						nocons `robust'											///
						y(`Yss') yname(`nameY')									///
						d(`d') dnames(`nameD') dvtnames(`dvtnames') 			///
						z(`z') znames(`nameZ') zvtnames(`zvtnames')				///
						mname(`mname') spec(ss) rep(`m') title(`title') `norep'
			}
	
		}
		
		// aggregate across resamplings
		if `nreps' > 1 {
			// numbered specifications
			forvalues i = 1/`ncombos' {
				local title "DDML model, specification `i' (mean)"
				qui _ddml_reg, mname(`mname') spec(`i') medmean(mn) title(`title')
				local title "DDML model, specification `i' (median)"
				qui _ddml_reg, mname(`mname') spec(`i') medmean(md) title(`title')
			}
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean)"
				qui _ddml_reg, mname(`mname') spec(ss) medmean(mn) title(`title')
				local title "Shortstack DDML model (median)"
				qui _ddml_reg, mname(`mname') spec(ss) medmean(md) title(`title')
			}
		}
		
		// estimation complete
		mata: `mname'.ncombos = `ncombos'
		mata: (`mname'.estAA).put(("nmat","all"),`nmat')
		mata: (`mname'.estAA).put(("bmat","all"),`bmat')
		mata: (`mname'.estAA).put(("semat","all"),`semat')
	}
	
	************** REPORT RESULTS **************
	
	if "`show'"=="all" {
		forvalues m=1/`nreps' {
			// all combos including optimal model
			forvalues i=1/`ncombos' {
				di
				_ddml_reg, mname(`mname') spec(`i') rep(`m') replay
				di
			}
		}
	}
		
	if "`show'"=="optimal" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			di
			_ddml_reg, mname(`mname') spec(`optspec`m'') rep(`m') replay
			di
		}
	}
	
	if `ssflag' {
		if "`show'"=="shortstack" | "`show'"=="all" {
			forvalues m=1/`nreps' {
				di
				_ddml_reg, mname(`mname') spec(ss) rep(`m') replay
				di
			}
		}
	}
	
	*** Summary results ***
	
	di
	di as text "Summary DDML estimation results:"
	di as text "spec  r" %14s "Y learner" _c
	forvalues j=1/`numeqnD' {
		di as text %14s "D learner" %10s "b" %10s "SE" _c
	}
	if "`model'"=="ivhd" {
		forvalues j=1/`numeqnD' {
			di as text %14s "DH learner" _c
		}
	}
	forvalues j=1/`numeqnZ' {
		di as text %14s "Z learner" _c
	}
	di
	forvalues m=1/`nreps' {
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
			local rcmd stata _ddml_reg, mname(`mname') spec(`i') rep(`m') replay
			di %6s "{`rcmd':`specrep'}" _c
			di %14s "`yt'" _c
			forvalues j=1/`numeqnD' {
				local vt : word `j' of `dtlist'
				mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
				mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
				di %14s "`vt'" _c
				di %10.3f `b' _c
				local pse (`: di %6.3f `se'')
				di %10s "`pse'" _c
			}
			forvalues j=1/`numeqnZ' {
				local vt : word `j' of `ztlist'
				di %14s "`vt'" _c
			}
			if "`model'"=="ivhd" {
				forvalues j=1/`numeqnD' {
					local vt : word `j' of `ztlist'
					di %14s "`vt'" _c
				}
			}
			di
		}
		if `ssflag' {
			qui _ddml_reg, mname(`mname') spec(ss) rep(`m') replay
			local specrep `: di "ss" %3.0f `m''
			// pad out to 6 spaces
			local specrep = "  " + "`specrep'"
			local rcmd stata _ddml_reg, mname(`mname') spec(ss) rep(`m') replay
			di %6s "{`rcmd':`specrep'}" _c
			di %14s "[shortstack]" _c
			forvalues j=1/`numeqnD' {
				di %14s "[ss]" _c
				di %10.3f el(e(b),1,`j') _c
				local pse (`: di %6.3f sqrt(el(e(V),`j',`j'))')
				di %10s "`pse'" _c
			}
			if "`model'"=="ivhd" {
				forvalues j=1/`numeqnD' {
					di %14s "[ss]" _c
				}
			}
			forvalues j=1/`numeqnZ' {
				di %14s "[ss]" _c
			}
			di
		}
	}
	if `nreps' > 1 {
		di as text "Mean/median:"
		foreach medmean in mn md {
			forvalues i=1/`ncombos' {
				qui _ddml_reg, mname(`mname') spec(`i') rep(`medmean') replay
				mata: st_local("yt",abbrev(`nmat'[`i',1],13))
				mata: st_local("dtlist",invtokens(abbrev(tokens(`nmat'[`i',2]),13)))
				mata: st_local("ztlist",`nmat'[`i',3])
				if "`ztlist'"~="" {
					mata: st_local("ztlist",invtokens(abbrev(tokens(`nmat'[`i',3]),13)))
				}
				di " " _c
				local specrep `: di %3.0f `i' %3s "`medmean'"'
				// pad out to 6 spaces
				local specrep = (6-length("`specrep'"))*" " + "`specrep'"
				local rcmd stata _ddml_reg, mname(`mname') spec(`i') rep(`medmean') replay
				di %6s "{`rcmd':`specrep'}" _c
				di %14s "`yt'" _c
				forvalues j=1/`numeqnD' {
					local vt : word `j' of `dtlist'
					di %14s "`vt'" _c
					di %10.3f el(e(b),1,`j') _c
					local pse (`: di %6.3f sqrt(el(e(V),`j',`j'))')
					di %10s "`pse'" _c
				}
				forvalues j=1/`numeqnZ' {
					local vt : word `j' of `ztlist'
					di %14s "`vt'" _c
				}
				if "`model'"=="ivhd" {
					forvalues j=1/`numeqnD' {
						local vt : word `j' of `ztlist'
						di %14s "`vt'" _c
					}
				}
				di
			}
			if `ssflag' {
				qui _ddml_reg, mname(`mname') spec(ss) rep(`medmean') replay
				local specrep `: di "ss" %3s "`medmean'"'
				// pad out to 6 spaces
				local specrep = "  " + "`specrep'"
				local rcmd stata _ddml_reg, mname(`mname') spec(ss) rep(`medmean') replay
				di %6s "{`rcmd':`specrep'}" _c
				di %14s "[shortstack]" _c
				forvalues j=1/`numeqnD' {
					di %14s "[ss]" _c
					di %10.3f el(e(b),1,`j') _c
					local pse (`: di %6.3f sqrt(el(e(V),`j',`j'))')
					di %10s "`pse'" _c
				}
				if "`model'"=="ivhd" {
					forvalues j=1/`numeqnD' {
						di %14s "[ss]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
					di %14s "[ss]" _c
				}
				di
			}
		}
	}
	
	// post selected estimates; rep is the resample number (default=1)
	if "`post'"=="opt" {
		di
		_ddml_reg, mname(`mname') spec(`optspec`rep'') rep(`rep') replay
		di
	}
	else if "`post'"=="shortstack" {
		di
		_ddml_reg, mname(`mname') spec(ss) rep(`rep') replay
		di	
	}
	else {
		// post macro denotes the specification to post
		di
		_ddml_reg, mname(`mname') spec(`post') rep(`rep') replay
		di
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


