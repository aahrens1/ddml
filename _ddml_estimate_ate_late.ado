*** ddml estimation: interactive model
* notes:
* add check somewhere that only a single D is allowed.

program _ddml_estimate_ate_late, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								post(string)	/// specification to post/display
								REP(integer 1)	/// resampling iteration to post/display
								avplot			///
								debug			///
								* ]
	
	// if post not specified, post optimal model
	if "`post'"=="" {
		local post "opt"
	}

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	
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
	if (`crossfitted'==0) {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}
	
	// ssflag is a model characteristic but is well-defined only if every equation has multiple learners.
	// will need to check this...
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
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
				local title Optimal `title'
			}
			// code works for both ATE and LATE
			qui _ddml_ate_late if `touse',	///
				yvar(`nameY')				///
				dvar(`nameD')				///
				zvar(`nameZ')				///
				y0tilde(`y0')				///
				y1tilde(`y1')				///
				dtilde(`d')					///
				d0tilde(`d0')				///
				d1tilde(`d1')				///
				ztilde(`z')					///
				spec(`i')					///
				rep(`m')					///
				mname(`mname')				///
				title(`title')				///
				touse(`touse')
			
			mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(b)")
			mata: `semat'[(`m'-1)*`ncombos'+`i',.] = sqrt(diagonal(st_matrix("e(V)"))')
			
		}
		
		if `ssflag' {
			
			// code works for both ATE and LATE
			local title "Shortstack DDML model`stext'"
			qui _ddml_ate_late if `touse',	///
				yvar(`nameY')				///
				dvar(`nameD')				///
				zvar(`nameZ')				///
				y0tilde(`Y0ss')				///
				y1tilde(`Y1ss')				///
				dtilde(`Dss')				///
				d0tilde(`D0ss')				///
				d1tilde(`D1ss')				///
				ztilde(`Zss')				///
				spec(ss)					///
				rep(`m')					///
				mname(`mname')				///
				title(`title')				///
				touse(`touse')
		
		}
	}
	
	// aggregate across resamplings
	if `nreps' > 1 {
		// numbered specifications
		forvalues i = 1/`ncombos' {
			local title "DDML model, specification `i' (mean)"
			qui _ddml_ate_late, mname(`mname') spec(`i') medmean(mean) title(`title')
			local title "DDML model, specification `i' (median)"
			qui _ddml_ate_late, mname(`mname') spec(`i') medmean(median) title(`title')
		}
		// shortstack
		local title "Shortstack DDML model (mean)"
		qui _ddml_ate_late, mname(`mname') spec(ss) medmean(mean) title(`title')
		local title "Shortstack DDML model (median)"
		qui _ddml_ate_late, mname(`mname') spec(ss) medmean(median) title(`title')
	}

	if "`show'"=="all" {
		forvalues m=1/`nreps' {
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			// all combos including optimal model
			forvalues i=1/`ncombos' {
				di
				_ddml_ate_late, mname(`mname') spec(`i') rep(`m') replay
				di
			}
			// shortstack
			di
			_ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
			di
		}
	}
		
	if "`show'"=="optimal" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			di
			_ddml_ate_late, mname(`mname') spec(`optspec`m'') rep(`m') replay
			di
		}
	}
		
	if "`show'"=="shortstack" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			di
			_ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
			di
		}
	}

	di
	di as text "Summary DDML estimation results:"
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
			mata: st_local("yt0",`nmat'[`i',1])
			mata: st_local("yt1",`nmat'[`i',2])
			if "`optspec`m''"=="`i'" {
				di "*" _c
			}
			else {
				di " " _c
			}
			local specrep `: di %3.0f `i' %3.0f `m''
			// pad out to 6 spaces
			local specrep = (6-length("`specrep'"))*" " + "`specrep'"
			local rcmd stata _ddml_ate_late, mname(`mname') spec(`i') rep(`m') replay
			di %6s "{`rcmd':`specrep'}" _c
			di %14s "`yt0'" _c
			di %14s "`yt1'" _c
			if `ateflag' {
				mata: st_local("dt",`nmat'[`i',3])
				di %14s "`dt'" _c
			}
			else {
				mata: st_local("dt0",`nmat'[`i',3])
				mata: st_local("dt1",`nmat'[`i',4])
				di %14s "`dt0'" _c
				di %14s "`dt1'" _c
			}
			mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
			mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
			di %10.3f `b' _c
			local pse (`: di %6.3f `se'')
			di %10s "`pse'" _c
			if ~`ateflag' {
				mata: st_local("zt",`nmat'[`i',5])
				di %14s "`zt'" _c
			}
			di
		}

		if `ssflag' {
			qui _ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
			local specrep `: di "ss" %3.0f `m''
			// pad out to 6 spaces
			local specrep = "  " + "`specrep'"
			local rcmd stata _ddml_ate_late, mname(`mname') spec(ss) rep(`m') replay
			di %6s "{`rcmd':`specrep'}" _c
			di %14s "[shortstack]" _c
			di %14s "[ss]" _c
			di %14s "[ss]" _c
			if ~`ateflag' {
				di %14s "[ss]" _c
			}
			di %10.3f el(e(b),1,1) _c
			local pse (`: di %6.3f sqrt(el(e(V),1,1))')
			di %10s "`pse'" _c
			if ~`ateflag' {
				di %14s "[ss]" _c
			}
			di
		}

	}
	if `nreps' > 1 {
		di as text "Mean/median:"
		foreach medmean in mean median {
			if "`medmean'"=="mean" {
				local mm mn
			}
			else {
				local mm md
			}
			forvalues i=1/`ncombos' {
				qui _ddml_ate_late, mname(`mname') spec(`i') rep(`medmean') replay
				mata: st_local("yt",abbrev(`nmat'[`i',1],13))
				mata: st_local("dtlist",invtokens(abbrev(tokens(`nmat'[`i',2]),13)))
				mata: st_local("ztlist",`nmat'[`i',3])
				if "`ztlist'"~="" {
					mata: st_local("ztlist",invtokens(abbrev(tokens(`nmat'[`i',3]),13)))
				}
				di " " _c
				local specrep `: di %3.0f `i' %3s "`mm'"'
				// pad out to 6 spaces
				local specrep = (6-length("`specrep'"))*" " + "`specrep'"
				local rcmd stata _ddml_ate_late, mname(`mname') spec(`i') rep(`medmean') replay
				di %6s "{`rcmd':`specrep'}" _c
				di %14s "`yt'" _c
				forvalues j=1/`numeqnD' {
					local vt : word `j' of `dtlist'
					di %14s "`vt'" _c
					di %10.3f el(e(b),1,`j') _c
					local pse (`: di %6.3f sqrt(el(e(V),`j',`j'))')
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
				qui _ddml_ate_late, mname(`mname') spec(ss) rep(`medmean') replay
				local specrep `: di "ss" %3s "`mm'"'
				// pad out to 6 spaces
				local specrep = "  " + "`specrep'"
				local rcmd stata _ddml_ate_late, mname(`mname') spec(ss) rep(`medmean') replay
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
		_ddml_ate_late, mname(`mname') spec(`optspec`rep'') rep(`rep') replay
		di
	}
	else if "`post'"=="shortstack" {
		di
		_ddml_ate_late, mname(`mname') spec(ss) rep(`rep') replay
		di	
	}
	else {
		// post macro denotes the specification to post
		di
		_ddml_ate_late, mname(`mname') spec(`post') rep(`rep') replay
		di
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

