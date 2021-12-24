*** ddml estimation: partial linear model
program _ddml_estimate_partial2, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	///
								clear			/// deletes all tilde-variables (to be implemented)
								post(string)	/// specification to post/display
								REP(integer 1)	/// resampling iteration to post/display
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

	// locals used below
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	local numeqnD	: word count `nameD'
	mata: st_local("nreps",strofreal(`mname'.nreps))

	if (`crossfitted'==0) {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}
	
	// ssflag is a model characteristic but is well-defined only if every equation has multiple learners.
	// will need to check this...
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	// get varlists
	_ddml_make_varlists, mname(`mname')
	local yvars `r(yvars)'
	local dvars `r(dvars)'
	
	// obtain all combinations
	_ddml_allcombos `yvars' - `dvars' ,				///
		`debug'										///
		addprefix("")
	
	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2
	local ylist `r(ystr)'
	local Dlist `r(dstr)'
	
	tempname nmat bmat semat
	mata: `nmat' = J(`ncombos',2,"")
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
	
	*** shortstack names
	if `ssflag' {
		local Yss `nameY'_ss
		foreach var of varlist `nameD' {
			local Dss `Dss' `var'_ss
		}
	}
	
	forvalues m=1/`nreps' {
		
		// reset locals
		local Dopt
		
		*** retrieve best model
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
			local Dopt `Dopt' `oneDopt'
		}

		// text used in output below
		if `nreps'>1 {
			local stext " (sample=`m')"
		}
		
		forvalues i = 1/`ncombos' {
			mata: st_local("y",`nmat'[`i',1])
			mata: st_local("d",`nmat'[`i',2])
			// check if opt for this resample; note === for D so order doesn't matter
			local isopt
			local isYopt : list Yopt == y
			local isDopt : list Dopt === d
			if `isYopt' & `isDopt' {
				local optspec`m' = `i'
				local isopt *
				local DMLtitle "Optimal DDML model"
			}
			else {
				local DMLtitle "DDML"
			}
			// add resample suffixes and estimate
			local y `y'_`m'
			add_suffix `d', suffix("_`m'")
			local d `s(vnames)'
			qui _ddml_reg if `touse' , nocons `robust' y(`y') d(`d') yname(`nameY') dnames(`nameD')
			mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(b)")
			mata: `semat'[(`m'-1)*`ncombos'+`i',.] = sqrt(diagonal(st_matrix("e(V)"))')
			estimates store `mname'_`i'_`m', title("Model `mname', specification `i' resample `m'`isopt'")
		}
		
		if `ssflag' {
			local y `Yss'_`m'
			add_suffix `Dss', suffix("_`m'")
			local d `s(vnames)'
			qui _ddml_reg if `touse' , nocons `robust' y(`y') d(`d') yname(`nameY') dnames(`nameD')
			estimates store `mname'_ss_`m', title("Model `mname', shorstack resample `m'")
		}

	}
	
	if "`show'"=="all" {
		forvalues m=1/`nreps' {
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			// all combos including optimal model
			forvalues i=1/`ncombos' {
				qui estimates restore `mname'_`i'_`m'
				if "`optspec`m''"=="`i'" {
					local DMLtitle "Optimal DDML model"
				}
				else {
					local DMLtitle "DDML"
				}
				di
				di as text "`DMLtitle'`stext':"
				_ddml_reg	// uses replay
				di
			}
			// shortstack
			qui estimates restore `mname'_ss_`m'
			di
			di as text "Shortstack DDML model`stext':"
			_ddml_reg	// uses replay
			di
		}
	}
		
	if "`show'"=="optimal" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			qui estimates restore `mname'_`optspec`m''_`m'
			di
			di as text "Optimal DDML specification, resample `m': `optspec`m''"
			_ddml_reg
			di
		}
	}
		
	if "`show'"=="shortstack" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			qui estimates restore `mname'_ss_`m'
			di
			di as text "Shortstack DDML model, resample `m':"
			_ddml_reg
			di
		}
	}
	
	di
	di as text "Summary DDML estimation results:"
	di as text "spec  r" %14s "Y learner" _c
	forvalues j=1/`numeqnD' {
		di as text %14s "D learner" %10s "b" %10s "SE" _c
	}
	di
	forvalues m=1/`nreps' {
		forvalues i=1/`ncombos' {
			mata: st_local("yt",`nmat'[`i',1])
			mata: st_local("dtlist",`nmat'[`i',2])
			if "`optspec`m''"=="`i'" {
				di "*" _c
			}
			else {
				di " " _c
			}
			local specrep `: di %3.0f `i' %3.0f `m''
			// pad out to 6 spaces
			local specrep = (6-length("`specrep'"))*" " + "`specrep'"
			di %6s "{stata estimates replay `mname'_`i'_`m':`specrep'}" _c
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
			di
		}
		if `ssflag' {
			qui estimates restore `mname'_ss_`m'
			local specrep `: di "ss" %3.0f `m''
			// pad out to 6 spaces
			local specrep = "  " + "`specrep'"
			di %6s "{stata estimates replay `mname'_ss_`m':`specrep'}" _c			
			di %14s "[shortstack]" _c
			forvalues j=1/`numeqnD' {
				di %14s "[ss]" _c
				di %10.3f el(e(b),1,`j') _c
				local pse (`: di %6.3f el(e(V),`j',`j')')
				di %10s "`pse'" _c
			}
			di
		}
	}
	
	// post selected estimates; rep is the resample number (default=1)
	if "`post'"=="opt" {
		qui estimates restore `mname'_`optspec`rep''_`rep'
		di
		di as text "Optimal DDML specification, resample `rep': `optspec`rep''"
		_ddml_reg
		di
	}
	else if "`post'"=="shortstack" {
		qui estimates restore `mname'_ss_`rep'
		di
		di as text "Shortstack DDML model, resample `rep':"
		_ddml_reg
		di	
	}
	else {
		// post macro should be an integer denoting the specification
		qui estimates restore `mname'_`post'_`rep'
		di
		di as text "DDML specification `post', resample `rep':"
		_ddml_reg
		di
	}
	
	// temp Mata object no longer needed
	cap mata: mata drop `nmat' `bmat' `semat'

/*
	// do for each specified resamples
	tempname bagg vagg 
	foreach m in `replist' {
		
		// reset locals
		local Dopt
		local Dss
		
		*** retrieve best model
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
			local Dopt `Dopt' `oneDopt'
		}
		*** shortstack names
		if `ssflag' {
			local Yss `nameY'_ss
			foreach var of varlist `nameD' {
				local Dss `Dss' `var'_ss
			}
		}

		// text used in output below
		if `nreps'>1 {
			local stext " (sample=`m')"
		}
		
		// should show have other option values?
		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
				tokenize `ylist' , parse("-")
				local y ``i''_`m'
				// remove extraneous space before the "_"
				local y : subinstr local y " " "", all
				tokenize `Dlist' , parse("-")
				local d ``i''
				if ("`y'"!="`Yopt'"|"`d'"!="Dopt") { // omit if opt model to show at the end
					add_suffix `d', suffix("_`m'")
					local d `s(vnames)'
					// do_regress is OLS/IV but with original varnames
					do_regress `y' `d' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
					di
					di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
					di as text "E[y|X] = " as res "`y'"
					di as text "E[D|X] = " as res "`d'"
					ereturn display				
				}
		    }
		}

		// estimate best model
		// local nodisp 
		// if `nreplist'>1 local nodisp qui
		`nodisp' {
			add_suffix `Yopt' `Dopt', suffix("_`m'")
			do_regress `s(vnames)' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
			di
			if `ncombos' > 1 {
				di as text "Optimal DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X] = " as res "`Yopt'"
			di as text "E[D|X] = " as res "`Dopt'"
			ereturn display
		}
		
		if `ssflag' {
			add_suffix `Yss' `Dss', suffix("_`m'")
			do_regress `s(vnames)' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
			di
			di as text "Shortstack DML model`stext':" _c
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X] = " as res "`Yss'"
			di as text "E[D|X] = " as res "`Dss'"
			ereturn display
		}
		
		/*
		*** aggregate over resampling iterations if there is more than one
		if `nreplist'>1 {
			tempname bi vi
			mat `bi' = e(b)
			mat `vi' = e(V)
			if `m'==1 {
				mat `bagg' = 1/`nreps' * `bi'
				mat `vagg' = 1/`nreps' * `vi'
			} 
			else {
				mat `bagg' = `bagg' + 1/`nreps' * `bi'
				mat `vagg' = `vagg' + 1/`nreps' * `vi'				
			}
			if `m'==`nreps' { // save & display on last iteration
				local N = e(N)
				ereturn clear
				di as text "Aggregate DML model over `nreplist' repetitions:"
				ereturn post `bagg' `vagg', depname(`nameY') obs(`N') esample(`touse')
				ereturn display
			}
		}
		*/
	}

	_ddml_ereturn, mname(`mname')
		
	// plot
	//if ("`avplot'"!="") {
	//   twoway (scatter `s(vnames)') (lfit `s(vnames)')
	//}
*/

end

// adds model name prefixes to list of varnames
program define add_prefix, sclass
	syntax anything , prefix(name)

	// anything is a list of to-be-varnames that need prefix added to them
	foreach vn in `anything' {
		local vnames `vnames' `prefix'`vn' 
	}
	
	sreturn local vnames `vnames'
end

// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	syntax anything , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end


