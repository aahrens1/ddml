*** ddml estimation: partial linear IV model

program _ddml_estimate_iv, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								REP(integer 0)	/// resampling iteration for the estimation, 0=use all
								avplot			///
								debug			///
								* ]

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
	
	// ssflag is a model characteristic but is well-defined only if every equation has multiple learners.
	// will need to check this...
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
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
	
	// replist empty => do all
	// replist = integer; do for specified resample iteration
	mata: st_local("numreps",strofreal(`mname'.nreps))
	if `rep'==0 { // use all
		numlist "1/`numreps'"
		local replist "`r(numlist)'"
	}
	else {
		numlist "`rep'"
		local replist "`r(numlist)'"
	}
	local nreplist : word count `replist'

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
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("oneZopt",return_learner_item(`eqn',"opt","`m'"))
			local Zopt `Zopt' `oneZopt'
		}
		*** shortstack names
		if `ssflag' {
			local Yss `nameY'_ss
			foreach var of varlist `nameD' {
				local Dss `Dss' `var'_ss
			}
			foreach var of varlist `nameZ' {
				local Zss `Zss' `var'_ss
			}
		}

		// text used in output below
		if `numreps'>1 {
			local stext " (sample=`m')"
		}

		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
				tokenize `ylist' , parse("-")
				local y ``i''_`m'
				// remove extraneous space before the "_"
				local y : subinstr local y " " "", all
				tokenize `Dlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local d `s(vnames)'
				tokenize `Zlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local z `s(vnames)'
				// do_regress is OLS/IV but with original varnames
				if ("`y'"!="`Yopt'"|"`d'"!="`Dopt'"|"`z'"!="Zopt") { // omit if opt model to show at the end
					do_regress `y' `d' (`z') if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
					di
					di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
					di as text "E[y|X] = " as res "`y'"
					di as text "E[D|X] = " as res "`d'"
					di as text "E[Z|X] = " as res "`z'"
					ereturn display
				}
			}
		}

		// estimate best model
		// local nodisp 
		// if `nreplist'>1 local nodisp qui
		`nodisp' {
			add_suffix `Yopt' `Dopt', suffix("_`m'")
			local YD `s(vnames)'
			add_suffix `Zopt', suffix("_`m'")
			local IVlist `s(vnames)'
			do_regress `YD' (`IVlist') if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
			// display
			if `ncombos' > 1 {
				di as text "Optimal DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X] = " as res "`Yopt'"
			di as text "E[D|X] = " as res "`Dopt'"
			di as text "E[Z|X] = " as res "`Zopt'"
			ereturn display
		}
		
		if `ssflag' {
			add_suffix `Yss' `Dss', suffix("_`m'")
			local YD `s(vnames)'
			add_suffix `Zss', suffix("_`m'")
			local IVlist `s(vnames)'
			do_regress `YD' (`IVlist') if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
			di as text "Shortstack DML model`stext':" _c
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X] = " as res "`Yss'"
			di as text "E[D|X] = " as res "`Dss'"
			di as text "E[Z|X] = " as res "`Zss'"
			ereturn display
		}
		
		/*
		*** aggregate over resampling iterations if there is more than one
		if `nreplist'>1 {
			tempname bi vi
			mat `bi' = e(b)
			mat `vi' = e(V)
			if `m'==1 {
				mat `bagg' = 1/`numreps' * `bi'
				mat `vagg' = 1/`numreps' * `vi'
			} 
			else {
				mat `bagg' = `bagg' + 1/`numreps' * `bi'
				mat `vagg' = `vagg' + 1/`numreps' * `vi'				
			}
			if `m'==`numreps' {
				local N = e(N)
				ereturn clear
				ereturn post `bagg' `vagg', depname(`nameY') obs(`N') esample(`touse')
				ereturn display
			}
		}
		*/
	}

end


// does OLS/IV and reports with substitute yname and dnames
program define do_regress, eclass
	syntax anything [if] [in] , [ yname(name) dnames(namelist) * ]

	marksample touse

	qui reg `anything' if `touse', `options'

	tempname b
	tempname V
	mat `b' = e(b)
	mat `V' = e(V)
	matrix colnames `b' = `dnames'
	matrix rownames `b' = `yname'
 	matrix colnames `V' = `dnames'
	matrix rownames `V' = `dnames'
	local N = e(N)
	ereturn clear
	ereturn post `b' `V', depname(`yname') obs(`N') esample(`touse')

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
