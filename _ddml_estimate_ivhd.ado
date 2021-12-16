*** ddml estimation: partial linear IV model

program _ddml_estimate_ivhd, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								REP(integer 0)	/// resampling iteration for the estimation, 0=use all
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
	// not needed?
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
	local DHlist `r(zstr)' 

	// replist empty => do for first resample
	// replist = "all" do for all resamples
	mata: st_local("numreps",strofreal(`mname'.nreps))
	if "`replist'"=="" {
		local replist 1
	}
	else if "`replist'"=="all" {
		numlist "1/`numreps'"
		local replist "`r(numlist)'"
	}
	else {
		numlist "`replist'"
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
		// mata: `eqn' = (*(`mname'.peqnAA)).get("`nameY'")
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
			local Dopt `Dopt' `oneDopt'
			mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
			local DHopt `DHopt' `oneDHopt'
		}
		*** shortstack names
		if `ssflag' {
			local Yss `nameY'_ss
			foreach var of varlist `nameD' {
				local Dss `Dss' `var'_ss
				local DHss `DHss' `var'_ss_h
			}
		}

		// text used in output below
		if `numreps'>1 {
			local stext " (sample=`m')"
		}

		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
			   	tokenize `ylist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local y `s(vnames)'
				tokenize `Dlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local d `s(vnames)'
				tokenize `DHlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local dh `s(vnames)'
				_ddml_optiv, yvar(`nameY') ytilde(`y')				///
					dvar(`nameD') dhtilde(`dh') dtilde(`d')			///
					touse(`touse') `debug'
				di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
				di as text "E[Y|X]   = " as res "`y'"
				di as text "E[D|X]   = " as res "`d'"
				di as text "E[D^|X,Z] = " as res "`dh'"
				ereturn display
			}
		}

		// estimate best model
		// local nodisp 
		// if `nreplist'>1 local nodisp qui
		`nodisp' {
			_ddml_optiv, yvar(`nameY') ytilde(`Yopt'_`m')				///
				dvar(`nameD') dhtilde(`DHopt'_`m') dtilde(`Dopt'_`m')	///
				touse(`touse') `debug'
			if `ncombos' > 1 {
				di as text "Optimal DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[Y|X]   = " as res "`Yopt'"
			di as text "E[D|X]   = " as res "`Dopt'"
			di as text "E[D^|X,Z] = " as res "`DHopt'"
			ereturn display
		}
		
		if `ssflag' {
			_ddml_optiv, yvar(`nameY') ytilde(`Yss'_`m')				///
				dvar(`nameD') dhtilde(`DHss'_`m') dtilde(`Dss'_`m')		///
				touse(`touse') `debug'
			if `ncombos' > 1 {
				di as text "Shortstack DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[Y|X]   = " as res "`Yss'"
			di as text "E[D|X]   = " as res "`Dss'"
			di as text "E[D^|X,Z] = " as res "`DHss'"
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


// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	syntax anything , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end
