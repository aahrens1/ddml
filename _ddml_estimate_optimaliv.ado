*** ddml estimation: partial linear IV model

program _ddml_estimate_optimaliv, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								replist(string)	/// list of resamplings to estimate
								debug			///
								* ]

	// what does this do?
	if ("`show'"=="") {
		local show opt 
	}

	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

	//mata: `mname'.nameDtilde
	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
	mata: st_local("DHtilde",invtokens(`mname'.nameDHtilde))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameY",invtokens(`mname'.nameY))

	// get varlists
	_ddml_make_varlists, mname(`mname')
	// obtain all combinations
	_ddml_allcombos `r(yvars)' - `r(dvars)' - `r(zvars)' ,	///
		`debug'												///
		zpos(`r(zpos_start)') 								///
		addprefix("")

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2
	local ylist `r(ystr)'
	local Dlist `r(dstr)'
	local DHlist `r(zstr)' 

	// replist empty => do for first resample
	// replist = "all" do for all resamples
	mata: st_local("numreps",strofreal(cols(`mname'.idFold)))
	// subtract 1 (first col is id variable)
	local numreps = `numreps' - 1
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

	// do for each specified resamples
	foreach m in `replist' {
		// text used in output below
		if `numreps'>1 {
			local stext " (sample=`m')"
		}
		forvalues i = 1(2)`tokenlen' {
			if "`show'"=="all" {
			   	tokenize `ylist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local y `s(vnames)'
				tokenize `Dlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local d `s(vnames)'
				tokenize `DHlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local dh `s(vnames)'
				di
				di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
				di as text "E[Y|X]   = " as res "`y'"
				di as text "E[D|X]   = " as res "`d'"
				di as text "E[D|X,Z] = " as res "`dh'"
				_ddml_optiv, yvar(`nameY') ytilde(`y')				///
					dvar(`nameD') dhtilde(`dh') dtilde(`d')			///
					touse(`touse') `debug'
			}
		}
	
		mata: st_local("Yopt",`mname'.nameYopt[`m'])
		add_suffix `Yopt', suffix("_`m'")
		local Yopt `s(vnames)'
		mata: st_local("Dopt",invtokens(`mname'.nameDopt[`m',.]))
		add_suffix `Dopt', suffix("_`m'")
		local Dopt `s(vnames)'
		mata: st_local("DHopt",invtokens(`mname'.nameDHopt[`m',.]))
		add_suffix `DHopt', suffix("_`m'")
		local DHopt `s(vnames)'
		
		di
		if `ncombos' > 1 {
			di as text "Optimal DML model`stext':" _c
		}
		else {
			di as text "DML`stext':" _c
		}
		di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
		di as text "E[Y|X]   = " as res "`Yopt'"
		di as text "E[D|X]   = " as res "`Dopt'"
		di as text "E[D|X,Z] = " as res "`DHopt'"
		_ddml_optiv, yvar(`nameY') ytilde(`Yopt')				///
			dvar(`nameD') dhtilde(`DHopt') dtilde(`Dopt')		///
			touse(`touse') `debug'

	}

	/*
	// display
	tempname b
	tempname V 
	mat `b' = e(b)
	mat `V' = e(V)
	matrix colnames `b' = `nameD'
	matrix rownames `b' = `nameY'
 	matrix colnames `V' = `nameD'
	matrix rownames `V' = `nameD'
	local N = e(N)
	ereturn clear
	ereturn post `b' `V', depname(`Yopt') obs(`N') esample(`touse')
	if "`robust'"~="" {
		ereturn local vcetype	robust
	}
	ereturn display
	*/
	
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
