*** ddml estimation: interactive model

program _ddml_estimate_interactive, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								replist(string)	/// list of resamplings to estimate
								avplot			///
								debug			///
								* ]

	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

	// what does this do?
	if ("`show'"=="") {
		local show opt 
	}
	
	//mata: `mname'.nameDtilde
	//mata: st_local("Ztilde",invtokens(`mname'.nameZtilde))
	mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameY",invtokens(`mname'.nameY))

	_ddml_allcombos `Ytilde'- `Ytilde' - `Dtilde',		///
		/* putlast(`Y0opt' `Y1opt' `Dopt') */			///
		addprefix("")									///
		`debug'  

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2 -1
	local y0list `r(colstr1)'
	local y1list `r(colstr2)'
	local Dlist `r(colstr3)'

	if "`debug'"~="" {
		di
		di "y0list: `y0list'"
		di "y1list: `y1list'"
		di "Dist:   `Dlist'"
		di
	}
	
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

		forvalues i = 1(2)`tokenlen' {
			// if ("`show'"=="all"|`i'==`tokenlen') {
			if "`show'"=="all" {
				tokenize `y0list' , parse("-")
				local y0 ``i''
				tokenize `y1list' , parse("-")
				local y1 ``i''
				tokenize `Dlist' , parse("-")
				local d ``i''
				di
				di as res "DML (sample = `m') with Y0=`y0', Y1=`y1', D=`d':"
				_ddml_ate,				///
					yvar(`nameY')		///
					dvar(`nameD')		///
					y0tilde(`y0'_`m')	///
					y1tilde(`y1'_`m')	///
					dtilde(`d'_`m')		///
					touse(`touse')
			}
		}

		mata: st_local("Y0opt",`mname'.nameY0opt[`m'])
		mata: st_local("Y1opt",`mname'.nameY1opt[`m'])
		mata: st_local("Dopt",`mname'.nameDopt[`m'])
		//mata: st_local("Zopt",`mname'.nameZopt[`m'])
		di
		di as res "Optimal model: DML (sample = `m') with Y0=`Y0opt', Y1=`Y1opt', D=`Dopt':"
		_ddml_ate,					///
			 yvar(`nameY')			///
			 dvar(`nameD')			///
			 y0tilde(`Y0opt'_`m')	///
			 y1tilde(`Y1opt'_`m')	///
			 dtilde(`Dopt'_`m')		///
			 touse(`touse')
	
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
			ereturn local vcetype   robust
		}
		ereturn display
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
