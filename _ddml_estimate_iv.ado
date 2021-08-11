*** ddml estimation: partial linear IV model

program _ddml_estimate_iv, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								replist(string)	/// list of resamplings to estimate
								avplot			///
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

	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
	mata: st_local("Ztilde",invtokens(`mname'.nameZtilde))
	mata: st_local("nameY",invtokens(`mname'.nameY))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))

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
	local Zlist `r(zstr)' 

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
				do_regress `y' `d' (`z') if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
				di
				di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
				di as text "E[y|X] = " as res "`y'"
				di as text "E[D|X] = " as res "`d'"
				di as text "E[Z|X] = " as res "`z'"
				ereturn di
				
			}
		}
	
		//mata: `mname'.nameDtilde
		mata: st_local("Yopt",`mname'.nameYopt[`m'])
		mata: st_local("Dopt",invtokens(`mname'.nameDopt[`m',.]))
		mata: st_local("Zopt",invtokens(`mname'.nameZopt[`m',.]))
		// do_regress is OLS/IV but with original varnames
		add_suffix `Yopt' `Dopt', suffix("_`m'")
		local YD `s(vnames)'
		add_suffix `Zopt', suffix("_`m'")
		local IVlist `s(vnames)'
		do_regress `YD' (`IVlist') if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')

		// display
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
		di as text "E[Z|X] = " as res "`Zopt'"
		ereturn display
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

/*
mata:

string scalar mat_to_varlist(string matrix inmat)
{
	r = rows(inmat)
	for (i=1;i<=r;i++) {

		if (i==1) {
			str = invtokens(inmat[i,]) 
		}
		else {
			str = str + " | " + invtokens(inmat[i,]) 
		}
	} 
	return(str)
}
end
*/
