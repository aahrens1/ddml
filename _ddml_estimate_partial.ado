*** ddml estimation: partial linear model
program _ddml_estimate_partial, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								replist(string)	/// list of resamplings to estimate
								* ]
	
	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

	// locals used below
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))

	// get varlists
	_ddml_make_varlists, mname(`mname')
	// obtain all combinations
	_ddml_allcombos `r(yvars)' - `r(dvars)' ,				///
		`debug'												///
		addprefix("")

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2
	local ylist `r(ystr)'
	local Dlist `r(dstr)'

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
		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
				tokenize `ylist' , parse("-")
				local y ``i''_`m'
				// remove extraneous space before the "_"
				local y : subinstr local y " " "", all
				tokenize `Dlist' , parse("-")
				local d ``i''
				add_suffix `d', suffix("_`m'")
				local d `s(vnames)'
				// do_regress is OLS/IV but with original varnames
				do_regress `y' `d' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
				di
				di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
				di as text "E[y|X] = " as res "`y'"
				di as text "E[D|X] = " as res "`d'"
				ereturn di
		     }
		}
	
		*** estimate best model
	   	mata: st_local("Yopt",`mname'.nameYopt[`m'])
   		mata: st_local("Dopt",invtokens(`mname'.nameDopt[`m',.]))
		// do_regress is OLS/IV but with original varnames
		add_suffix `Yopt' `Dopt', suffix("_`m'")
		do_regress `s(vnames)' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
	
		// plot
		if ("`avplot'"!="") {
		   twoway (scatter `s(vnames)') (lfit `s(vnames)')
		}
	
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
		ereturn display
	}
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

mata:

// returns an empty eqnStruct
struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
