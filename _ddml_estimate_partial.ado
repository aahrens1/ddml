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

	/*make_varlists, mname(`mname')

	local ylist		`s(ylist)'
	local Dlist		`s(Dlist)'
	local ncombos	= s(ncombos)
	local tokenlen	= `ncombos'*2*/

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
		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
				tokenize `ylist' , parse("-")
				local y ``i''
				tokenize `Dlist' , parse("-")
				local d ``i''
				//add_prefix `y' `d', prefix("")
				add_suffix `y' `d', suffix("_`m'")
				// do_regress is OLS but with original varnames
				// do_regress `y' `d' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
				do_regress `s(vnames)' if `touse' , nocons `robust' yname(`nameY') dnames(`nameD')
				di
				di as res "DML (sample=`m') with Y=`y' and D=`d' (N=`e(N)'):"
				ereturn di
		     }
		}
	
		*** estimate best model
	   	mata: st_local("Yopt",`mname'.nameYopt[`m'])
   		mata: st_local("Dopt",invtokens(`mname'.nameDopt[`m',.]))
		//add_prefix `Yopt' `Dopt', prefix("")
		// do_regress is OLS but with original varnames
		add_suffix `Yopt' `Dopt', suffix("_`m'")
		// do_regress `Yopt' `Dopt' if `touse' , nocons `robust' yname(`Yopt') dnames(`Dopt')
		do_regress `s(vnames)' if `touse' , nocons `robust' yname(`Yopt') dnames(`Dopt')
	
		// plot
		if ("`avplot'"!="") {
		   twoway (scatter `s(vnames)') (lfit `s(vnames)')
		}
	
		// display
		di
		if `ncombos' > 1 {
			di as res "Optimal model: DML (sample=`m') with optimal Y=`Yopt' and optimal D=`Dopt' (N=`e(N)'):"
		}
		else {
			di as res "DML (sample=`m') with Y=`Yopt' and D=`Dopt' (N=`e(N)'):"
		}
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

// does OLS and reports with substitute yname and dnames
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
