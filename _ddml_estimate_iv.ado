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


	/* not in use (?)
	_ddml_make_varlists, mname(`mname')
	if ("`debug'"!="") {
		return list
	}
	*/
	mata: st_local("Ztilde",invtokens(`mname'.nameZtilde))
	mata: st_local("Dtilde",invtokens(`mname'.nameDtilde))
	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameY",invtokens(`mname'.nameY))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))

	/*
	// seems to require _ddml_make_varlists?
	_ddml_allcombos `r(eq)' ,								///
		/* putlast(`Yopt' `Dopt' `Zopt') */					///
		`debug'												///
		dpos_end(`r(dpos_end)')								///
		zpos_start(`r(zpos_start)') zpos_end(`r(zpos_end)')	///
		addprefix("")
	*/
	_ddml_allcombos `Ytilde' - `Dtilde' - `Ztilde' ,	///
		`debug'											///
		addprefix("")

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2
	local ylist `r(colstr1)'
	local Dlist `r(colstr2)'
	local Zlist `r(colstr3)'

	/*
	local ylist `r(ystr)'
	local Dlist `r(dstr)'
	local Zlist `r(zstr)' 
	*/

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
				tokenize `ylist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local y `s(vnames)'
				tokenize `Dlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local d `s(vnames)'
				tokenize `Zlist' , parse("-")
				add_suffix ``i'', suffix("_`m'")
				local z `s(vnames)'
				di
				di as res "DML (sample = `m') with Y=`y', D=`d', Z=`z':"
			   	ivreg2 `y' (`d'=`z') if `touse', nocons `robust' noheader nofooter
			 }
		}
	
		//mata: `mname'.nameDtilde
		mata: st_local("Yopt",`mname'.nameYopt[`m'])
		mata: st_local("Dopt",invtokens(`mname'.nameDopt[`m']))
		mata: st_local("Zopt",invtokens(`mname'.nameZopt[`m']))
		di
		di as res "Optimal model: DML (sample = `m') with Y=`Yopt', D=`Dopt', Z=`Zopt':"
	   	ivreg2 `Yopt'_`m' (`Dopt'_`m'=`Zopt'_`m') if `touse', nocons `robust' noheader nofooter
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
