*** ddml estimation: interactive model
* notes:
* add check somewhere that only a single D is allowed.

program _ddml_estimate_interactive, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								REP(integer 0)	/// resampling iteration for the estimation, 0=use all
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
	mata: st_local("Y0tilde",invtokens(`mname'.nameY0tilde))
	mata: st_local("Y1tilde",invtokens(`mname'.nameY1tilde))
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameY",invtokens(`mname'.nameY))

	_ddml_allcombos `Y0tilde'- `Y1tilde' - `Dtilde',		///
		addprefix("`mname'_")										///
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

		*** retrieve best model
		mata: st_local("Y0opt",`mname'.nameY0opt[`m'])
		mata: st_local("Y1opt",`mname'.nameY1opt[`m'])
		mata: st_local("Dopt",`mname'.nameDopt[`m'])

		// text used in output below
		if `numreps'>1 {
			local stext " (sample=`m')"
		}

		if "`show'"=="all" {
			forvalues i = 1(2)`tokenlen' {
				tokenize `y0list' , parse("-")
				local y0 ``i''
				tokenize `y1list' , parse("-")
				local y1 ``i''
				tokenize `Dlist' , parse("-")
				local d ``i''
				_ddml_ate,				///
					yvar(`nameY')		///
					dvar(`nameD')		///
					y0tilde(`y0'_`m')	///
					y1tilde(`y1'_`m')	///
					dtilde(`d'_`m')		///
					touse(`touse')
				di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
				di as text "E[y|X,D=0] = " as res "`y0'_`m'"
				di as text "E[y|X,D=1] = " as res "`y1'_`m'"
				di as text "E[D|X]     = " as res "`d'_`m'"
				ereturn display
			}
		}

		// estimate best model
		local nodisp 
		if `nreplist'>1 local nodisp qui
		`nodisp' {
			_ddml_ate,					///
				 yvar(`nameY')			///
				 dvar(`nameD')			///
				 y0tilde(`Y0opt'_`m')	///
				 y1tilde(`Y1opt'_`m')	///
				 dtilde(`Dopt'_`m')		///
				 touse(`touse')
			if `ncombos' > 1 {
				di as text "Optimal DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X,D=0] = " as res "`Y0opt'_`m'"
			di as text "E[y|X,D=1] = " as res "`Y1opt'_`m'"
			di as text "E[D|X]     = " as res "`Dopt'_`m'"
			ereturn display
		}

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
			if `m'==`numreps' { // save & display on last iteration
				local N = e(N)
				ereturn clear
				di as text "Aggregate DML model over `nreplist' repetitions:"
				ereturn post `bagg' `vagg', depname(`nameY') obs(`N') esample(`touse')
				ereturn display
			}
		}
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
