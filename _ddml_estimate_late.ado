*** ddml estimation: LATE model
* notes:
* add check somewhere that only a single D and a single Z are allowed.

program _ddml_estimate_late, eclass sortpreserve

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
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	
	// ssflag is a model characteristic but is well-defined only if every equation has multiple learners.
	// will need to check this...
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	// mata: `eqn' = (*(`mname'.peqnAA)).get(`mname'.nameY)
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	foreach var in `vtlistY' {
		local Y0tilde `Y0tilde' `var'0
		local Y1tilde `Y1tilde' `var'1
	}
	// mata: `eqn' = (*(`mname'.peqnAA)).get(`mname'.nameD)
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameD)
	mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
	foreach var in `vtlistD' {
		local D0tilde `D0tilde' `var'0
		local D1tilde `D1tilde' `var'1
	}
	// mata: `eqn' = (*(`mname'.peqnAA)).get(`mname'.nameZ)
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameZ)
	mata: st_local("Ztilde",invtokens(`eqn'.vtlist))

	if ("`debug'"!="") {
		di "`Ytilde'"
		di "`Ztilde'"
		di "`Dtilde'"
	}

	_ddml_allcombos `Y0tilde' - `Y1tilde' - `D0tilde' - `D1tilde' - `Ztilde' ,	///
		`debug'																	///
		addprefix("")

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2
	local y0list `r(colstr1)'
	local y1list `r(colstr2)'
	local d0list `r(colstr3)'
	local d1list `r(colstr4)'
	local Zlist `r(colstr5)' 

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
		// mata: `eqn' = (*(`mname'.peqnAA)).get("`nameY'")
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("Y0opt",return_learner_item(`eqn',"opt0","`m'"))
		mata: st_local("Y1opt",return_learner_item(`eqn',"opt1","`m'"))
		local Y0opt `Y0opt'0
		local Y1opt `Y1opt'1
		// mata: `eqn' = (*(`mname'.peqnAA)).get("`nameD'")
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		mata: st_local("D0opt",return_learner_item(`eqn',"opt0","`m'"))
		mata: st_local("D1opt",return_learner_item(`eqn',"opt1","`m'"))
		local D0opt `D0opt'0
		local D1opt `D1opt'1
		// mata: `eqn' = (*(`mname'.peqnAA)).get("`nameZ'")
		mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
		mata: st_local("Zopt",return_learner_item(`eqn',"opt","`m'"))
		
		*** shortstack names
		if `ssflag' {
			local Y0ss `nameY'_ss0
			local Y1ss `nameY'_ss1
			local D0ss `nameD'_ss0
			local D1ss `nameD'_ss1
			local Zss `nameZ'_ss
		}

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
				tokenize `d0list' , parse("-")
				local d0 ``i''
				tokenize `d1list' , parse("-")
				local d1 ``i''
				tokenize `Zlist' , parse("-")
				local z ``i''
				// omit if this is the optimal model, which will be estimated at the end
				if ("`Y0opt'"!="`y0'_`m'"|"`Y1opt'"!="`y1'_`m'"|"`D1opt'"!="`d1'_`m'"|"`D0opt'"!="`d0'_`m'"|"`Zopt'"!="`z'_`m'") {
					_ddml_late, yvar(`nameY') y0tilde(`y0'_`m') y1tilde(`y1'_`m')	///
							dvar(`nameD') d0tilde(`d0'_`m') d1tilde(`d1'_`m')	///
							zvar(`nameZ') ztilde(`z'_`m')						///
							touse(`touse')
					di as text "DML`stext':" _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
					di as text "E[y|X,Z=0] = " as res "`y0'_`m'"
					di as text "E[y|X,Z=1] = " as res "`y1'_`m'"
					di as text "E[D|X,Z=0] = " as res "`d0'_`m'"
					di as text "E[D|X,Z=1] = " as res "`d1'_`m'"
					di as text "E[Z|X]     = " as res "`z'_`m'"			
					ereturn display						
				}
			}
		}
	
		*** optimal model
		// local nodisp 
		// if `nreplist'>1 local nodisp qui
		`nodisp' {
			_ddml_late, yvar(`nameY') y0tilde(`Y0opt'_`m') y1tilde(`Y1opt'_`m')	///
						dvar(`nameD') d0tilde(`D0opt'_`m') d1tilde(`D1opt'_`m')	///
						zvar(`nameZ') ztilde(`Zopt'_`m')						///
						touse(`touse')  
			if `ncombos' > 1 {
				di as text "Optimal DML model`stext':" _c
			}
			else {
				di as text "DML`stext':" _c
			}
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X,Z=0] = " as res "`Y0opt'_`m'"
			di as text "E[y|X,Z=1] = " as res "`Y1opt'_`m'"
			di as text "E[D|X,Z=0] = " as res "`D0opt'_`m'"
			di as text "E[D|X,Z=1] = " as res "`D1opt'_`m'"
			di as text "E[Z|X]     = " as res "`Zopt'_`m'"
			ereturn display
		}
		
		if `ssflag' {
			_ddml_late, yvar(`nameY') y0tilde(`Y0ss'_`m') y1tilde(`Y1ss'_`m')	///
						dvar(`nameD') d0tilde(`D0ss'_`m') d1tilde(`D1ss'_`m')	///
						zvar(`nameZ') ztilde(`Zss'_`m')						///
						touse(`touse')  
			di as text "Shortstack DML model`stext':" _c
			di as text _col(52) "Number of obs   =" _col(70) as res %9.0f e(N)
			di as text "E[y|X,Z=0] = " as res "`Y0ss'_`m'"
			di as text "E[y|X,Z=1] = " as res "`Y1ss'_`m'"
			di as text "E[D|X,Z=0] = " as res "`D0ss'_`m'"
			di as text "E[D|X,Z=1] = " as res "`D1ss'_`m'"
			di as text "E[Z|X]     = " as res "`Zss'_`m'"
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
 