*** ddml cross-fitting

* notes:
* why is vtype field needed?
* crossfitting code in separate program
* reporting code in separate subroutine below
* debug reporting in separate subroutine below
* number of resamples set here in reps(.) option
* in eqn struct, (*p).idVtilde is a problem ... but is it needed? currently commented out.
* add noisily option

program _ddml_crossfit, eclass sortpreserve

	syntax [anything] ,								/// 
							[						///
							NOIsily					///
							debug					/// 
							Robust					///
							mname(name)				///
							]

	// no checks included yet

	local debugflag		= "`debug'"~=""
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()

	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens((`mname'.nameD)))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("reps",strofreal(`mname'.nreps))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	local touse `mname'_sample
	di as text "Model: `model'"
	
	// fold IDs
	forvalues m=1/`reps' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	di as text "Fold IDs: `fidlist'"
	
	// equations and learners
	mata: `eqn' = (*(`mname'.peqnAA)).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	di as text "Y eqn learners (`numlnrY'): `vtlistY'"
	local vtlist `vtlistY'
	if `numeqnD' {
		di as text "D equations (`numeqnD'): `nameD'"
		foreach var of varlist `nameD' {
			di as text _col(5) "D equation `var':"
			mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			di as text _col(10) "learners: `vtlistD'"
			local vtlist `vtlist' `vtlistD'
		}
	}
	if `numeqnZ' {
		di as text "Z equations (`numeqnZ'): `nameZ'"
		foreach var of varlist `nameZ' {
			di as text _col(5) "Z equation `var':"
			mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			di as text _col(10) "learners: `vtlistZ'"
			local vtlist `vtlist' `vtlistZ'
		}
	}
	di as text "All learners: `vtlist'"
	
	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname'
	}
	else {
		// Y equation, all learners
		mata: `eqn' = (*(`mname'.peqnAA)).get("`nameY'")
		if "`model'"=="interactive" {
			local treatvar	`nameD'
		}
		else if ("`model'"=="late") {
			local treatvar	`nameZ'
		}
		else {
			// clear local
			local treatvar
		}
		mata: st_local("ssname",`eqn'.shortstack)
		di as text "Cross-fitting Y equation: `nameY'"
		crossfit if `touse',						///
			ename(`eqn') noreplace					///
			foldvar(`fidlist')						///
			shortstack(`ssname')					///
			treatvar(`treatvar')					///
			resid `noisily'
		// D equations
		if `numeqnD' {
			if ("`model'"=="late") {
				local treatvar	`nameZ'
			}
			else {
				// clear local
				local treatvar
			}
			foreach var of varlist `nameD' {
				mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
				mata: st_local("ssname",`eqn'.shortstack)
				di as text "Cross-fitting D equation: `var'"
				// All learners for each D eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					treatvar(`treatvar')					///
					resid `noisily'
			}
		}
		// Z equations
		if `numeqnZ' {
			foreach var of varlist `nameZ' {
				mata: `eqn' = (*(`mname'.peqnAA)).get("`var'")
				mata: st_local("ssname",`eqn'.shortstack)
				di as text "Cross-fitting Z equation: `var'"
				// All learners for each Z eqn
				crossfit if `touse',						///
					ename(`eqn') noreplace					///
					foldvar(`fidlist')						///
					shortstack(`ssname')					///
					resid `noisily'
			}
		}
	}
	
	// set flag on model struct
	mata: `mname'.crossfitted = 1
	
	// report results by equation type with resamplings grouped together
	di
	di as res "Reporting crossfitting results:"
	_ddml_describe `mname'
	
	/*		
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsDH",strofreal(cols(`mname'.nameDHtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	di "Model: `model'"
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	di "Number of Y estimating equations: `numeqnsY'"
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
		di "Number of D estimating equations: `numeqnsD'"
	}
	if `numeqnsDH' {
		mata: st_local("listDH",invtokens(`mname'.nameDH))
		mata: st_local("listDHtilde",invtokens(`mname'.nameDHtilde))
		di "Number of DH estimating equations: `numeqnsDH'"
	}
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		mata: st_local("listZtilde",invtokens(`mname'.nameZtilde))
		di "Number of Z estimating equations: `numeqnsZ'"
	}
	// initialize opt lists with void matrices of correct dimension
	mata: st_local("numvarsD",strofreal(cols(`mname'.nameD)))
	mata: st_local("numvarsZ",strofreal(cols(`mname'.nameZ)))
	mata: `mname'.nameYopt = J(0,1,"")
	mata: `mname'.nameDopt = J(0,`numvarsD',"")
	mata: `mname'.nameDHopt = J(0,`numvarsD',"")
	mata: `mname'.nameZopt = J(0,`numvarsZ',"")
	// folds and resamplings
	mata: st_local("reps",strofreal(`mname'.nreps))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()


	if `debugflag' {
		*** report estimates for full sample for debugging purposes
		report_debugging `mname'
	}
	else {
		forvalues m=1/`reps' {
		
			di
			di as text "Starting cross-fitting (sample = `m')"

			*** initialize tilde variables
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("cvdone",strofreal(`eqn'.crossfitted))
				if (`cvdone'!=1) {
					/// macro m is resampling counter
					cap drop `vtilde'_`m'
					qui gen double `vtilde'_`m'=.
				}
			}

			*** do cross-fitting
			*** loop over equations
			forvalues i=1/`numeqns' {
	
				// initialize prior to calling crossfit
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("eststring",`eqn'.eststring)
				mata: st_local("eqntype",`eqn'.eqntype)

				if ("`model'"=="optimaliv") {
					mata: st_local("vtilde_h",`eqn'.Vtilde_h)
					mata: st_local("eststring_h",`eqn'.eststring_h)					
				}

				// seems to be unused
				// mata: st_local("vtype",`eqn'.vtype)
				local touse `mname'_sample

				// request residuals unless optimal IV model & deq or dheq
				if ~("`model'"=="optimaliv"&("`eqntype'"=="deq"|"`eqntype'"=="dheq")) {
					local resid 
					mata: `eqn'.resid = 0
				}
				else {
					local resid resid
					mata: `eqn'.resid = 1
				}
				
				di as text "Cross-fitting equation `i' (`vname', `vtilde')" _c
				
				crossfit if `touse',					///
					eststring(`eststring')				///
					kfolds(`kfolds')					///
					foldvar(`mname'_fid_`m')			/// macro m is resampling counter
					vtilde(`vtilde'_`m')				///
					vtildeh(`vtilde_h'_`m')				/// LIE only
					eststringh(`eststring_h')			/// LIE only
					vname(`vname')						///
					`resid'

				// store MSE and sample size; also set eqn crossfitted flag = 1
				// assumes needed results from crossfit are in r(.) macros
				mata: add_to_eqn(`mname',`i')
				if ("`eqntype'"=="deq"&"`model'"=="optimaliv") {
					mata: add_to_eqn_h(`mname',`i')	
				}	
			}
	
			// for each equation: save names of tilde vars with smallest MSE
			// in each resample m and for each variable in nameY/listD/etc.
			// subroutine will handle case of empty lists
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(dheq) vlist(`listDH') m(`m') model(`model')
			_ddml_crossfit_update_optlist `mname', etype(zeq) vlist(`listZ') m(`m') model(`model')
		
		}	// end crossfitting block
	}	// end resampling block
	
	// report results by equation type with resamplings grouped together
	di
	di as res "Reporting crossfitting results:"
	_ddml_crossfit_report `mname'
	*/
	
end

program report_debugging

	syntax name(name=mname)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("eqnlist",invtokens(`mname'.eqnlistNames))
	
	*** report estimates for full sample for debugging purposes
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		local do_eqn	: list posof "`vtilde'" in eqnlist

		if `do_eqn' {
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			local 0 "`eststring'"
			syntax [anything] , [*]
			local est_main `anything'
			local est_options `options'
			if "`foldlist'"=="" {
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				di as res "(full sample, for debugging; no crossfit)"
				// estimate
				`est_main' if `mname'_sample, `est_options'
			}
			else {
				di
				di as res "Estimating equation `i', `vname'/`vtilde':"
				foreach k of numlist `foldlist' {
					di
					di as res "Fold=`k' (for debugging; no crossfit)"
					// estimate
					`est_main' if `mname'_fid!=`k' & `mname'_sample, `est_options'
				}
			}
		}
	}	

end

/*
mata:

// function to set crossfit dummy indicating whether crossfit has been done already
void set_crossfit(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar cf)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	(*p).crossfitted = cf
}

end
*/
