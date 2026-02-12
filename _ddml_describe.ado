*! ddml v1.5.0
*! last edited: 6feb2026
*! authors: aa/ms


program define _ddml_describe, rclass
	version 16
	syntax name(name=mname), [LEARNers CROSSfit ESTimates SAMple STACKing cvc rsq mse rmse n stweights ssweights psweights all mean median ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// used for storing matrices before returning
	tempname vmat
	
	local allflag	= "`all'"~=""
	if `allflag' {
		local learners		learners
		local crossfit		crossfit
		local estimates		estimates
		local sample		sample
		local stacking		stacking
		local cvc			cvc
		local rsq			rsq
		local mse			mse
		local rmse			rmse
		local n				n
		local stweights		stweights
		local ssweights		ssweights
		local psweights		psweights
	}
	local lflag		= "`learners'"~=""
	local cflag		= "`crossfit'"~=""
	local eflag		= "`estimates'"~=""
	local sflag		= "`sample'"~=""
	local stwflag	= "`stacking'"~=""
	local cvcflag	= "`cvc'"~=""
	local rsqflag	= "`rsq'"~=""
	local mseflag	= "`mse'"~=""
	local rmseflag	= "`rmse'"~=""
	local nflag		= "`n'"~=""
	local stflag	= "`stweights'"~=""	| `stwflag'
	local ssflag	= "`ssweights'"~=""	| `stwflag'
	local psflag	= "`psweights'"~=""	| `stwflag'
	
	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))	// flag for crossfitting results available
	mata: st_local("estimated",strofreal(`mname'.estimated))		// flag for estimation results available
	mata: st_local("ncombos",strofreal(`mname'.ncombos))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	// fold IDs
	forvalues m=1/`nreps' {
		local fidlist `fidlist' `mname'_fid_`m'
		local rslist `rslist' `mname'_sample_`m'
	}
	
	// syntax check - options require prior estimation
	local optlist eflag
	local optprob =0
	if `estimated'==0 {
		foreach opt in `optlist' {
			if ``opt'' {
				local optprob =1
				local `opt' = 0
			}
		}
		if `optprob'	di as res "some selected option(s) require prior estimation; not reported"
	}

	// syntax check - options required prior crossfitting
	local optlist cflag stwflag cvcflag rsqflag mseflag rmseflag nflag stflag ssflag psflag
	local optprob =0
	if `crossfitted'==0 {
		foreach opt in `optlist' {
			if ``opt'' {
				local optprob =1
				local `opt' = 0
			}
		}
		if `optprob'	di as res "some selected option(s) require prior crossfitting; not reported"
	}
	
	// syntax check - incompatible options
	if "`mean'"~="" & "`median'"~="" {
		di as err "incompatible options - mean and median"
		exit 198
	}
	else if "`mean'`median'"=="" {
		// default is mean
		local mean mean
	}
	local medmean `mean'`median'

	// basic info about equations and learners - always displayed
	di
	di as text "Model:" _col(25) as res "`model', crossfit folds k=" `kfolds' ", resamples r=" `nreps'
	di as text "Mata global (mname):" _col(25) as res "`mname'"
	mata: st_local("fclustvar",`mname'.fclustvar)
	if "`fclustvar'"~="" {
		di as res _col(25) "Folds respect clustering by `fclustvar'"
	}
	if "`nameY'"~="" {
		di as text "Dependent variable (Y):" _col(25) as res "`nameY'"
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
		local numlnrY : word count `vtlistY'
		di as text _col(2) "`nameY' learners:" as text _col(25) as res "`vtlistY'"
	}
	if `numeqnD' {
		di as text "D equations (`numeqnD'):" _col(25) as res "`nameD'"
		foreach var in `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			local numlnrD : word count `vtlistD'
			local comboD `comboD' * `numlnrD'
			di as text _col(2) "`var' learners:" _col(25) as res "`vtlistD'"
		}
	}
	if `numeqnZ' {
		di as text "Z equations (`numeqnZ'):" _col(25) as res "`nameZ'"
		foreach var in `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			local numlnrZ : word count `vtlistZ'
			local comboZ `comboZ' * `numlnrZ'
			di as text _col(2) "`var' learners:" _col(25) as res "`vtlistZ'"
		}
	}
	// report number of specifications
	if "`model'"=="interactive" | "`model'"=="interactiveiv" {
		local comboY `numlnrY' * `numlnrY'
	}
	else {
		local comboY `numlnrY'
	}
	if "`model'"=="interactiveiv" | "`model'"=="fiv" {
		local comboD `comboD' `comboD'
	}
	
	if `allflag' {
		di as text "Specifications:" _col(25) as res `comboY' `comboD' `comboZ' " possible specs" _c
		if `nreps' > 1 {
			di as res " * " `nreps' " crossfit splits = " `comboY' `comboD' `comboZ' * `nreps'
		}
		else {
			di
		}
	}
	
	// sample and folds in detail
	if `sflag' {
		di
		di as text "ID:" _col(25) as res "`mname'_id"
		di as text "Full sample indic.:" _col(25) as res "`mname'_sample" _c
		cap count if `mname'_sample
		if _rc==0 {
			di as res " (N=`r(N)')"
		}
		else {
			di as res " (N=n.a.)"
		}
		if "`fclustvar'"~="" {
			qui tab `fclustvar' if `mname'_sample
			di as text "Cluster variable:" _col(25) as res "`fclustvar' (N_clust=`r(r)')"
		}
		// di as text "Number of resamples =" _col(25) as res %3.0f `nreps'
		// di as text "Number of folds     =" _col(25) as res %3.0f `kfolds'
		di as text "Fold ID:" _col(25) as res _c
		forvalues m=1/`nreps' {
			local fid : word `m' of `fidlist'
			di as res %~12s "`fid'" _c
		}
		di
		di as text "Fold sample indic.:" _col(25) as res _c
		forvalues m=1/`nreps' {
			local rs : word `m' of `rslist'
			di as res %~12s "`rs'" _c
		}
		di
		di as text "Estimation N:" _col(25) as res _c
		forvalues m=1/`nreps' {
			local rs : word `m' of `rslist'
			cap count if `rs'
			if _rc==0 {
				local nobs "`: di %2.0f r(N)'"
			}
			else {
				local nobs "(n.a.)"
			}
			di as res %~12s "`nobs'" _c
		}
		di
	}
	
	// learners in detail
	if `lflag' {
		di
		di as text "Y learners (detail):"
		desc_learners `mname', vname(`nameY') etype(yeq)
		if `numeqnD' {
			di as text "D learners (detail):"
			foreach var in `nameD' {
				desc_learners `mname', vname(`var') etype(deq)
			}
		}
		if `numeqnZ' {
			di as text "Z learners (detail):"
			foreach var in `nameZ' {
				desc_learners `mname', vname(`var') etype(zeq)
			}
		}
	}
		
	// crossfit results in detail
	if `cflag' {
		di
		di as text "Crossfit results (detail):"
		desc_crossfit `mname', vname(`nameY') etype(yeq) header
		mat `vmat' = r(cfresults)
		return mat cfresults_`nameY' = `vmat'
		// should always be a D eqn
		if `numeqnD' {
			foreach var in `nameD' {
				desc_crossfit `mname', vname(`var') etype(deq)
				mat `vmat' = r(cfresults)
				return mat cfresults_`var' = `vmat'
			}
		}
		if `numeqnZ' {
			local znum 1
			foreach var in `nameZ' {
				desc_crossfit `mname', vname(`var') etype(zeq)
				mat `vmat' = r(cfresults)
				return mat cfresults_`var' = `vmat'
			}
		}
	}
	
	// stacking results in detail
	if `stwflag' {	
		di
		di as text "Stacking results (detail):"
		
		// always a Y eqn
		tempname stweights_Y1 ssweights_Y1 psweights_Y1
		desc_stacking `mname', vname(`nameY') etype(yeq) `medmean'
		mat `vmat' = r(stweights)
		return mat stweights_`nameY' = `vmat'
		mat `vmat' = r(ssweights)
		return mat ssweights_`nameY' = `vmat'
		mat `vmat' = r(psweights)
		return mat psweights_`nameY' = `vmat'
		// should always be a D eqn
		if `numeqnD' {
			foreach var in `nameD' {
				tempname stweights_D`dnum' ssweights_D`dnum' psweights_D`dnum'
				desc_stacking `mname', vname(`var') etype(deq) `medmean'
				mat `vmat' = r(stweights)
				return mat stweights_`var' = `vmat'
				mat `vmat' = r(ssweights)
				return mat ssweights_`var' = `vmat'
				mat `vmat' = r(psweights)
				return mat psweights_`var' = `vmat'
			}
		}
		// may or may not be a Z eqn
		if `numeqnZ' {
			local znum 1
			foreach var in `nameZ' {
				desc_stacking `mname', vname(`var') etype(zeq) `medmean'
				mat `vmat' = r(stweights)
				return mat stweights_`var' = `vmat'
				mat `vmat' = r(ssweights)
				return mat ssweights_`var' = `vmat'
				mat `vmat' = r(psweights)
				return mat psweights_`var' = `vmat'
			}
		}
	}	// end stacking block

	// CVC results in detail
	if `cvcflag' & `crossfitted' {
		di
        di as res `"CVC test ({browse "https://doi.org/10.1080/01621459.2019.1672556":Lei 2020}):"'
        di as res "H0: " as text "Learner has the lowest predictive risk among all candidate learners."
        di as res "HA: " as text "There is another learner with lower predictive risk."
        // bootstrap num can vary across vars or reps
        local allvars `nameY' `nameD' `nameZ'
        forvalues m=1/`nreps' {
        	foreach var in `allvars' {
				_ddml_extract bootnum, mname(`mname') vname(`var') key1(`var') key2("cvc_bootnum") key3("`m'") stata
				local prevbootnum `thisbootnum'
				local thisbootnum `r(bootnum)'
				// list of bootnums by rep and variable
				local cvcbootnum `cvcbootnum' `thisbootnum'
				// set flag
				if "`prevbootnum'"~="" & "`prevbootnum'"~="`thisbootnum'"	local cvcbootnum_msg "(avg; varies by rep/variable)"
			}
		}
		// if msg is blank, it's the same bootnum for every rep/variable
		if "`cvcbootnum_msg'"==""	local cvcbootnum = `thisbootnum'
		else						mata: st_local("cvcbootnum",strofreal(mean(strtoreal(tokens("`cvcbootnum'")'))))
		di as res %3.0f `cvcbootnum' as text " bootstrap reps `cvcbootnum_msg'"
        di
        di as res "CVC test p-values:"
		desc_values `mname', vname(`nameY') etype(yeq) cvc `medmean'
		mat `vmat' = r(values)
		return mat cvc_`nameY' = `vmat'
		return scalar cvcbootnum = `cvcbootnum'
		// should always be a D eqn
		if `numeqnD' {
			foreach var in `nameD' {
				desc_values `mname', vname(`var') etype(deq) cvc `medmean'
				mat `vmat' = r(values)
				return mat cvc_`var' = `vmat'
			}
		}
		if `numeqnZ' {
			foreach var in `nameZ' {
				desc_values `mname', vname(`var') etype(zeq) cvc `medmean'
				mat `vmat' = r(values)
				return mat cvc_`var' = `vmat'
			}
		}
	}

	// R-sq/MSE/RMSE/N results in detail
	if (`rsqflag' | `mseflag' | `rmseflag' | `nflag') & `crossfitted' {
		foreach stat in `rsq' `rmse' `mse' `n' {
			di
			if "`stat'"=="rsq"		di as res "R-sq by learner:"
			if "`stat'"=="mse"		di as res "MSE by learner:"
			if "`stat'"=="rmse"		di as res "RMSE by learner:"
			if "`stat'"=="n"		di as res "Sample size by learner:"
			desc_values `mname', vname(`nameY') etype(yeq) `stat' `medmean'
			mat `vmat' = r(values)
			return mat `stat'_`nameY' = `vmat'
			// should always be a D eqn
			if `numeqnD' {
				foreach var in `nameD' {
					desc_values `mname', vname(`var') etype(deq) `stat' `medmean'
					mat `vmat' = r(values)
					return mat `stat'_`var' = `vmat'
				}
			}
			if `numeqnZ' {
				foreach var in `nameZ' {
					desc_values `mname', vname(`var') etype(zeq) `stat' `medmean'
					mat `vmat' = r(values)
					return mat `stat'_`var' = `vmat'
				}
			}
		}
	}

	// estimate results in detail; notable option since _ddml_estimate routines would otherwise output this
	if `eflag' & ("`model'"=="interactive" | "`model'"=="interactiveiv") & `estimated' {
		di
		_ddml_estimate_ate_late `mname', `options' replay notable
	}
	else if `eflag' & `estimated' {
		di
		_ddml_estimate_linear `mname', `options' replay notable
	}
	
	// clear this global from Mata
	mata: mata drop `eqn'


end

prog define desc_stacking, rclass

	syntax name(name=mname), vname(string) etype(string) [ mean median ]	// etype is yeq, deq or zeq (not dheq)

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	mata: st_local("psflag",strofreal(`mname'.psflag))

	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	mata: st_local("poolstack",invtokens(`eqn'.poolstack))
	
	if ("`model'"=="interactive" | "`model'"=="interactiveiv") & "`etype'"=="yeq" {
		local dh d
	}
	if ("`model'"=="interactiveiv") & "`etype'"=="deq" {
		local dh d
	}
	if ("`model'"=="fiv") & "`etype'"=="deq" {
		local dh h
	}
	
	local medmean `mean'`median'

	tempname wmat
	
	di
	di as text "Conditional expectation: " as res "`vname'"
	
	// standard stacking
	if `stdflag' {
		local vtilde `vtlist'
		extract_weights `mname', vname(`vname') etype(`etype') stweights `medmean'
		mat `wmat' = r(weights)
		di
		di as text "Standard stacking weights: " _c
		if `wmat'[1,1] ~= . {
			di as res "`vtilde'"
			display_values `wmat', dh(`dh') nreps(`nreps') `medmean'
			di as text "(nb: std stacking weights above are `medmean's across `kfolds' folds)"
		}
		else {
			// don't display varname if weights not available
			di as res "(no weights available)"
		}
		return matrix stweights=`wmat'
	}

	// short-stacking
	if `ssflag' {
		extract_weights `mname', vname(`vname') etype(`etype') ssweights `medmean'
		mat `wmat' = r(weights)
		di
		di as text "Short-stacking weights: " as res "`shortstack'_ss"
		display_values `wmat', dh(`dh') nreps(`nreps') `medmean'
		return matrix ssweights=`wmat'
	}

	// pooled stacking
	if `psflag' {
		extract_weights `mname', vname(`vname') etype(`etype') psweights `medmean'
		mat `wmat' = r(weights)
		di
		di as text "Pooled stacking weights: " _c
		if `wmat'[1,1] ~= . {
			di as res "`poolstack'_ps"
			display_values `wmat', dh(`dh') nreps(`nreps') `medmean'
		}
		else {
			// don't display varname if weights not available
			di as res "(no weights available)"
		}
		return matrix psweights=`wmat'
	}


end


prog define extract_weights, rclass

	syntax name(name=mname), vname(string) etype(string)		/// etype is yeq, deq or zeq (not dheq)
								[								///
								stweights ssweights psweights	///
								mean median						///
								 ]

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))

	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("pystackedmulti",strofreal(`eqn'.pystackedmulti))
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	mata: st_local("poolstack",invtokens(`eqn'.poolstack))
	
	if "`etype'"=="yeq" {
		local ydz y
	}
	else if "`etype'"=="deq" {
		local ydz d
	}
	else if "`etype'"=="zeq" {
		local ydz z
	}
	
	if "`stweights'"=="stweights"		{
		local stflag	=1
		local ssflag	=0
		local psflag	=0
		local key1		`vtlist'
		local key2		stack_weights
	}
	else if "`ssweights'"=="ssweights" {
		local stflag	=0
		local ssflag	=1
		local psflag	=0
		local key1		`shortstack'_ss
		local key2		ss_weights
	}
	else if "`psweights'"=="psweights" {
		local stflag	=0
		local ssflag	=0
		local psflag	=1
		local key1		`poolstack'_ps
		local key2		ps_weights
	}
	else {
		di as err "ddml extract error"
		exit 198
	}

	if "`median'"=="median"		local medmean	median
	else if "`mean'"=="mean"	local medmean	mean
	else {
		di as err "ddml describe error"
		exit 198
	}
	
	tempname wmat
	tempname val val_m val_0 val_1 val_h val_m_h

	if ("`model'"=="interactive" & "`ydz'"=="y") 	|	///
		("`model'"=="interactiveiv" & "`ydz'"=="y")	|	///
		("`model'"=="interactiveiv" & "`ydz'"=="d")		///
		{
		
		// interactive, treat=0 and treat=1
				
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`key1') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		forvalues m=1/`nreps' {
			forvalues i=0/1 {
				_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'`i') key3(`m') stata
				mat `val_m' = r(val)
				if `stflag' {
					// for std stacking, need to convert into means/medians across folds
					if "`medmean'"=="mean"	mata: st_matrix("`val_m'", mean(st_matrix("`val_m'")')')
					else					mata: st_matrix("`val_m'", ddml_median(st_matrix("`val_m'")')')
				}
				mat `val_`i'' = nullmat(`val_`i'') , `val_m'
			}
		}
		forvalues i=0/1 {
			mata: st_matrix("`val_`i''",(mean(st_matrix("`val_`i''")')' , st_matrix("`val_`i''")))
			mat `val_`i'' = J(`nlearners',1,`i') , `val_`i''
			mata: st_matrix("`val_`i''",(runningsum(J(`nlearners',1,1)), st_matrix("`val_`i''")))
			// rownames
			mat rownames `val_`i'' = `lnames'
		}
		mat `val' = `val_0' \ `val_1'
		// colnames
		if "`model'"=="interactive"		local cnames learner D=0/1 mean_weight
		else							local cnames learner Z=0/1 mean_weight
		forvalues i=1/`nreps' {
			local cnames `cnames' rep_`i'
		}
		mat colnames `val' = `cnames'
	}
	else if "`ydz'"=="d" & "`model'"=="fiv" {
		// fiv, with/without h		
		cap mat drop `val'
		
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est_h) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		forvalues m=1/`nreps' {
			_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2') key3(`m') stata
			mat `val_m' = r(val)
			_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'_h) key3(`m') stata
			mat `val_m_h' = r(val)
			if `stflag' {
				// for std stacking, need to convert into means/medians across folds
				if "`medmean'"=="mean" {
					mata: st_matrix("`val_m'", mean(st_matrix("`val_m'")')')
					mata: st_matrix("`val_m_h'", mean(st_matrix("`val_m_h'")')')
				}
				else {
					mata: st_matrix("`val_m'", ddml_median(st_matrix("`val_m'")')')
					mata: st_matrix("`val_m_h'", ddml_median(st_matrix("`val_m_h'")')')
				}
			}
			mat `val' = nullmat(`val') , `val_m'
			mat `val_h' = nullmat(`val_h') , `val_m_h'
		}
		// add learner number, h=0/1 and mean as first column
		if "`medmean'"=="mean" {
			mata: st_matrix("`val'", (runningsum(J(`nlearners',1,1)), J(`nlearners',1,0), mean(st_matrix("`val'")')' , st_matrix("`val'")))
			mata: st_matrix("`val_h'", (runningsum(J(`nlearners',1,1)), J(`nlearners',1,1), mean(st_matrix("`val_h'")')' , st_matrix("`val_h'")))
		}
		else {
			mata: st_matrix("`val'", (runningsum(J(`nlearners',1,1)), J(`nlearners',1,0), ddml_median(st_matrix("`val'")')' , st_matrix("`val'")))
			mata: st_matrix("`val_h'", (runningsum(J(`nlearners',1,1)), J(`nlearners',1,1), ddml_median(st_matrix("`val_h'")')' , st_matrix("`val_h'")))
		}
		// rownames
		mat rownames `val' = `lnames'
		mat rownames `val_h' = `lnames'
		// append
		mat `val' = `val' \ `val_h'
		// colnames
		local cnames learner h=0/1 mean_weight
		forvalues i=1/`nreps' {
			local cnames `cnames' rep_`i'
		}
		mat colnames `val' = `cnames'
	}
	else {
		// standard
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`key1') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
	
		// get values for each resample; rows are learners, columns are resamples
		// or folds/resamples if standard stacking weights
		forvalues m=1/`nreps' {
			_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2') key3(`m') stata
			mat `val_m' = r(val)
			if `stflag' {
				// for std stacking, need to convert into means across folds
				if "`medmean'"=="mean"		mata: st_matrix("`val_m'", mean(st_matrix("`val_m'")')')
				else						mata: st_matrix("`val_m'", ddml_median(st_matrix("`val_m'")')')
			}
			mat `val' = nullmat(`val') , `val_m'
		}
		// add learner number and mean as first column
		if "`medmean'"=="mean"		mata: st_matrix("`val'",(runningsum(J(`nlearners',1,1)), mean(st_matrix("`val'")')' , st_matrix("`val'")))
		else						mata: st_matrix("`val'",(runningsum(J(`nlearners',1,1)), ddml_median(st_matrix("`val'")')' , st_matrix("`val'")))
		// rownames and colnames
		mat rownames `val' = `lnames'
		local cnames learner mean_weight
		forvalues i=1/`nreps' {
			local cnames `cnames' rep_`i'
		}
		mat colnames `val' = `cnames'
	}

	return mat weights	= `val'
	
end


prog define desc_values, rclass

	syntax name(name=mname), vname(string) etype(string)	/// etype is yeq, deq or zeq (not dheq)
								[							///
								 cvc rsq mse rmse n			///
								 mean median				///
								 ]

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	mata: st_local("psflag",strofreal(`mname'.psflag))

	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	mata: st_local("poolstack",invtokens(`eqn'.poolstack))
	mata: st_local("pystackedmulti",strofreal(`eqn'.pystackedmulti))
	
	if "`etype'"=="yeq" {
		local ydz y
	}
	else if "`etype'"=="deq" {
		local ydz d
	}
	else if "`etype'"=="zeq" {
		local ydz z
	}
	
	if "`median'"=="median"		local medmean	median
	else if "`mean'"=="mean"	local medmean	mean
	else {
		di as err "ddml describe error"
		exit 198
	}
	
	// pystackedmulti means 1 call to pystacked with multiple learners, so vtlist is vtilde
	if `pystackedmulti'			local vtilde `vtlist'
	
	// cvc saved under vname; others saved under learner
	if "`cvc'"=="cvc"			local key1 `vname'
	else						local key1 `vtilde'
	
	if "`cvc'"=="cvc"			local key2 cvc_pval
	else if "`rsq'"=="rsq"		local key2 R-sq
	else if "`mse'"=="mse"		local key2 MSE
	else if "`rmse'"=="rmse"	local key2 RMSE
	else if "`n'"=="n"			local key2 N
	else {
		di as err "ddml extract error"
		exit 198
	}
	if "`n'"=="n"				local nodecimal nodecimal

	tempname wmat
	tempname val val_0 val_1 val_st val_st_0 val_st_1 val_ss val_ss_0 val_ss_1 val_ps val_ps_0 val_ps_1
	
	di
	di as text "Conditional expectation " _c
	if "`ydz'"=="y" {
		di as text "E[Y|X]: " _c
	}
	else if "`ydz'"=="d" {
		di as text "E[D|Z,X]: " _c
	}
	else if "`ydz'"=="z" {
		di as text "E[Z|X]: " _c
	}
	di as res "`vname'"

	if ("`model'"=="interactive" & "`ydz'"=="y") 	|	///
		("`model'"=="interactiveiv" & "`ydz'"=="y")	|	///
		("`model'"=="interactiveiv" & "`ydz'"=="d")		///
		{
		
		// interactive, treat=0 and treat=1
				
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		// get values for each resample; rows are learners, columns are resamples
		forvalues m=1/`nreps' {
			forvalues i=0/1 {
				if "`cvc'"=="cvc" {
				// single matrix of cvc results for all learners
					_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'`i') key3(`m') stata
					mat `val_`i'' = nullmat(`val_`i'') , r(val)
				}
				else if `pystackedmulti' {
					// pystacked learners all in one matrix
					_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'_L`i') key3(`m') stata
					mat `val_`i'' = nullmat(`val_`i'') , r(val)
				}
				else {
					// loop through learners
					tempname val_j
					foreach vt in `vtlist' {
						_ddml_extract val, mname(`mname') vname(`vname') key1(`vt') key2(`key2'`i') key3(`m') stata
						mat `val_j' = nullmat(`val_j') \ r(val)
					}
					mat `val_`i'' = nullmat(`val_`i''), `val_j'
				}
			}
		}
		// add rows for stacking learners unless cvc
		if `stdflag' & "`cvc'"=="" {
			forvalues i=0/1 {
				forvalues m=1/`nreps' {
					_ddml_extract val, mname(`mname') vname(`vname') key1(`vtlist') key2(`key2'`i') key3(`m') stata
					mat `val_st_`i'' = nullmat(`val_st_`i'') , r(val)
				}
				mat `val_`i'' = `val_`i'' \ `val_st_`i''
			}
			local lnames `lnames' Std_stacking
		}
		if `ssflag' & "`cvc'"=="" {
			forvalues i=0/1 {
				forvalues m=1/`nreps' {
					_ddml_extract val, mname(`mname') vname(`vname') key1(`shortstack'_ss) key2(`key2'`i') key3(`m') stata
					mat `val_ss_`i'' = nullmat(`val_ss_`i'') , r(val)
				}
				mat `val_`i'' = `val_`i'' \ `val_ss_`i''
			}
			local lnames `lnames' Short_stacking
		}
		if `psflag' & "`cvc'"=="" {
			forvalues i=0/1 {
				forvalues m=1/`nreps' {
					_ddml_extract val, mname(`mname') vname(`vname') key1(`poolstack'_ps) key2(`key2'`i') key3(`m') stata
					mat `val_ps_`i'' = nullmat(`val_ps_`i'') , r(val)
				}
				mat `val_`i'' = `val_`i'' \ `val_ps_`i''
			}
			local lnames `lnames' Pooled_stacking
		}
		// add mean/median as first column
		local nrows : list sizeof lnames
		forvalues i=0/1 {
			if "`medmean'"=="median" {
				mata: st_matrix("`val_`i''",(ddml_median(st_matrix("`val_`i''")')' , st_matrix("`val_`i''")))
			}
			else {
				mata: st_matrix("`val_`i''",(mean(st_matrix("`val_`i''")')' , st_matrix("`val_`i''")))
			}
			mat `val_`i'' = J(`nrows',1,`i') , `val_`i''
			mata: st_matrix("`val_`i''",(runningsum(J(`nrows',1,1)), st_matrix("`val_`i''")))
			mat rownames `val_`i'' = `lnames'
		}
		mat `val' = `val_0' \ `val_1'
		display_values `val', dh(d) nreps(`nreps') `medmean' `nodecimal'
	}
	else {
		// standard
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtilde') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		// get values for each resample; rows are learners, columns are resamples
		forvalues m=1/`nreps' {
			if "`cvc'"=="cvc" {
				// single matrix of cvc results for all learners
				_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2') key3(`m') stata
				mat `val' = nullmat(`val') , r(val)
			}
			else if `pystackedmulti' {
				// pystacked learners all in one matrix
				_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'_L) key3(`m') stata
				mat `val' = nullmat(`val') , r(val)
			}
			else {
				// loop through learners
				tempname val_j
				foreach vt in `vtlist' {
					_ddml_extract val, mname(`mname') vname(`vname') key1(`vt') key2(`key2') key3(`m') stata
					mat `val_j' = nullmat(`val_j') \ r(val)
				}
				mat `val' = nullmat(`val'), `val_j'
			}
		}
		// add rows for stacking learners unless cvc
		if `stdflag' & "`cvc'"=="" {
			forvalues m=1/`nreps' {
				_ddml_extract val, mname(`mname') vname(`vname') key1(`vtilde') key2(`key2') key3(`m') stata
				mat `val_st' = nullmat(`val_st') , r(val)
			}
			mat `val' = `val' \ `val_st'
			local lnames `lnames' Std_stacking
		}
		if `ssflag' & "`cvc'"=="" {
			forvalues m=1/`nreps' {
				_ddml_extract val, mname(`mname') vname(`vname') key1(`shortstack'_ss) key2(`key2') key3(`m') stata
				mat `val_ss' = nullmat(`val_ss') , r(val)
			}
			mat `val' = `val' \ `val_ss'
			local lnames `lnames' Short_stacking
		}
		if `psflag' & "`cvc'"=="" {
			forvalues m=1/`nreps' {
				_ddml_extract val, mname(`mname') vname(`vname') key1(`poolstack'_ps) key2(`key2') key3(`m') stata
				mat `val_ps' = nullmat(`val_ps') , r(val)
			}
			mat `val' = `val' \ `val_ps'
			local lnames `lnames' Pooled_stacking
		}
		// add mean/median as first column
		local nrows : list sizeof lnames
		if "`medmean'"=="median" {
			mata: st_matrix("`val'",(runningsum(J(`nrows',1,1)), ddml_median(st_matrix("`val'")')' , st_matrix("`val'")))
		}
		else {
			mata: st_matrix("`val'",(runningsum(J(`nrows',1,1)), mean(st_matrix("`val'")')' , st_matrix("`val'")))
		}
		mat rownames `val' = `lnames'
		display_values `val', nreps(`nreps') `medmean' `nodecimal'
	}
	if "`ydz'"=="d" & "`model'"=="fiv" {
		
		cap mat drop `val'
		
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est_h) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		// get values for each resample; rows are learners, columns are resamples
		forvalues m=1/`nreps' {
			if "`cvc'"=="cvc" {
				// single matrix of cvc results for all learners
				_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'_h) key3(`m') stata
				mat `val' = nullmat(`val') , r(val)
			}
			else if `pystackedmulti' {
				// pystacked learners all in one matrix
				_ddml_extract val, mname(`mname') vname(`vname') key1(`key1') key2(`key2'_L_h) key3(`m') stata
				mat `val' = nullmat(`val') , r(val)
			}
			else {
				// loop through learners
				tempname val_j
				foreach vt in `vtlist' {
					_ddml_extract val, mname(`mname') vname(`vname') key1(`vt') key2(`key2'_h) key3(`m') stata
					mat `val_j' = nullmat(`val_j') \ r(val)
				}
				mat `val' = nullmat(`val'), `val_j'
			}
		}
		// add mean/median as first column
		if "`medmean'"=="median" {
			mata: st_matrix("`val'",(runningsum(J(`nlearners',1,1)), ddml_median(st_matrix("`val'")')' , st_matrix("`val'")))
		}
		else {
			mata: st_matrix("`val'",(runningsum(J(`nlearners',1,1)), mean(st_matrix("`val'")')' , st_matrix("`val'")))
		}
		mat rownames `val' = `lnames'
		di as text "Conditional expectation E[D|X]: " as res "`vname'"
		display_values `val', nreps(`nreps') `medmean' `nodecimal'
	}
	
	// column names
	local cnames learner `medmean'
	forvalues m=1/`nreps' {
		local cnames `cnames' rep_`m'
	}
	mat colnames `val' = `cnames'
	
	return matrix values	= `val'

end


prog define display_values

	syntax name(name=wmat), [ dh(string) nreps(integer 1) mean median NOdecimal ]
	
	local rnames : rownames `wmat'
	local nrows = rowsof(`wmat')

	if "`mean'"=="mean"				local medmean "  Mean"
	else if "`median'"=="median"	local medmean "Median"
	else {
		di as err "internal ddml describe error
		exit 198
	}

	// check if values available
	if el(`wmat',1,1)==. {
		di as res "(not available)"
		exit
	}
	
	mata: st_numscalar("r(max)",max(abs(st_matrix("`wmat'"))))
	if "`nodecimal'"=="nodecimal"	local fmt 10.0fc
	else if r(max) < 10^3			local fmt 10.3f
	else if r(max) >= 10^6			local fmt 10.2e
	else							local fmt 10.0fc
	
	if "`dh'"=="" {
		di as text "Learner" _col(30) "`medmean'" as res "   |" _c
		local mloc=26
		local dhcol=0
	}
	else if "`dh'"=="d" {
		di as text "Learner" _col(23) "Treated" _col(38) "Mean" as res "   |" _c
		local mloc=32
		local dhcol=1
	}
	else if "`dh'"=="h" {
		di as text "Learner" _col(22) "Cond exp" _col(38) "Mean" as res "   |" _c
		local mloc=32
		local dhcol=1
	}
	forvalues j=1/`nreps' {
		di as text %10s "Rep `j'" _c
	}
	di

	forvalues i=1/`nrows' {
		local lname : word `i' of `rnames'
		di as res %2.0f el(`wmat',`i',1) ". " _c
		di as res "`lname'" _c
		if "`dh'"=="d" {
			di as res _col(26) %1.0f el(`wmat',`i',2) _c
		}
		else if "`dh'"=="h" {
			if el(`wmat',`i',2)==0 {
				local cetext E[D|X,Z]
			}
			else {
				local cetext E[D^|X]
			}
			di as res _col(22) %8s "`cetext'" _c
		}
		di as res _col(`mloc') %`fmt' el(`wmat',`i',2+`dhcol') "   |" _c
		forvalues j=1/`nreps' {
			di as res %`fmt' el(`wmat',`i',2+`dhcol'+`j') _c
		}
		di
	}

end



prog define desc_learners, rclass

	syntax name(name=mname), vname(string) etype(string)	// etype is yeq, deq or zeq (not dheq)

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	mata: st_local("psflag",strofreal(`mname'.psflag))
	
	if ("`etype'"=="deq") & ("`model'"=="fiv") {
		// includes both deq and dheq
		local heqn	= 1
	}
	else {
		local heqn	= 0
	}

	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))

	foreach vt in `vtlist' {
		di as res _col(2) "Learner:" _col(15) "`vt'"
		mata: st_local("estring", return_learner_item(`eqn',"`vt'","estring"))
		// remove tabs and extraneous spaces
		local estring = subinstr("`estring'","	"," ",.)
		local estring = strtrim(stritrim("`estring'"))
		if `heqn' {
			di as res _col(15) "est cmd (D): `estring'"
		}
		else {
			di as res _col(15) "est cmd: `estring'"
		}
		if `heqn' {
			mata: st_local("estring_h", return_learner_item(`eqn',"`vt'","estring_h"))
			di as res _col(15) "est cmd (H): `estring_h'"
		}
	}

	// clear this global from Mata
	mata: mata drop `eqn'

end

prog define desc_crossfit, rclass

	syntax name(name=mname), vname(string) etype(string) [ header ]	// etype is yeq, deq or zeq (not dheq)
	
	local vnlen = 15
	local lnrlen = 20
	
	local showheader = "`header'"~=""
	local vnabbrev = abbrev("`vname'",`vnlen')

	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	mata: st_local("psflag",strofreal(`mname'.psflag))
	local nostackflag = ~`stdflag' & ~`ssflag' & ~`psflag'
	
	// used below to indicate set of crossfitting results to report
	local pairs		= 0
	local heqn		= 0
	if ("`etype'"=="yeq") & ("`model'"=="interactive" | "`model'"=="interactiveiv") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="interactiveiv") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="fiv") {
		// includes both deq and dheq
		local heqn	= 1
	}
	// extra spacing required
	if ("`model'"=="interactive" | "`model'"=="interactiveiv") {
		local tv = 3
	}
	else {
		local tv = 0
	}

	if `showheader' {
		local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8+4+8
		di as text _col(`c') "MSE by fold:" _c
		di
		di as text "Cond. exp." _c
		local c = `vnlen'+3
		di as text _col(`c') "Learner" _c
		local c = `vnlen'+2+`lnrlen'+1
		di as text _col(`c') "rep" _c
		if `pairs' {
			local c = `vnlen'+2+`lnrlen'+1+`tv'+1
			di as text _col(`c') "tv" _c
		}
		local c = `vnlen'+2+`lnrlen'+1+`tv'+5+1
		di as text _col(`c') "R-sq" _c
		local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8+4
		di as text _col(`c') "MSE" _c
		forvalues k=1/`kfolds' {
			di as text "        `k' " _c
		}
		di
	}
	
	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	mata: st_local("poolstack",invtokens(`eqn'.poolstack))

	tempname cfresults cfresults_r cfresults0_r cfresults1_r
	// initialize
	local rnames
	local firstrow = 1
	mat `cfresults_r' = .
	mat `cfresults0_r' = .
	mat `cfresults1_r' = .
	foreach vt in `vtlist' {
		// for pystacked stacking and non-pystacked
		if `crossfitted' & (`stdflag' | `nostackflag') {
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse_folds
					mata: st_local("rsq", strofreal(return_result_item(`eqn',"`vt'","R-sq","`m'")))
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`vt'","MSE","`m'")))
					mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`vt'","MSE_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					local lnrname `vt'_`m'
					local lnrabbrev = abbrev("`lnrname'",`lnrlen')
					local c = `vnlen'+2
					di as res _col(`c') "`lnrabbrev'" _c
					local c = `vnlen'+2+`lnrlen'+1
					di _col(`c') %2.0f `m' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5
					di _col(`c') %5.2f `rsq' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
					di _col(`c') %8.2f `mse' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse_folds',1,`k') _c
					}
					di
					// for r(.)
					local rnames `rnames' `lnrname'
					mat `cfresults_r' = `m', `rsq', `mse', `mse_folds'
					mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
				}
			}
			else {
				forvalues m=1/`nreps' {
					tempname mse0_folds mse1_folds
					mata: st_local("rsq0", strofreal(return_result_item(`eqn',"`vt'","R-sq0","`m'")))
					mata: st_local("rsq1", strofreal(return_result_item(`eqn',"`vt'","R-sq1","`m'")))
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`vt'","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`vt'","MSE1","`m'")))
					mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`vt'","MSE0_folds","`m'"))
					mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`vt'","MSE1_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					forvalues i=0/1 {
						local lnrname `vt'`i'_`m'
						local lnrabbrev = abbrev("`lnrname'",`lnrlen')
						local c = `vnlen'+2
						di as res _col(`c') "`lnrabbrev'" _c
						local c = `vnlen'+2+`lnrlen'+1
						di _col(`c') %2.0f `m' _c
						local c = `vnlen'+2+`lnrlen'+1+4
						di _col(`c') %2.0f `i' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4
						di _col(`c') %5.2f `rsq`i'' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4+8
						di _col(`c') %8.2f `mse`i'' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse`i'_folds',1,`k') _c
						}
						di
						// for r(.)
						local rnames `rnames' `lnrname'
						mat `cfresults`i'_r' = `m', `rsq`i'', `mse`i'', `mse`i'_folds'
					}
					mat `cfresults' = nullmat(`cfresults') \ `cfresults0_r' \ `cfresults1_r'
				}
			}
		}
		// fiv model, for pystacked stacking and non-pystacked
		if `heqn' & `crossfitted' & (`stdflag' | `nostackflag') {
			forvalues m=1/`nreps' {
				tempname mse_h_folds
				mata: st_local("rsq_h", strofreal(return_result_item(`eqn',"`vt'","R-sq_h","`m'")))
				mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`vt'","MSE_h","`m'")))
				mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`vt'","MSE_h_folds","`m'"))
				local lnrname `vt'_h_`m'
				local lnrabbrev = abbrev("`lnrname'",`lnrlen')
				local c = `vnlen'+2
				di as res _col(`c') "`lnrabbrev'" _c
				local c = `vnlen'+2+`lnrlen'+1
				di _col(`c') %2.0f `m' _c
				local c = `vnlen'+2+`lnrlen'+1+`tv'+5
				di _col(`c') %5.2f `rsq_h' _c
				local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
				di _col(`c') %8.2f `mse_h' _c
				forvalues k=1/`kfolds' {
					di "  " %8.2f el(`mse_h_folds',1,`k') _c
				}
				di
				// for r(.)
				local rnames `rnames' `lnrname'
				mat `cfresults_r' = `cfresults_r' \ (`m', `rsq_h', `mse_h', `mse_h_folds')
				mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
			}
		}
	}

	// short-stacking
	if `ssflag' {
		if `crossfitted' {
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse_folds
					mata: st_local("rsq", strofreal(return_result_item(`eqn',"`shortstack'_ss","R-sq","`m'")))
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE","`m'")))
					mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					local lnrname `shortstack'_ss_`m'
					local lnrabbrev = abbrev("`lnrname'",`lnrlen')
					local c = `vnlen'+2
					di as res _col(`c') "`lnrabbrev'" _c
					local c = `vnlen'+2+`lnrlen'+1
					di _col(`c') %2.0f `m' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5
					di _col(`c') %5.2f `rsq' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
					di _col(`c') %8.2f `mse' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse_folds',1,`k') _c
					}
					di
					// for r(.)
					local rnames `rnames' `lnrname'
					mat `cfresults_r' = `m', `rsq', `mse', `mse_folds'
					mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
				}
				if `heqn' {
					forvalues m=1/`nreps' {
						tempname mse_h_folds
						mata: st_local("rsq", strofreal(return_result_item(`eqn',"`shortstack'_ss","R-sq_h","`m'")))
						mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE_h","`m'")))
						mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE_h_folds","`m'"))
						local c = `vnlen'+2
						local lnrname `shortstack'_h_ss_`m'
						local lnrabbrev = abbrev("`lnrname'",`lnrlen')
						di as res _col(`c') "`lnrabbrev'" _c
						local c = `vnlen'+2+`lnrlen'+1
						di _col(`c') %2.0f `m' _c
						local c = `vnlen'+2+`lnrlen'+1+`tv'+5
						di _col(`c') %5.2f `rsq_h' _c
						local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
						di _col(`c') %8.2f `mse_h' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse_h_folds',1,`k') _c
						}
						di
						// for r(.)
						local rnames `rnames' `lnrname'
						mat `cfresults_r' = `cfresults_r' \ (`m', `rsq_h', `mse_h', `mse_h_folds')
						mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
					}
				}
			}
			else {
				forvalues m=1/`nreps' {
					tempname mse0_folds mse1_folds
					mata: st_local("rsq0", strofreal(return_result_item(`eqn',"`shortstack'_ss","R-sq0","`m'")))
					mata: st_local("rsq1", strofreal(return_result_item(`eqn',"`shortstack'_ss","R-sq1","`m'")))
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`shortstack'_ss","MSE1","`m'")))
					mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE0_folds","`m'"))
					mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`shortstack'_ss","MSE1_folds","`m'"))
					forvalues i=0/1 {
						local lnrname `shortstack'_ss`i'_`m'
						local lnrabbrev = abbrev("`lnrname'",`lnrlen')
						local c = `vnlen'+2
						di as res _col(`c') "`lnrabbrev'" _c
						local c = `vnlen'+2+`lnrlen'+1
						di _col(`c') %2.0f `m' _c
						local c = `vnlen'+2+`lnrlen'+1+4
						di _col(`c') %2.0f `i' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4
						di _col(`c') %5.2f `rsq`i'' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4+8
						di _col(`c') %8.2f `mse`i'' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse`i'_folds',1,`k') _c
						}
						di
						// for r(.)
						local rnames `rnames' `lnrname'
						mat `cfresults`i'_r' = `m', `rsq`i'', `mse`i'', `mse`i'_folds'
					}
					mat `cfresults' = nullmat(`cfresults') \ `cfresults0_r' \ `cfresults1_r'
				}
			}
		}
	}

	// pooled stacking
	if `psflag' {
		if `crossfitted' {
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse_folds
					mata: st_local("rsq", strofreal(return_result_item(`eqn',"`poolstack'_ps","R-sq","`m'")))
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`poolstack'_ps","MSE","`m'")))
					mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`poolstack'_ps","MSE_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					local lnrname `poolstack'_ps`i'_`m'
					local lnrabbrev = abbrev("`lnrname'",`lnrlen')
					local c = `vnlen'+2
					di as res _col(`c') "`lnrabbrev'" _c
					local c = `vnlen'+2+`lnrlen'+1
					di _col(`c') %2.0f `m' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5
					di _col(`c') %5.2f `rsq' _c
					local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
					di _col(`c') %8.2f `mse' _c
					forvalues k=1/`kfolds' {
						di "  " %8.2f el(`mse_folds',1,`k') _c
					}
					di
					// for r(.)
					local rnames `rnames' `lnrname'
					mat `cfresults_r' = `m', `rsq', `mse', `mse_folds'
					mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
				}
				if `heqn' {
					forvalues m=1/`nreps' {
						tempname mse_h_folds
						mata: st_local("rsq", strofreal(return_result_item(`eqn',"`poolstack'_ps","R-sq_h","`m'")))
						mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`poolstack'_ps","MSE_h","`m'")))
						mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`poolstack'_ps","MSE_h_folds","`m'"))
						local c = `vnlen'+2
						local lnrname `poolstack'_h_`m'
						local lnrabbrev = abbrev("`lnrname'",`lnrlen')
						di as res _col(`c') "`lnrabbrev'" _c
						local c = `vnlen'+2+`lnrlen'+1
						di _col(`c') %2.0f `m' _c
						local c = `vnlen'+2+`lnrlen'+1+`tv'+5
						di _col(`c') %5.2f `rsq_h' _c
						local c = `vnlen'+2+`lnrlen'+1+`tv'+5+8
						di _col(`c') %8.2f `mse_h' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse_h_folds',1,`k') _c
						}
						di
						// for r(.)
						local rnames `rnames' `lnrname'
						mat `cfresults_r' = `cfresults_r' \ (`m', `rsq_h', `mse_h', `mse_h_folds')
						mat `cfresults' = nullmat(`cfresults') \ `cfresults_r'
					}
				}
			}
			else {
				forvalues m=1/`nreps' {
					tempname mse0_folds mse1_folds
					mata: st_local("rsq0", strofreal(return_result_item(`eqn',"`poolstack'_ps","R-sq0","`m'")))
					mata: st_local("rsq1", strofreal(return_result_item(`eqn',"`poolstack'_ps","R-sq1","`m'")))
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`poolstack'_ps","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`poolstack'_ps","MSE1","`m'")))
					mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`poolstack'_ps","MSE0_folds","`m'"))
					mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`poolstack'_ps","MSE1_folds","`m'"))
					forvalues i=0/1 {
						local c = `vnlen'+2
						local lnrname `poolstack'_ps`i'_`m'
						local lnrabbrev = abbrev("`lnrname'",`lnrlen')
						di as res _col(`c') "`lnrabbrev'" _c
						local c = `vnlen'+2+`lnrlen'+1
						di _col(`c') %2.0f `m' _c
						local c = `vnlen'+2+`lnrlen'+1+4
						di _col(`c') %2.0f `i' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4
						di _col(`c') %5.2f `rsq`i'' _c
						local c = `vnlen'+2+`lnrlen'+1+4+4+8
						di _col(`c') %8.2f `mse`i'' _c
						forvalues k=1/`kfolds' {
							di "  " %8.2f el(`mse`i'_folds',1,`k') _c
						}
						di
						// for r(.)
						local rnames `rnames' `lnrname'
						mat `cfresults`i'_r' = `m', `rsq`i'', `mse`i'', `mse`i'_folds'
					}
					mat `cfresults' = nullmat(`cfresults') \ `cfresults0_r' \ `cfresults1_r'
				}
			}
		}
	}
	
	// clear this global from Mata
	mata: mata drop `eqn'
	
	mat rownames `cfresults' = `rnames'
	local cnames  rep Rsq MSE
	forvalues k=1/`kfolds' {
		local cnames `cnames' MSE_`k'
	}
	mat colnames `cfresults' = `cnames'
	
	return mat cfresults	= `cfresults'
	
end


