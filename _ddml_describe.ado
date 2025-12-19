*! ddml v1.5.0
*! last edited: 18dec2025
*! authors: aa/ms

program define _ddml_describe, rclass
	version 16
	syntax name(name=mname), [LEARNers CROSSfit ESTimates SAMple STACKing cvc all *]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	local allflag	= "`all'"~=""
	local lflag		= "`learners'"~=""	| `allflag'
	local cflag		= "`crossfit'"~=""	| `allflag'
	local eflag		= "`estimates'"~=""	| `allflag'
	local sflag		= "`sample'"~=""	| `allflag'
	local stflag	= "`stacking'"~=""	| `allflag'
	local cvcflag	= "`cvc'"~=""		| `allflag'
	
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
				local N "`: di %2.0f r(N)'"
			}
			else {
				local N "(n.a.)"
			}
			di as res %~12s "`N'" _c
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
	if `cflag' & `crossfitted' {
		di
		di as text "Crossfit results (detail):"
		desc_crossfit `mname', vname(`nameY') etype(yeq) header
		tempname cfresults_Y1 
		mat `cfresults_Y1' = r(cfresults)
		// should always be a D eqn
		if `numeqnD' {
			local dnum 1
			foreach var in `nameD' {
				tempname cfresults_D`dnum'
				desc_crossfit `mname', vname(`var') etype(deq)
				mat `cfresults_D`dnum'' = r(cfresults)
				local ++dnum
			}
		}
		if `numeqnZ' {
			local znum 1
			foreach var in `nameZ' {
				tempname cfresults_Z`znum'
				desc_crossfit `mname', vname(`var') etype(zeq)
				mat `cfresults_Z`znum'' = r(cfresults)
				local ++znum
			}
		}
	}
	else if `cflag' {
		di
		di as text "No crossfitting results to display."
	}
	
	// stacking results in detail
	if `stflag' & `crossfitted' {	
		di
		di as text "Stacking results (detail):"
		desc_stacking `mname', vname(`nameY') etype(yeq)
		tempname stweights_Y1 ssweights_Y1 psweights_Y1
		mat `stweights_Y1'=r(stweights)
		mat `ssweights_Y1'=r(ssweights)
		mat `psweights_Y1'=r(psweights)
		// should always be a D eqn
		if `numeqnD' {
			local dnum 1
			foreach var in `nameD' {
				tempname stweights_D`dnum' ssweights_D`dnum' psweights_D`dnum'
				desc_stacking `mname', vname(`var') etype(deq)
				mat `stweights_D`dnum''=r(stweights)
				mat `ssweights_D`dnum''=r(ssweights)
				mat `psweights_D`dnum''=r(psweights)
				local ++dnum
			}
		}
		if `numeqnZ' {
			local znum 1
			foreach var in `nameZ' {
				desc_stacking `mname', vname(`var') etype(zeq)
			}
		}
	}
	else if `stflag' {
		di
		di as text "No stacking results to display."
	}	// end stacking block

	// CVC results in detail
	if `cvcflag' & `crossfitted' {
		di
        di as res `"CVC test ({browse "https://doi.org/10.1080/01621459.2019.1672556":Lei 2020}):"'
        di as res "H0: " as text "Learner has the lowest predictive risk among all candidate learners."
        di as res "HA: " as text "There is another learner with lower predictive risk."
        di
        di as res "CVC test p-values:"
		desc_cvc `mname', vname(`nameY') etype(yeq)
		// should always be a D eqn
		if `numeqnD' {
			foreach var in `nameD' {
				desc_cvc `mname', vname(`var') etype(deq)
			}
		}
		if `numeqnZ' {
			foreach var in `nameZ' {
				desc_cvc `mname', vname(`var') etype(zeq)
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
	else if `eflag' {
		di
		di as text "No estimation results to display."
	}
	
	// clear this global from Mata
	mata: mata drop `eqn'

	// return results
	if `cflag' & `crossfitted' {
		// return results in r(.) macros
		return mat cfresults_Y1			= `cfresults_Y1'
		forvalues j=1/`numeqnD' {
			return mat cfresults_D`j'	= `cfresults_D`j''
		}
		forvalues j=1/`numeqnZ' {
			return mat cfresults_Z`j'	= `cfresults_Z`j''
		}
	}
	if `stflag' & `crossfitted' {
		// return results in r(.) macros
		return mat stweights_Y1			= `stweights_Y1'
		return mat ssweights_Y1			= `ssweights_Y1'
		return mat psweights_Y1			= `psweights_Y1'
		forvalues j=1/`numeqnD' {
			return mat stweights_D`j'	= `stweights_D`j''
			return mat ssweights_D`j'	= `ssweights_D`j''
			return mat psweights_D`j'	= `psweights_D`j''
		}
	}
end

prog define desc_stacking, rclass

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

	tempname wmat
	
	di
	di as text "Conditional expectation: " as res "`vname'"
	

	// standard stacking
	local vtilde `vtlist'
	qui _ddml_extract, mname(`mname') show(stweights)
	mat `wmat' = r(`vtilde'_w_mn)
	di
	di as text "Standard stacking weights: " _c
	if `wmat'[1,1] ~= . {
		di as res "`vtilde'"
		display_weights `wmat', dh(`dh') nreps(`nreps')
		di as text "(nb: std stacking weights above are means across `kfolds' folds)"
	}
	else {
		// don't display varname if weights not available
		di as res "(no weights available)"
	}
	return matrix stweights=`wmat'

	if `ssflag' {
		qui _ddml_extract, mname(`mname') show(ssweights)
		mat `wmat' = r(`shortstack'_ss)
		di
		di as text "Short-stacking weights: " as res "`shortstack'_ss"
		display_weights `wmat', dh(`dh') nreps(`nreps')
		return matrix ssweights=`wmat'
	}

	// pooled stacking
	qui _ddml_extract, mname(`mname') show(psweights)
	mat `wmat' = r(`poolstack'_ps)
	di
	di as text "Pooled stacking weights: " _c
	if `wmat'[1,1] ~= . {
		di as res "`poolstack'_ps"
		display_weights `wmat', dh(`dh') nreps(`nreps')
	}
	else {
		// don't display varname if weights not available
		di as res "(no weights available)"
	}
	return matrix psweights=`wmat'


end


prog define desc_cvc, rclass

	syntax name(name=mname), vname(string) etype(string)	// etype is yeq, deq or zeq (not dheq)

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
	
	if "`etype'"=="yeq" {
		local yd y
	}
	else if "`etype'"=="deq" {
		local yd d
	}

	tempname wmat
	
	di
	di as text "Conditional expectation " _c
	if "`yd'"=="y" {
		di as text "E[Y|X]: " _c
	}
	else if "`yd'"=="d" {
		di as text "E[D|Z,X]: " _c
	}
	di as res "`vname'"
	

	if ("`model'"=="interactive" & "`yd'"=="y") 	|	///
		("`model'"=="interactiveiv" & "`yd'"=="y")	|	///
		("`model'"=="interactiveiv" & "`yd'"=="d")		///
		{
		tempname cvc_pval cvc0_pval cvc1_pval
		
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		forvalues m=1/`nreps' {
			_ddml_extract cvc0_pval, mname(`mname') vname(`vname') key1(`vname') key2(cvc0_pval) key3(`m') stata
			mat `cvc0_pval' = nullmat(`cvc0_pval') \ r(cvc0_pval)
			_ddml_extract cvc1_pval, mname(`mname') vname(`vname') key1(`vname') key2(cvc1_pval) key3(`m') stata
			mat `cvc1_pval' = nullmat(`cvc1_pval') \ r(cvc1_pval)
		}
		forvalues i=0/1 {
			mat `cvc`i'_pval' = `cvc`i'_pval''
			mata: st_matrix("`cvc`i'_pval'",(mean(st_matrix("`cvc`i'_pval'")')' , st_matrix("`cvc`i'_pval'")))
			mat `cvc`i'_pval' = J(`nlearners',1,`i') , `cvc`i'_pval'
			mata: st_matrix("`cvc`i'_pval'",(runningsum(J(`nlearners',1,1)), st_matrix("`cvc`i'_pval'")))
			mat rownames `cvc`i'_pval' = `lnames'
		}
		mat `cvc_pval' = `cvc0_pval' \ `cvc1_pval'
		display_weights `cvc_pval', dh(d) nreps(`nreps')
	}
	else {
		tempname cvc_pval
		
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		forvalues m=1/`nreps' {
			_ddml_extract cvc_pval, mname(`mname') vname(`vname') key1(`vname') key2(cvc_pval) key3(`m') stata
			mat `cvc_pval' = nullmat(`cvc_pval') \ r(cvc_pval)
		}
		mat `cvc_pval' = `cvc_pval''
		mata: st_matrix("`cvc_pval'",(runningsum(J(`nlearners',1,1)), mean(st_matrix("`cvc_pval'")')' , st_matrix("`cvc_pval'")))
		mat rownames `cvc_pval' = `lnames'
		display_weights `cvc_pval', nreps(`nreps')
	}
	if "`yd'"=="d" & "`model'"=="fiv" {
		tempname cvc_h_pval
		
		if `pystackedmulti' {
			_ddml_extract lnames, mname(`mname') vname(`vname') key1(`vtlist') key2(stack_base_est_h) stata
			local lnames = r(lnames)
		}
		else {
			local lnames `vtlist'
		}
		local nlearners : list sizeof lnames
		
		forvalues m=1/`nreps' {
			_ddml_extract cvc_h_pval, mname(`mname') vname(`vname') key1(`vname') key2(cvc_h_pval) key3(`m') stata
			mat `cvc_h_pval' = nullmat(`cvc_h_pval') \ r(cvc_h_pval)
		}
		mat `cvc_h_pval' = `cvc_h_pval''
		mata: st_matrix("`cvc_h_pval'",(runningsum(J(`nlearners',1,1)), mean(st_matrix("`cvc_h_pval'")')' , st_matrix("`cvc_h_pval'")))
		mat rownames `cvc_h_pval' = `lnames'
		di as text "Conditional expectation E[D|X]: " as res "`vname'"
		display_weights `cvc_h_pval', nreps(`nreps')
	}

end


prog define display_weights

	syntax name(name=wmat), [ dh(string) nreps(integer 1) ]
	
	local rnames : rownames `wmat'
	local nrows = rowsof(`wmat')

	// check if weights available
	if el(`wmat',1,1)==. {
		di as res "(no weights available)"
		exit
	}
	
	if "`dh'"=="" {
		di as text "Learner" _col(27) "Mean" as res "   |" _c
		local mloc=26
		local dhcol=0
	}
	else if "`dh'"=="d" {
		di as text "Learner" _col(23) "Treated" _col(33) "Mean" as res "   |" _c
		local mloc=32
		local dhcol=1
	}
	else if "`dh'"=="h" {
		di as text "Learner" _col(22) "Cond exp" _col(33) "Mean" as res "   |" _c
		local mloc=32
		local dhcol=1
	}
	forvalues j=1/`nreps' {
		di as text %8s "Rep `j'" _c
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
		di as res _col(`mloc') %5.3f el(`wmat',`i',2+`dhcol') "   |" _c
		forvalues j=1/`nreps' {
			di as res %8.3f el(`wmat',`i',2+`dhcol'+`j') _c
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

	foreach vtilde in `vtlist' {
		di as res _col(2) "Learner:" _col(15) "`vtilde'"
		mata: st_local("estring", return_learner_item(`eqn',"`vtilde'","estring"))
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
			mata: st_local("estring_h", return_learner_item(`eqn',"`vtilde'","estring_h"))
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
	foreach vtilde in `vtlist' {
// ?
		// for pystacked stacking and non-pystacked
		if `crossfitted' /* & `stdflag' */ {
			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse_folds
					mata: st_local("rsq", strofreal(return_result_item(`eqn',"`vtilde'","R-sq","`m'")))
					mata: st_local("mse", strofreal(return_result_item(`eqn',"`vtilde'","MSE","`m'")))
					mata: st_matrix("`mse_folds'", return_result_item(`eqn',"`vtilde'","MSE_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					local lnrname `vtilde'_`m'
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
					mata: st_local("rsq0", strofreal(return_result_item(`eqn',"`vtilde'","R-sq0","`m'")))
					mata: st_local("rsq1", strofreal(return_result_item(`eqn',"`vtilde'","R-sq1","`m'")))
					mata: st_local("mse0", strofreal(return_result_item(`eqn',"`vtilde'","MSE0","`m'")))
					mata: st_local("mse1", strofreal(return_result_item(`eqn',"`vtilde'","MSE1","`m'")))
					mata: st_matrix("`mse0_folds'", return_result_item(`eqn',"`vtilde'","MSE0_folds","`m'"))
					mata: st_matrix("`mse1_folds'", return_result_item(`eqn',"`vtilde'","MSE1_folds","`m'"))
					if `firstrow' {
						di as res "`vnabbrev'" _c
						local firstrow = 0
					}
					forvalues i=0/1 {
						local lnrname `vtilde'`i'_`m'
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
		if `heqn' & `crossfitted' /* & `stdflag' */ {
			forvalues m=1/`nreps' {
				tempname mse_h_folds
				mata: st_local("rsq_h", strofreal(return_result_item(`eqn',"`vtilde'","R-sq_h","`m'")))
				mata: st_local("mse_h", strofreal(return_result_item(`eqn',"`vtilde'","MSE_h","`m'")))
				mata: st_matrix("`mse_h_folds'", return_result_item(`eqn',"`vtilde'","MSE_h_folds","`m'"))
				local lnrname `vtilde'_h_`m'
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


