*! v20dec2021 ereturn after final estimation

program _ddml_ereturn, eclass

syntax , [mname(name) ]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))

	ereturn local model `model'
	ereturn scalar kfolds = `kfolds'
	ereturn scalar nreps = `nreps'

	ereturn_learners `mname', vname("`nameY'") etype(yeq)
	//ereturn_learners `mname', vname("`nameD'") etype(deq)
	//ereturn_learners `mname', vname("`nameY'") etype(yeq)

end

prog define ereturn_learners, eclass

	syntax name(name=mname), vname(string) etype(string) // etype is yeq, deq or zeq (not dheq)
	
	tempname eqn
	mata: `eqn' = init_eStruct()

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("model",`mname'.model)
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	
	// used below to indicate set of crossfitting results to report
	local pairs		= 0
	local heqn		= 0
	if ("`etype'"=="yeq") & ("`model'"=="interactive" | "`model'"=="late") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="late") {
		local pairs	= 1
	}
	if ("`etype'"=="deq") & ("`model'"=="ivhd") {
		// includes both deq and dheq
		local heqn	= 1
	}

	// return MSPE and pystacked weights
	mata: `eqn' = (`mname'.eqnAA).get("`vname'")
	mata: st_local("vtlist",invtokens(`eqn'.vtlist))
	mata: st_local("ssvname",invtokens(`eqn'.shortstack))
	local firstrow = 1
	foreach vtilde in `vtlist' {
			
			mata: st_local("cmd",return_learner_item(`eqn',"`vtilde'","cmd"))
			local get_pysw = "`cmd'"=="pystacked"

			if `pairs'==0 {
				forvalues m=1/`nreps' {
					tempname mse mse_folds mse_tmp mse_folds_tmp pysw pysw_tmp
					if (`m'==1) {
						mata: `mse'=return_result_item(`eqn',"`vtilde'","MSE","`m'")
						mata: `mse_folds'=return_result_item(`eqn',"`vtilde'","MSE_folds","`m'")
						if `get_pysw' {
							mata: `pysw'=mean(return_result_item(`eqn',"`vtilde'","stack_weights","`m'")')
						}
					}
					else {
						mata: `mse_tmp'=return_result_item(`eqn',"`vtilde'","MSE","`m'")
						mata: `mse_folds_tmp'=return_result_item(`eqn',"`vtilde'","MSE_folds","`m'")
						mata: `mse'=(`mse',`mse_tmp')	
						mata: `mse_folds'=(`mse_folds'\`mse_folds_tmp')		
						if `get_pysw' {
							mata: `pysw_tmp'=mean(return_result_item(`eqn',"`vtilde'","stack_weights","`m'")')
							mata: `pysw'=(`pysw'\`pysw_tmp')									
						}
					}
				}
				mata: st_matrix("`mse'",`mse')
				mata: st_matrix("`mse_folds'",`mse_folds')
				ereturn matrix `vtilde'_mse = `mse'
				ereturn matrix `vtilde'_folds_mse = `mse_folds' 
				if `get_pysw' {
					mata: st_matrix("`pysw'",`pysw')
					ereturn matrix `vtilde'_pysw = `pysw'
				}									
			}
			else {
				forvalues j=0(1)1 {
					forvalues m=1/`nreps' {
						tempname mse mse_folds mse_tmp mse_folds_tmp pysw pysw_tmp
						if (`m'==1) {
							mata: `mse'=return_result_item(`eqn',"`vtilde'","MSE`j'","`m'")
							mata: `mse_folds'=return_result_item(`eqn',"`vtilde'","MSE`j'_folds","`m'")
							if `get_pysw' {
								mata: `pysw'=mean(return_result_item(`eqn',"`vtilde'","stack_weights","`m'")')
							}
						}
						else {
							mata: `mse_tmp'=return_result_item(`eqn',"`vtilde'","MSE","`m'")
							mata: `mse_folds_tmp'=return_result_item(`eqn',"`vtilde'","MSE_folds","`m'")
							mata: `mse'=(`mse',`mse_tmp')	
							mata: `mse_folds'=(`mse_folds'\`mse_folds_tmp')		
							if `get_pysw' {
								mata: `pysw_tmp'=mean(return_result_item(`eqn',"`vtilde'","stack_weights","`m'")')
								mata: `pysw'=(`pysw'\`pysw_tmp')									
							}
						}
					}
					mata: st_matrix("`mse'",`mse')
					mata: st_matrix("`mse_folds'",`mse_folds')
					ereturn matrix `vtilde'`j'_mse = `mse'
					ereturn matrix `vtilde'`j'_folds_mse = `mse_folds' 
					if `get_pysw' {
						mata: st_matrix("`pysw'",`pysw')
						ereturn matrix `vtilde'`j'_pysw = `pysw'
					}
				}									
			}
			
	}

	// clear this global from Mata
	mata: mata drop `eqn'
	
end
