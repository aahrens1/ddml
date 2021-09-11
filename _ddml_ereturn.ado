

program _ddml_ereturn, eclass

syntax , [mname(name) ]
	
		// ereturn stacking weights
		mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
		
		tempname eqn
		mata: `eqn' = init_eqnStruct()

		forvalues i=1/`numeqns' {
		
			// initialize prior to calling crossfit
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			//mata: st_local("vname",`eqn'.Vname)
			mata: st_local("cmd",`eqn'.command)
			if ("`cmd'"=="pystacked") {
				tempname wmat wmat_h wmat0 wmat1
				mata: st_local("wmat_exists",strofreal(`eqn'.stack_weights!=.))
				mata: st_local("wmat_h_exists",strofreal(`eqn'.stack_weights_h!=.))
				mata: st_local("wmat0_exists",strofreal(`eqn'.stack_weights0!=.))
				mata: st_local("wmat1_exists",strofreal(`eqn'.stack_weights1!=.))
				mata: st_matrix("`wmat'",`eqn'.stack_weights')
				mata: st_matrix("`wmat_h'",`eqn'.stack_weights_h')
				mata: st_matrix("`wmat0'",`eqn'.stack_weights0')
				mata: st_matrix("`wmat1'",`eqn'.stack_weights1')
				if (`wmat_exists'==1) ereturn mat `vtilde'_wmat = `wmat'
				if (`wmat_h_exists'==1) ereturn mat `vtilde'_wmath = `wmat_h'
				if (`wmat0_exists'==1) ereturn mat `vtilde'_wmat0 = `wmat0'
				if (`wmat1_exists'==1) ereturn mat `vtilde'_wmat1 = `wmat1'
			}
		}	 

end


mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}


end