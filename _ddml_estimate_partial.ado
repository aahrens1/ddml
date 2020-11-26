*** ddml estimation: partial linear model
program _ddml_estimate_partial, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] , /// 
								[  ///
								ROBust ///
								show(string) /// dertermines which to post
								clear /// deletes all tilde-variables (to be implemented)
								avplot ///
								* ]

	// locals used below
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
   	mata: st_local("Yopt",`mname'.nameYopt)
   	mata: st_local("Dopt",`mname'.nameDopt)

	make_varlists, mname(`mname')

	local ylist		`s(ylist)'
	local Dlist		`s(Dlist)'
	local ncombos	= s(ncombos)
	local tokenlen	= `ncombos'*2 -1

	local j = 1
	forvalues i = 1(2)`tokenlen' {
		tokenize `ylist' , parse("-")
		local y ``i''
		tokenize `Dlist' , parse("-")
		local d ``i''
		add_prefix `y' `d', prefix("`mname'_")
		// do_regress is OLS but with original varnames
		do_regress `s(vnames)' , nocons `robust' yname(`nameY') dnames(`nameD')
		di
		di as res "DML with Y=`y' and D=`d' (N=`e(N)'):"
		ereturn di
        local j= `j'+1
     }

	*** estimate best model
	add_prefix `Yopt' `Dopt', prefix("`mname'_")
	// do_regress is OLS but with original varnames
	do_regress `s(vnames)' , nocons `robust' yname(`Yopt') dnames(`Dopt')

	// plot
	if ("`avplot'"!="") {
	   twoway (scatter `s(vnames)') (lfit `s(vnames)')
	}

	// display
	di
	if `ncombos' > 1 {
		di as res "Optimal model: DML with optimal Y=`Yopt' and optimal D=`Dopt' (N=`e(N)'):"
	}
	else {
		di as res "DML with Y=`Yopt' and D=`Dopt' (N=`N'):"
	}
	ereturn display

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

// does OLS and reports with substitute yname and dnames
program define do_regress, eclass
	syntax anything, [ yname(name) dnames(namelist) * ]

	qui reg `anything' , `options'

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
	ereturn post `b' `V', depname(`yname') obs(`N')

end

program define make_varlists, sclass

	syntax [anything], mname(name)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	// locals used below
	mata: st_local("numD",strofreal(cols(`mname'.nameD)))
 	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))

	mata: st_local("numeqnsD",strofreal(cols(`mname'.eqnlistD)))
	// vn_list has list of all variables, original names
	// vno_list has list of all corresponding orthogonalized variables
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vname",`eqn'.vname)
		mata: st_local("vtilde",`eqn'.vtilde)
		local vn_list `vn_list' `vname'
		local vno_list `vno_list' `vtilde'
	}

	// get list of unique original varnames
	local vn_uniq	: list uniq vn_list

	// loop through unique original varnames and create varlists
	// varlists are stored in Stata macros
	// each of these macros is a list of all orthogonalized vars for a given original var
	// macro orthog_lists is a macro with these macro names in it
	tempname v_list o_list
	mata: `v_list' = tokens("`vn_list'")
	mata: `o_list' = tokens("`vno_list'")

	foreach vn in `vn_uniq' {
		tempname o_index
		mata: `o_index' = `v_list' :== "`vn'"
		// select all orthogonalized variables for unique original varname vn
		mata: st_local("o_sublist", invtokens(select(`o_list',`o_index')))
		tempname t_list
		local `t_list' `o_sublist'
		local orthog_lists `orthog_lists' `t_list'
	}

	// di "numD: `numD'"
	// di "numeqnsD: `numeqnsD'"
	// di "unique original varnames: `vn_uniq'"
	// di "orthog_lists: `orthog_lists'"
	foreach vl in `orthog_lists' {
		// di "`vl' = ``vl''"
		local Dtilde `Dtilde' - ``vl''
	}
	// starts with a -
	// di "Dtilde=`Dtilde'"
	local dpos_end = `numD'+1
	_ddml_allcombos `Ytilde' `Dtilde', dpos_start(2) dpos_end(`dpos_end')

	sreturn local ncombos	`r(ncombos)'
	sreturn local ylist		`r(ystr)'
	sreturn local Dlist		`r(dstr)'

	// clear from Mata
	mata: mata drop `eqn'
	mata: mata drop `v_list' `o_list' `o_index'

end

mata:

// returns an empty eqnStruct
struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
