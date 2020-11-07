program define allcombos
	version 13
	
	// typical call: allcombos myest.eqnlistD
	// and macro meqnlist has "myest.eqnlistD" in it
	// note that myest is a global Mata object
	args meqnlist

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_numscalar("r(numeqns)",cols(`meqnlist'))
	local numeqnsD	= `r(numeqns)'
	// vn_list has list of all variables, original names
	// vno_list has list of all corresponding orthogonalized variables
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`meqnlist'[1,`i'])
		mata: st_local("vname",`eqn'.vname)
		mata: st_local("vtilde",`eqn'.vtilde)
		local vn_list `vn_list' `vname'
		local vno_list `vno_list' `vtilde'
	}
	// clear this global from Mata
	mata: mata drop `eqn'

	// get list of unique original varnames
	local vn_uniq	: list uniq vn_list
	// di "unique original varnames: `vn_uniq'"

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

	// di "orthog_lists: `orthog_lists'"

	mata: setup_recursion("`orthog_lists'")
	
	mata: mata drop `v_list' `o_list' `o_index'

end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void setup_recursion(string scalar omacros)
{

	omacros_t = tokens(omacros)
	olists = J(0,1,"")
	oplists = J(cols(omacros_t),1,NULL)
	for (i=1; i<=cols(omacros_t); i++) {
		olists = (olists \ st_local(omacros_t[i]))
		orow = tokens(st_local(omacros_t[i]))
		oplists[i] = &(tokens(st_local(omacros_t[i])))
	}
	
	// "oplists"
	// oplists
	// showlist(oplists)

	st_global("r(oplists)",do_recursion(oplists))

}
string scalar do_recursion(pointer matrix oplists)
{

	singletons = 0
	for (i=1; i<=rows(oplists); i++) {
		singletons = singletons + cols(*(oplists[i]))
	}
	if (singletons==rows(oplists)) {

		eststring = ""
		for (i=1; i<=rows(oplists); i++) {
			eststring = eststring + " " + *(oplists[i])
		}
		return(eststring)
	}
	else {
	
		// loop through
		newoplists = oplists
		breakloop = 0
		allstring = ""
		for (i=1; i<=rows(oplists); i++) {
			thisrow = *(oplists[i])
			if (cols(thisrow)>1) {
				for (j=1; j<=cols(thisrow); j++) {
					newoplists[i] = &(thisrow[1,j])
					allstring = allstring + do_recursion(newoplists)
					if (j<cols(thisrow)) {
						allstring = allstring + ", "
					}
				}
				breakloop = 1
				break
			}
			if (breakloop) break
		}
	
		return(allstring)
	}
}

void showlist(pointer matrix oplists)
{
""
	"contents"
	for (i=1; i<=rows(oplists); i++) {
		*(oplists[i])
	}

}

end

