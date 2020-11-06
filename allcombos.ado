program define allcombos
	version 13

	args varnames orthog
	
	// get list of unique original varnames
	local vn_uniq	: list uniq varnames
	di "unique original varnames: `vn_uniq'"

	// loop through unique original varnames and create varlists
	tempname v_list o_list vo_list
	mata: `v_list' = tokens("`varnames'")
	mata: `o_list' = tokens("`orthog'")
	mata: `vo_list' = (tokens("`varnames'") \ tokens("`orthog'"))

	foreach vn in `vn_uniq' {
		di "variable `vn':"
		tempname o_index
		mata: `o_index' = `v_list' :== "`vn'"
		mata: st_global("r(o_sublist)", invtokens(select(`o_list',`o_index')))
		local o_sublist `r(o_sublist)'		
		di "o_sublist: `o_sublist'"
		tempname t_list
		local `t_list' `o_sublist'
		local orthog_lists `orthog_lists' `t_list'
	}
	di
	di "orthog_lists: `orthog_lists'"
	foreach vlist in `orthog_lists' {
		di "vlist = `vlist'"
		di "eval vlist = ``vlist''"
	}

	di
	mata: setup_recursion("`vn_uniq'", "`orthog_lists'")
	
	mata: mata drop `v_list' `o_list' `vo_list' `o_index'

end

mata:

void setup_recursion(string scalar vmacro, string scalar omacros)
{
	"vmacro"
	vmacro
	vn_key = tokens(vmacro)'
	"vn_key"
	vn_key
	"omacros"
	omacros
	omacros_t = tokens(omacros)
	olists = J(0,1,"")
	oplists = J(rows(vn_key),1,NULL)
	for (i=1; i<=rows(vn_key); i++) {
		olists = (olists \ st_local(omacros_t[i]))
		orow = tokens(st_local(omacros_t[i]))
		"orow"
		orow
		oplists[i] = &(tokens(st_local(omacros_t[i])))
	}
	"oplists"
	oplists

	do_recursion(oplists)

}
void do_recursion(pointer matrix oplists)
{
// ""
//"entering do_recursion"

//	showlist(oplists)
	singletons = 0
	for (i=1; i<=rows(oplists); i++) {
		singletons = singletons + cols(*(oplists[i]))
	}
	if (singletons==rows(oplists)) {
""
		"end recursion"
		eststring = ""
		for (i=1; i<=rows(oplists); i++) {
			eststring = eststring + " " + *(oplists[i])
		}
		eststring
	}
	
	// loop through
	newoplists = oplists
	breakloop = 0
	for (i=1; i<=rows(oplists); i++) {
		thisrow = *(oplists[i])
		if (cols(thisrow)>1) {
			for (j=1; j<=cols(thisrow); j++) {
				newoplists[i] = &(thisrow[1,j])
				do_recursion(newoplists)
			}
			breakloop = 1
			break
		}
		if (breakloop) break
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

/*
	// loop through unique original varnames
	local orthog_uniq
	foreach vn in `vn_uniq' {
		di "variable `vn':"
		local i	: list posof "`vn'" in varnames
		di "position: `i'"
		local ovn : word `i' of `orthog'
		di "corresponding orthog var: `ovn'"
		local orthog_uniq `orthog_uniq' `ovn'
	}
	di "list of D vars for estimation: `orthog_uniq'"
	
	// remove first duplicate
	local vn_dups	: list dups varnames
	di "all vars with multiple orthogonalisations: `vn_dups'"
	tokenize `vn_dups'
	local first_dup `1'
	di "first duplicate: `first_dup'"
	local i	: list posof "`first_dup'" in varnames
	di "position of first duplicate: `i'"
	local vn_count : list sizeof varnames
	di "number of entries in vnames: `vn_count'"
	tokenize `vnames'
	forvalues j=1/`vn_count' {
		if `j'~=`i' {
			local new_varnames `new_varnames' ``j''
		}
	}
	tokenize `orthog'
	forvalues j=1/`vn_count' {
		if `j'~=`i' {
			local new_orthog `new_orthog' ``j''
		}
	}
	di "new varname without first duplicate: `new_varnames'"
	di "new orthog without first duplicate: `new_orthog'"
*/


// Internal version of matchnames
// Sample syntax:
// matchnames "`varlist'" "`list1'" "`list2'"
// takes list in `varlist', looks up in `list1', returns entries in `list2', called r(names)
program define matchnames, rclass
	version 11.2
	args	varnames namelist1 namelist2

	local k1 : word count `namelist1'
	local k2 : word count `namelist2'

	if `k1' ~= `k2' {
		di as err "namelist error"
		exit 198
	}
	foreach vn in `varnames' {
		local i : list posof `"`vn'"' in namelist1
		if `i' > 0 {
			local newname : word `i' of `namelist2'
		}
		else {
* Keep old name if not found in list
			local newname "`vn'"
		}
		local names "`names' `newname'"
	}
	local names	: list clean names
	return local names "`names'"
end
