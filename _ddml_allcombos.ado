*! ddml v1.4.2
*! last edited: 8aug2023
*! authors: aa/ms

program define _ddml_allcombos, rclass
	version 16
	
	syntax anything , [ putlast(string) ///
						debug ///  
						ypos(int 1) /// position of Y variables
						DPOS_start(int 2) /// position of D variables
						ZPOS_start(int 0) /// position of Z variables
						sep(string) ///
						addprefix(string) ///
						]
	if ("`sep'"=="") {
		local sep -
	}

	tokenize "`anything'" , parse("`sep'")

	// obtain all combinations
	tempname out
	mata: st_rclear()
	mata: `out' = get_combos("`anything'","`sep'")
	return scalar ncombos = `r(ncombos)'
	local ncols = `r(ncols)'

	// determine end position for D and Z vars
	if `zpos_start'>0 {
		local dpos_end = `zpos_start'-1
		local zpos_end = `ncols'
	}
	else {
		local dpos_end = `ncols'
	}

	// put one specific order at the end (intended for optimal model)
	mata: `out' = put_last(`out',"`putlast'")
	if ("`debug'"!="") {
		di as text "_ddml_all_combos:"
		mata: `out'
	}

	// save all Y in one list separated by `sep'
	mata: mat_to_string(`out'[,`ypos'],"`sep'","`addprefix'")
	return local ystr `r(str)'

	// save D variables in list separated by `sep'
	if (`dpos_start'>0) {
		mata: mat_to_string(`out'[,`dpos_start'..`dpos_end'],"`sep'","`addprefix'")
		return local dstr `r(str)'
	}
	// save Z variables in list separated by `sep'
	if (`zpos_start'>0) {
		mata: mat_to_string(`out'[,`zpos_start'..`zpos_end'],"`sep'","`addprefix'")
		return local zstr `r(str)'
	}

	// one string per column
	mata: mat_to_colstring(`out',"`sep'","`addprefix'")
	return scalar nvars = `r(k)'
	forvalues i = 1(1)`r(k)' {
		return local colstr`i' `r(colstr`i')'
	}

	// clear
	mata: mata drop `out'
end

mata: 

// put one specific combination at the last place 
// intended for optimal combination, which should be estimated at the end
string matrix put_last(string matrix mat,
						string scalar last)
{
	last = tokens(last)

	// check if cols match; otherwise do nothing
	if (cols(last)>0 & cols(last)==cols(mat)) {

		// find match
		is = rowsum(mat :== last) :== cols(mat)

		// check if there is one match; otherwise do nothing
		if (sum(is)==1) {

			// split in two and row bind
			mat0 = select(mat,is:==0)
			mat1 = select(mat,is:==1)
			mat = (mat0\mat1)
		}
	} 

	return(mat)
}

// obtain full matrix of all combinations
string matrix get_combos(string scalar input,string scalar sep)
{

	input = tokens(input,sep)

	for (i=1; i<=cols(input); i=i+2) {

		vars = input[1,i]
		vars = ustrtrim(vars)

		if (i==1) {

			out = tokens(vars)'

		}
		else {

			// put next set of variables into mata vector
			x = tokens(vars)'
				
			// save dimensions
			orows = rows(out)
			xrows = rows(x)
				
			// duplicate rows			
			out = Jsort(out,xrows)
			x = J(orows,1,x)
				
			// column bind
			out = (out,x)
			
		}

	}

	st_numscalar("r(ncombos)",rows(out))
	st_numscalar("r(ncols)",cols(out))

	return(out)

}

// matrix to one string where combinations are seperated by "|"
void mat_to_string(string matrix inmat, string scalar sep,string scalar prefix)
{
	
	if (prefix!="") {
		inmat = prefix :+ inmat 
	}

	r = rows(inmat)
	for (i=1;i<=r;i++) {

		the_row = inmat[i,]

		// put in string
		if (i==1) {
			str = invtokens(the_row) 
		}
		else {
			str = str + " " + sep + " " + invtokens(the_row) 
		}
	} 

	st_global("r(str)",str)

}

// matrix to one string per column
void mat_to_colstring(string matrix inmat,string scalar sep,string scalar prefix)
{
	
	if (prefix!="") {
		inmat = prefix :+ inmat 
	}

	k = cols(inmat)
	st_numscalar("r(k)",k)
	for (j=1;j<=k;j++) {

		str = invtokens(inmat[,j]'," "+sep+" ") 
		st_global("r(colstr"+strofreal(j)+")",str)

	} 

}

// replicate elements of a vector, while maintaining order
string matrix Jsort(string matrix mat,
				real scalar rep)
{
	r = rows(mat)
	for (i=1;i<=r;i++) {
		if (i==1) {
			out = J(rep,1,mat[i,])
		} 
		else {
			out = (out\J(rep,1,mat[i,]))
		}
	}
	return(out)
}

end
		
