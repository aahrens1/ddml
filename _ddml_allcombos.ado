** last edited: 14 nov 2020

program define _ddml_allcombos, rclass
	version 13
	
	syntax anything , [ putlast(string) ///
						debug ///  
						ypos_start(int 1) ypos_end(int 1) /// position of Y variables
						dpos_start(int 2) dpos_end(int 0) /// position of D variables
						zpos_start(int 0) zpos_end(int 0) /// position of Z variables
						sep(string) ///
						addprefix(string) ///
						]
	if ("`sep'"=="") {
		local sep -
	}

	// why tokenize? not used below
	tokenize "`anything'" , parse("`sep'")

	// obtain all combinations
	tempname out
	mata: st_rclear()
	mata: `out' = get_combos("`anything'","`sep'")
	return scalar ncombos = `r(ncombos)'
	if ("`debug'"!="") {
		mata: "Matrix `out' :"
		mata: `out'
	}

	// put one specific order at the end (indended for optimal model)
	mata: `out' = put_last(`out',"`putlast'")
	if ("`debug'"!="") {
		mata: "Matrix `out' after putlast:"
		mata: `out'
	}

	// save all in one list separated by |
	mata: mat_to_string(`out'[,`ypos_start'..`ypos_end'],"`sep'","`addprefix'")
	return local ystr `r(str)'

	// save D variables in list separated by |
	if (`dpos_start'!=0 & `dpos_end'!=0) {
		mata: mat_to_string(`out'[,`dpos_start'..`dpos_end'],"`sep'","`addprefix'")
		return local dstr `r(str)'
	}
	// save Z variables in list separated by |
	if (`zpos_start'!=0 & `zpos_end'!=0) {
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

// replicate elements of a vector, while mainting order
string matrix Jsort(string matrix mat,
				real scalar rep
)
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
		