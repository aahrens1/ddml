** last edited: 14 nov 2020

program define _ddml_allcombos, rclass
	version 13
	
	syntax anything , [putlast(string) debug]
	tokenize `anything' , parse("|")

	tempname out
	mata: `out' = get_combos("`anything'")
	mata: `out' = put_last(`out',"`putlast'")
	mata: st_rclear()
	mata: mat_to_string(`out')
	mata: mat_to_colstring(`out')

	if ("`debug'"!="") {
		mata: `out'
	}

	forvalues i = 1(1)`r(k)' {
		return local colstr`i' `r(colstr`i')'
	}
	return local str `r(str)'
	return local k `r(k)'
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
string matrix get_combos(string scalar input)
{

	input = tokens(input,"|")

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
	
	return(out)

}

// matrix to one string where combinations are seperated by "|"
void mat_to_string(string matrix inmat)
{
	r = rows(inmat)
	for (i=1;i<=r;i++) {

		if (i==1) {
			str = invtokens(inmat[i,]) 
		}
		else {
			str = str + " | " + invtokens(inmat[i,]) 
		}
	} 

	st_global("r(str)",str)

}

// matrix to one string per column
void mat_to_colstring(string matrix inmat)
{
	k = cols(inmat)
	st_numscalar("r(k)",k)
	for (j=1;j<=k;j++) {

		str = invtokens(inmat[,j]') 
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
		