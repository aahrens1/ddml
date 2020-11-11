
cap program drop allcombos2
program define allcombos2
	version 13
	
	tokenize `0', parse(|)
	
	mata: out = get_combos("`0'")
	mata: out

end

mata: 

string matrix get_combos(
				string scalar input
)
{

	input = tokens(input,"|")
	input

	for (i=1; i<=cols(input); i=i+2) {

		vars = input[1,i]
		vars = ustrtrim(vars)
		vars

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
			out = J(xrows,1,out)
			x = J(orows,1,x)
				
			// sort existing output matrix
			sortindex = range(1,cols(out),1)'
			out = sort(out,sortindex)
				
			// column bind
			out = (out,x)
			
		}

	}
	
	out

	return(vars)

}

end
		