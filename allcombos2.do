
cap program drop allcombos2
program define allcombos2
	version 13
	
	tokenize `0', parse(|)
	
	local i = 1
	while ("``i''"!="") {
		
		if `i'== 1 {
			
			// put first set of variables into mata vector
			mata: out = tokens("``i''")' 
			
		}
		else {
			
			// put next set of variables into mata vector
			mata: x = tokens(ustrtrim("``i''"))'
			
			// save dimensions
			mata: orows = rows(out)
			mata: xrows = rows(x)
			
			// duplicate rows
			mata: out = J(xrows,1,out)
			mata: x = J(orows,1,x)
			
			// sort existing output matrix
			mata: sortindex = range(1,cols(out),1)'
			mata: out = sort(out,sortindex)
			
			// column bind
			mata: out = (out,x)
		}
		
		local i = `i'+2
	}
	

	mata: out

end

allcombos2 y1 y2 | d1 d2 
allcombos2 y1 y2 y3 | d1 d2 d3 | z1 z2
allcombos2 y1 y2 y3 | d1 d2 d3 | z1 z2 | x1 x2 x3


 
