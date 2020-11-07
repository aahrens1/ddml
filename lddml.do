* Locals used in whichddml; set when compiled
local stata_version `c(stata_version)'
local born_date `c(born_date)'
local current_date `c(current_date)'

version 13
mata:
mata clear

void whichddml()
{
""
"ddml ver xxx 3nov2020"
"compiled under Stata " + "`stata_version'" + " born " + "`born_date'"
"Mata library for ddml and related programs"
"authors AA/MS"
st_sclear()
st_global("s(stata_born_date)","`born_date'")
st_global("s(stata_version)","`stata_version'")
st_global("s(compiled_date)","`current_date")
}

struct ddmlStruct {
	string scalar		model
	string scalar		nameY
	string colvector	nameYtilde
	string scalar		nameYopt
	string scalar		nameY0opt
	string scalar		nameY1opt
	string colvector	nameD
	string matrix		nameDtilde
	string colvector	nameDopt
	pointer matrix		eqnlistY
	pointer matrix		eqnlistD
	pointer matrix		eqnlistZ
}

struct eqnStruct {
	string scalar		vname
	string scalar		vtilde
	string scalar		eststring
	string scalar		command
	real scalar			MSE
}

mata mlib create lddml, dir(PERSONAL) replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end
