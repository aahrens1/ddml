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
"ddml ver xxx 2dec2020"
"compiled under Stata " + "`stata_version'" + " born " + "`born_date'"
"Mata library for ddml and related programs"
"authors AA/MS"
st_sclear()
st_global("s(stata_born_date)","`born_date'")
st_global("s(stata_version)","`stata_version'")
st_global("s(compiled_date)","`current_date")
}

// uniquely identified by vname = dep var in the equation
struct eStruct {
	string scalar		vname			// name of variable to be orthogonalized
	real matrix			vtlist			// list of orthogonalized (learner) variables
	string scalar		shortstack		// name of shortstack variable
	real scalar			nlearners		// number of learners
	real scalar			nreps			// number of resamplings
	real scalar			lieflag			// =1 if LIE spec with two estimation strings etc.

	pointer(class AssociativeArray) scalar	plrnAA		// pointer to AssociativeArray with all learners (keys=vtilde,object)
	pointer(class AssociativeArray) scalar	presAA		// pointer to AssociativeArray with all learner results (keys=vtilde,object,rep)

}

// should perhaps make vname a required argument for this function
struct eStruct init_eStruct()
{
	struct eStruct scalar			d
	class AssociativeArray scalar	A2, A3

	d.vname			= ""
	d.vtlist		= J(1,0,"")
	d.shortstack	= ""
	d.nlearners		= 0
	d.nreps			= 0
	d.lieflag		= 0
		
	A2.reinit("string",2)
	A3.reinit("string",3)
	
	A2.notfound(NULL)
	A3.notfound(NULL)

	d.plrnAA		= &A2
	d.presAA		= &A3
	
	return(d)
}

// clear results from eStruct
void clear_equation_results(struct eStruct e)
{
	class AssociativeArray scalar	A2, A3
	A3.reinit("string",3)
	e.presAA		= &A3
}

struct mStruct {
	string scalar							model			// model; partial, iv, late, etc
	real colvector							id				// id variable (name in Stata will be modelname_id)
	real scalar								nreps			// number of resamplings
	real scalar								kfolds			// number of crossfitting folds
	pointer(class AssociativeArray) scalar	peqnAA			// pointer to AssociativeArray with all equations
	string scalar							nameY			// dependent variable 
	string colvector						nameD			// treatment variable(s)
	string colvector						nameZ			// instrument(s)
	real scalar								ssflag			// flag for shortstacking
	string scalar							strDatavars		// string with expanded names of Stata variables
	real scalar								crossfitted   	// =1 if crossvalidation has been done; 0 if not
}

struct mStruct init_mStruct()
{

	struct mStruct scalar			d
	class AssociativeArray scalar	A
	
	A.reinit("string",1)
	A.notfound(NULL)
	
	d.model			= ""
	d.id			= J(0,1,.)
	d.nreps			= 0
	d.kfolds		= 0
	d.peqnAA		= &A
	d.nameY			= ""
	d.nameD			= J(1,0,"")
	d.nameZ			= J(1,0,"")
	d.ssflag		= 0
	d.strDatavars	= ""
	d.crossfitted	= 0
	return(d)
}

// clear results from all eStructs in mStruct
void clear_model_results(struct mStruct d)
{
	// if not crossfitted, nothing to clear
	if (d.crossfitted==1) {
		struct eStruct scalar			e
		class AssociativeArray scalar	A3
		
		eqnlist = (d.nameY, d.nameD, d.nameZ)
		for (i=1; i<=cols(eqnlist); i++) {
			clear_equation_results((*(d.peqnAA)).get(eqnlist[i]))
		}
		d.crossfitted	= 0
	}
}


/*
// some of the string matrices are actually vectors
struct ddmlStruct {
	string scalar		model			// model; partial, iv, late, etc
	real colvector		id				// id variable (name in Stata will be modelname_id)
	real matrix			idFold			// col 1 = id, col 2 = fold identifier (can probably drop this)
	real matrix			idSample		// col 1 = id, col 2 = sample indicator (can probably drop this)
	real scalar			nreps			// number of resamplings
	real scalar			kfolds			// number of crossfitting folds
	string scalar		strDatavars		// string with expanded names of Stata variables
	string scalar		nameY			// dependent variable 
	string colvector	nameYtilde		// names of orthogonalized variables
	string colvector	nameY0tilde		// names of orthogonalized variables
	string colvector	nameY1tilde		// names of orthogonalized variables
	string colvector	nameYopt 		// name of optimal orthog. Y variable
	string colvector	nameY0opt		// name of optimal orthog. Y variable E[Y|D=0]
	string colvector	nameY1opt		// name of optimal orthog. Y variable E[Y|D=1]
	string colvector	nameD			// name of treatment variable(s)
	string matrix		nameDtilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameD0tilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameD1tilde		// names of orthogonalized treatment variables OR name of optimal instrument
	string matrix		nameDopt		// name of optimal orthog. D variable(s) (partial linear model)
	string matrix		nameD0opt		// name of optimal orthog. D variable(s) E[D|Z=0]
	string matrix		nameD1opt		// name of optimal orthog. D variable(s) E[D|Z=1]
	string scalar		nameDH			// treatment variable 
	string colvector	nameDHtilde		// E[D|X,Z]
	string scalar		nameDHopt 		// E[D|X,Z] optimal 
	string colvector	nameZ			// name of instrument(s)
	string matrix		nameZtilde		// names of orthogonalized instruments
	string matrix		nameZopt		// names of optimal orthog. instruments
	pointer rowvector	eqnlist			// list of all orthog equations
	string rowvector	eqnlistNames	// names of corresponding Vtilde variable
	real scalar 		crossfitted   	// =1 if crossvalidation has been done; 0 if not
}
*/

/*
// to add: boolean to indicate min MSE / optimal orthogonalized var
struct eqnStruct {
	string scalar		eqntype			// yeq, deq, zeq or dheq
	real scalar			interactive 	// 
	string scalar		Vname			//  
	string scalar		Vtilde			// also acts as equation name  
	string scalar		Vtilde0			//  
	string scalar		Vtilde1			//  
	string scalar 		Vtilde_h 		// intended for LIE
	real scalar 		resid 			// 
	real matrix			idVtilde		// col 1 = id, col 2 = orthogonalized
	string scalar		eststring
	string scalar 		eststring_h  	// secondary estimation string indended for LIE
	string scalar		command
	string scalar 		command_h 		// (intended for LIE)
	string scalar 		vtype 			// type of variable that is generated by -predict-
	real colvector		MSE
	real matrix			MSE_folds		// MSE by fold; col=fold, row=resample
	real colvector      MSE_h 			// (intended for LIE)
	real matrix			MSE_h_folds		// (intended for LIE)
	real colvector 		MSE0
	real colvector 		MSE1
	real matrix			MSE0_folds		// MSE by fold; col=fold, row=resample
	real matrix			MSE1_folds		// MSE by fold; col=fold, row=resample
	real colvector		N
	real matrix			N_folds			// sample size by fold; col=fold, row=resample
	real colvector		N_h				// (intended for LIE)
	real matrix         N_h_folds		// (intended for LIE)
	real colvector		N0
	real colvector		N1
	real matrix			N0_folds		// sample size by fold; col=fold, row=resample
	real matrix			N1_folds		// sample size by fold; col=fold, row=resample
	// possibly drop, or make a colvector corresponding to the resampling number
	real scalar 		crossfitted   	// =1 if crossvalidation has been done; 0 if not
	real matrix 		stack_weights   // weights from use of pystacked
	real matrix 		stack_weights_h   // weights from use of pystacked
	real matrix 		stack_weights0   // weights from use of pystacked
	real matrix 		stack_weights1   // weights from use of pystacked
}
*/

/*
struct ddmlStruct init_ddmlStruct()
{
	struct ddmlStruct scalar	d

	d.eqnlist		= J(1,0,NULL)
	d.nameY			= ""
	d.nameYtilde	= J(1,0,"")
	d.nameY0tilde	= J(1,0,"")
	d.nameY1tilde	= J(1,0,"")
	d.nameYopt		= ""
	d.nameY0opt		= ""
	d.nameY1opt		= ""
	d.nameD			= J(1,0,"")
	d.nameDtilde	= J(1,0,"")
	d.nameD0tilde	= J(1,0,"")
	d.nameD1tilde	= J(1,0,"")
	d.nameDopt		= J(1,0,"")
	d.nameD0opt		= J(1,0,"")
	d.nameD1opt		= J(1,0,"")
	d.nameDH		= J(1,0,"")
	d.nameDHtilde	= J(1,0,"")
	d.nameDHopt		= J(1,0,"")
	d.nameZ			= J(1,0,"")
	d.nameZtilde	= J(1,0,"")
	d.nameZopt		= J(1,0,"")
	return(d)
}
*/

struct ddmlStruct use_model(		string scalar fname)
{
	struct ddmlStruct scalar	m
	fh = fopen(fname,"r")
	m = fgetmatrix(fh,1)	// nonzero second argument required for "strict"
	fclose(fh)
	return(m)
}

/*
struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}
*/

// add item about learner to eStruct
void add_learner_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									transmorphic scalar s
									)
{
	(*(e.plrnAA)).put((key1,key2),s)
}

// retrieve item about learner from eStruct
transmorphic return_learner_item(	struct eStruct e,
									string scalar key1,
									string scalar key2
									)
{
	return((*(e.plrnAA)).get((key1,key2)))
}

// add result from learner/resample to eStruct
void add_result_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep,
									s						//  s can be anything
									)
{
	(*(e.presAA)).put((key1,key2,rep),s)
}

// retrieve result from learner/resample from eStruct
transmorphic return_result_item(	struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep
									)
{
	return((*(e.presAA)).get((key1,key2,rep)))
}

/*
void add_eqn(						struct mStruct m,
									real scalar posof,
									string scalar eqntype,
									string scalar vname,
									string scalar vtilde,
									string scalar estcmd,
									string scalar vtype,
									string scalar prefix,
									string scalar estcmd_h
									)
{

	struct eStruct scalar e
	class AssociativeArray scalar A
	pointer(struct eStruct) scalar p

	model = m.model
	
	// posof = 0 means a new variable that isn't in any list
	// posof > 0 means the variable has already been used and hence an eqn created
	//           value of posof is the position in the name list
	
	// if a new variable, update corresponding varlist
	if (posof==0) {
		if (eqntype=="yeq") {
			m.nameY = vname
		}
		else if (eqntype=="deq") {
			m.nameD = (m.nameD, (vname))
		}
		else if (eqntype="zeq") {
			m.nameZ = (m.nameZ, vname)
		}
	}
	
	// neq = cols(m.eqnlist)
	// set default vtilde
	//if (vtilde=="") {
	//	vtilde = vname+"t"+strofreal(neq+1)
	//}
	
	// retrieve equation entry from AA with key=vname
	e	= (*(m.peqnAA)).get(vname)
	// if no such entry is already there, create a new one
	if (e==NULL) {
		e = init_eStruct()
	}

	(*(m.peqnAA)).put(vname,e)
	
	// check
	"nameY and nameD"
	m.nameY
	m.nameD
	"keys in AA"
	(*(m.peqnAA)).keys()
	e	= ((*(m.peqnAA)).get(vname))
	"vtlist:"
	e.vtlist
	// "elist:"
	// e.elist
	
	


	// everything else
	//e.eqntype		= eqntype
	//e.Vname			= vname
	//e.Vtilde		= prefix+vtilde
	//e.eststring		= estcmd
	//e.command		= tokens(estcmd)[1,1]
	//e.vtype		 	= vtype
	//e.crossfitted	= 0
	//e.stack_weights   	= .
	//e.stack_weights_h 	= .
	//e.stack_weights0 	= .
	//e.stack_weights1 	= .

	// look for existing entry
	//newentry		= 1
	//if (neq==0) {
	//	m.eqnlist		= &e
	//	m.eqnlistNames	= vtilde
	//}
	//else {
	//	for (i=1;i<=neq;i++) {
	//		e0 = *(m.eqnlist[i])
	//		if (e0.Vtilde==vtilde) {
	//			// replace
	//			m.eqnlist[i]		= &e
	//			// unnecessary?
	//			m.eqnlistNames[i]	= vtilde
	//			newentry = 0
	//		}
	//	}
	//	if (newentry==1) {
	//		// new entry
	//		m.eqnlist		= (m.eqnlist, &e)
	//		m.eqnlistNames	= (m.eqnlistNames, vtilde)
	//	}
	//}

	// add to appropriate list of tilde variables if a new entry
	//if (newentry) {
	//	if (eqntype=="yeq") {
	//		m.nameYtilde	= (m.nameYtilde, vtilde)
	//		if (interactive==1) {
	//			m.nameY0tilde = (m.nameY0tilde, vtilde0)
	//			m.nameY1tilde = (m.nameY1tilde, vtilde1)
	//		}
	//	}
	//	else if (eqntype=="deq") {
	//		m.nameDtilde	= (m.nameDtilde, vtilde)
	//		if (estcmd_h!="") {
	//			m.nameDHtilde	= (m.nameDHtilde, vtilde_h)
	//		}
	//		if (interactive==1) {
	//			m.nameD0tilde = (m.nameD0tilde, vtilde0)
	//			m.nameD1tilde = (m.nameD1tilde, vtilde1)
	//		}
	//	}
	//	else if (eqntype=="dheq") {
	//		m.nameDHtilde	= (m.nameDHtilde, vtilde)
	//	}
	//	else if (eqntype=="zeq") {
	//		m.nameZtilde	= (m.nameZtilde, vtilde)
	//	}
	//}
	
	//st_global("r(newentry)",strofreal(newentry))
}
*/
/*
void add_eqn(						struct ddmlStruct m,
									string scalar eqntype,
									string scalar vname,
									string scalar vtilde,
									string scalar estcmd,
									string scalar vtype,
									string scalar prefix,
									string scalar estcmd_h,
									string scalar vtilde_h
									)
{
	struct eqnStruct scalar		e, e0

	neq = cols(m.eqnlist)
	model = m.model
	
	// set default vtilde
	if (vtilde=="") {
		vtilde = vname+"t"+strofreal(neq+1)
	}
	
	// specific to interactive equations
	interactive = 0
	if ((model=="interactive")&(eqntype=="yeq")) {
		interactive = 1
	}
	if ((model=="late")&(eqntype=="yeq"|eqntype=="deq")) {
		interactive = 1
	}
	if (interactive == 1) {
		vtilde0 		= vtilde+"0"
		vtilde1 		= vtilde+"1"
		e.Vtilde0		= prefix+vtilde0
		e.Vtilde1		= prefix+vtilde1
	}
	e.interactive 	= interactive

	// specific to dheq (LIE)
	if (eqntype=="dheq") {
		if (vtilde_h=="") {
			vtilde_h = vname+"h"+strofreal(neq+1)
		}
		e.Vtilde_h 		= prefix+vtilde_h
		e.eststring_h 	= estcmd_h
		e.command_h		= tokens(estcmd_h)[1,1]
	}

	// everything else
	e.eqntype		= eqntype
	e.Vname			= vname
	e.Vtilde		= prefix+vtilde
	e.eststring		= estcmd
	e.command		= tokens(estcmd)[1,1]
	e.vtype		 	= vtype
	e.crossfitted	= 0
	e.stack_weights   	= .
	e.stack_weights_h 	= .
	e.stack_weights0 	= .
	e.stack_weights1 	= .

	// look for existing entry
	newentry		= 1
	if (neq==0) {
		m.eqnlist		= &e
		m.eqnlistNames	= vtilde
	}
	else {
		for (i=1;i<=neq;i++) {
			e0 = *(m.eqnlist[i])
			if (e0.Vtilde==vtilde) {
				// replace
				m.eqnlist[i]		= &e
				// unnecessary?
				m.eqnlistNames[i]	= vtilde
				newentry = 0
			}
		}
		if (newentry==1) {
			// new entry
			m.eqnlist		= (m.eqnlist, &e)
			m.eqnlistNames	= (m.eqnlistNames, vtilde)
		}
	}

	// add to appropriate list of tilde variables if a new entry
	if (newentry) {
		if (eqntype=="yeq") {
			m.nameYtilde	= (m.nameYtilde, vtilde)
			if (interactive==1) {
				m.nameY0tilde = (m.nameY0tilde, vtilde0)
				m.nameY1tilde = (m.nameY1tilde, vtilde1)
			}
		}
		else if (eqntype=="deq") {
			m.nameDtilde	= (m.nameDtilde, vtilde)
			if (estcmd_h!="") {
				m.nameDHtilde	= (m.nameDHtilde, vtilde_h)
			}
			if (interactive==1) {
				m.nameD0tilde = (m.nameD0tilde, vtilde0)
				m.nameD1tilde = (m.nameD1tilde, vtilde1)
			}
		}
		else if (eqntype=="dheq") {
			m.nameDHtilde	= (m.nameDHtilde, vtilde)
		}
		else if (eqntype=="zeq") {
			m.nameZtilde	= (m.nameZtilde, vtilde)
		}
	}
	
	st_global("r(newentry)",strofreal(newentry))
}
*/

/*
void add_to_eqn(					struct ddmlStruct m,
									real scalar eqnumber)

{
	pointer(struct eqnStruct) scalar p

	cmd 			= st_global("r(cmd)")
	mse				= st_numscalar("r(mse)")
	mse_folds		= st_matrix("r(mse_folds)")
	n				= st_numscalar("r(N)")
	n_folds			= st_matrix("r(N_folds)")
	p				= m.eqnlist[1,eqnumber]
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= ((*p).N \ n)
	(*p).command	= cmd

	if (cmd == "pystacked") {
		(*p).stack_weights = st_matrix("r(pysw)")		 
	}

	// MSE by fold list should be initialized to void 0-by-k matrix
	// (otherwise concat fails because of conformability)
	(*p).MSE_folds	= ((*p).MSE_folds \ mse_folds)
	(*p).N_folds	= ((*p).N_folds \ n_folds)
	
	// set crossfitted flag = 1
	(*p).crossfitted	= 1

}
*/

/*
void add_to_eqn_h(					struct ddmlStruct m,
									real scalar eqnumber)
{
	pointer(struct eqnStruct) scalar p

	cmd 			= st_global("r(cmd_h)")
	mse_h			= st_numscalar("r(mse_h)")
	mse_h_folds		= st_matrix("r(mse_h_folds)")
	n_h				= st_numscalar("r(N_h)")
	n_h_folds		= st_matrix("r(N_h_folds)")
	p				= m.eqnlist[1,eqnumber]
	(*p).MSE_h		= ((*p).MSE_h \ mse_h)
	(*p).N_h		= ((*p).N_h \ n_h)
	(*p).command_h	= cmd

	if (cmd == "pystacked") {
		(*p).stack_weights_h = st_matrix("r(pysw_h)")		 
	}

	// MSE by fold list should be initialized to void 0-by-k matrix
	// (otherwise concat fails because of conformability)
	(*p).MSE_h_folds= ((*p).MSE_h_folds \ mse_h_folds)
	(*p).N_h_folds	= ((*p).N_h_folds \ n_h_folds)

}
*/

/*
void add_to_eqn01(					struct ddmlStruct m,
									real scalar eqnumber,
									real scalar Z)
{
	pointer(struct eqnStruct) scalar p
	
	Zstr = strofreal(Z)

	cmd 			= st_global("r(cmd)")
	mse				= st_numscalar("r(mse" + Zstr + ")")
	mse_folds		= st_matrix("r(mse" + Zstr + "_folds)")
	n				= st_numscalar("r(N" + Zstr + ")")
	n_folds			= st_matrix("r(N" + Zstr + "_folds)")

	p				= m.eqnlist[1,eqnumber]
	if (Z==0) {
		(*p).MSE0		= ((*p).MSE0 \ mse)
		(*p).N0			= ((*p).N0 \ n)
		// MSE by fold list should be initialized to void 0-by-k matrix
		// (otherwise concat fails because of conformability)
		(*p).MSE0_folds	= ((*p).MSE0_folds \ mse_folds)
		(*p).N0_folds	= ((*p).N0_folds \ n_folds)
		if (cmd == "pystacked") {
			(*p).stack_weights0 = st_matrix("r(pysw0)")		 
		}
	}
	else {
		(*p).MSE1		= ((*p).MSE1 \ mse)
		(*p).N1			= ((*p).N1 \ n)
		// MSE by fold list should be initialized to void 0-by-k matrix
		// (otherwise concat fails because of conformability)
		(*p).MSE1_folds	= ((*p).MSE1_folds \ mse_folds)
		(*p).N1_folds	= ((*p).N1_folds \ n_folds)
		if (cmd == "pystacked") {
			(*p).stack_weights1 = st_matrix("r(pysw1)")		 
		}
	}
	
	// set crossfitted flag = 1
	(*p).crossfitted	= 1

}
*/

mata mlib create lddml, dir(PERSONAL) replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end
