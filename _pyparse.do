program define parse_LassoIC, rclass
	syntax [anything], [criterion(string) ///
					NOCONStant ///
					normalize ///
					max_iter(int 500) ///
					eps(real -1)]
	local optstr 
	
	** criterion
	if "`criterion'"=="bic" {
		local optstr `optstr' 'criterion':'bic', 
	}
	else {
		local optstr `optstr' 'criterion':'aic', 
	}
	** intercept
	if "`noconstant'"!="" {
		local optstr `optstr' 'fit_intercept':False,
	}
	else {
		local optstr `optstr' 'fit_intercept':True,
	}
	** normalize
	if ("`normalize'"!="") {
		local optstr `optstr' 'normalize':True,
	}
	else {
		local optstr `optstr' 'normalize':False,
	}
	** max iterations
	if (`max_iter'>0) {
		local optstr `optstr' 'max_iter':`max_iter',
	}
	** eps
	if (`eps'>0) {
		local optstr `optstr' 'eps':`eps',
	} 
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end 

program define parse_rfReg, rclass
	syntax [anything] , [ ///
					n_estimators(int 100) ///
					criterion(string) ///
					max_depth(int -1) ///  
					min_samples_split(int 2) /// only int supported
					min_samples_leaf(int 1) ///
					min_weight_fraction_leaf(real 0) ///
					///max_features(string) ///
					max_leaf_nodes(int -1) ///
					min_impurity_decrease(real 0) ///
					NOBOOTStrap  ///
					oob_score  ///
					n_jobs(int 1) ///
					random_state(int -1) ///
					warm_start ///
					ccp_alpha(real 0) ///
					max_samples(int -1) ///
					]

	local optstr 
	** number of trees in the forest
	if `n_estimators'>0 {
		local optstr `optstr' 'n_estimators':`n_estimators',
	}
	** criterion
	if "`criterion'"=="mae"|"`criterion'"=="mse" {
		local optstr `optstr' 'criterion':'`criterion'',		
	}
	** max depth
	if `max_depth'>0 {
		local optstr `optstr' 'max_depth':`max_depth',
	} 
	else {
		local optstr `optstr' 'max_depth':None,
	}
	** min sample split
	if `min_samples_split'>0 {
		local optstr `optstr' 'min_samples_split':`min_samples_split',
	}
	** min samples leaf
	if `min_samples_leaf'>0 {
		local optstr `optstr' 'min_samples_leaf':`min_samples_leaf',
	}
	** min weight fraction leaf
	if `min_weight_fraction_leaf'>=0 {
		local optstr `optstr' 'min_weight_fraction_leaf':`min_weight_fraction_leaf',
	}
	** max features
	if "`max_features'"=="auto"|"`max_features'"=="sqrt"|"`max_features'"=="log2" {
		local optstr `optstr' 'max_features':'`max_features'',
	} 
	else {
		local optstr `optstr' 'max_features':'auto',
	}
	** max leaf nodes
	if `max_leaf_nodes'>0 {
		local optstr `optstr' 'max_leaf_nodes':`max_leaf_nodes',
	}	
	else {
		local optstr `optstr' 'max_leaf_nodes':None,
	}
	** min impurity decrease
	if `min_impurity_decrease'>=0 {
		local optstr `optstr' 'min_impurity_decrease':`min_impurity_decrease',
	}
	** bootstrap
	if "`nobootstrap'"!="" {
		local optstr `optstr' 'bootstrap':False,
	}
	else {
		local optstr `optstr' 'bootstrap':True,
	}
	** oob score
	if "`oob_score'"!="" {
		local optstr `optstr' 'oob_score':True,
	}
	else {
		local optstr `optstr' 'oob_score':False,
	}
	** n jobs
	local optstr `optstr' 'n_jobs':1,
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** warm start 
	if "`warm_start'"!="" {
		local optstr `optstr' 'warm_start':True,
	}
	else {
		local optstr `optstr' 'warm_start':False,
	}
	** ccp alpha
	if `ccp_alpha'>=0 {
		local optstr `optstr' 'ccp_alpha':`ccp_alpha',
	}
	** max samples
	if `max_samples'>0 {
		local optstr `optstr' 'max_samples':`max_samples',
	}
	else {
		local optstr `optstr' 'max_samples':None,
	}
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end 

program define parse_gradboostReg, rclass
	syntax [anything] , [ ///
					loss(string) ///
					learning_rate(real 0.1) ///
					n_estimators(int 100) ///
					subsample(real 1) ///
					criterion(string) ///
					min_samples_split(int 2) /// only int supported
					min_samples_leaf(int 1) ///
					min_weight_fraction_leaf(real 0) ///
					max_depth(int 3) ///  
					min_impurity_decrease(real 0) ///
					///init()
					///max_features(string) ///
					alpha(real 0.9) ///
					max_leaf_nodes(int -1) ///
					oob_score  ///
					n_jobs(int 1) ///
					validation_fraction(real 0.1) ///
					n_iter_no_change(int -1) ///
					tol(real 1e-4) ///
					random_state(int -1) ///
					warm_start ///
					ccp_alpha(real 0) ///
					]

	local optstr 
	** loss
	if "`loss'"=="ls"|"`loss'"=="lad"|"`loss'"=="huber"|"`loss'"=="huber" {
		local optstr `optstr' 'loss':'`loss'',
	} 
	else if "`loss'"=="" {
		local optstr `optstr' 'loss':'ls',
	}
	** learning rate
	if `learning_rate'>=0 {
		local optstr `optstr' 'learning_rate':`learning_rate',
	}
	** number of trees in the forest
	if `n_estimators'>0 {
		local optstr `optstr' 'n_estimators':`n_estimators',
	}
	** subsample
	if `subsample'>=0 {
		local optstr `optstr' 'subsample':`subsample',
	}
	** criterion
	if "`criterion'"=="mae"|"`criterion'"=="mse"|"`criterion'"=="friedman_mse" {
		local optstr `optstr' 'criterion':'`criterion'',		
	}
	else if "`criterion'"=="" {
		local optstr `optstr' 'criterion':'friedman_mse',	
	}
	** max depth
	if `max_depth'>0 {
		local optstr `optstr' 'max_depth':`max_depth',
	} 
	else {
		local optstr `optstr' 'max_depth':None,
	}
	** min sample split
	if `min_samples_split'>0 {
		local optstr `optstr' 'min_samples_split':`min_samples_split',
	}
	** min samples leaf
	if `min_samples_leaf'>0 {
		local optstr `optstr' 'min_samples_leaf':`min_samples_leaf',
	}
	** min weight fraction leaf
	if `min_weight_fraction_leaf'>=0 {
		local optstr `optstr' 'min_weight_fraction_leaf':`min_weight_fraction_leaf',
	}
	** max features
	if "`max_features'"=="auto"|"`max_features'"=="sqrt"|"`max_features'"=="log2" {
		local optstr `optstr' 'max_features':'`max_features'',
	} 
	else {
		local optstr `optstr' 'max_features':'auto',
	}
	** max leaf nodes
	if `max_leaf_nodes'>0 {
		local optstr `optstr' 'max_leaf_nodes':`max_leaf_nodes',
	}	
	else {
		local optstr `optstr' 'max_leaf_nodes':None,
	}
	** min impurity decrease
	if `min_impurity_decrease'>=0 {
		local optstr `optstr' 'min_impurity_decrease':`min_impurity_decrease',
	}
	** warm start 
	if "`warm_start'"!="" {
		local optstr `optstr' 'warm_start':True,
	}
	else {
		local optstr `optstr' 'warm_start':False,
	}
	** ccp alpha
	if `ccp_alpha'>=0 {
		local optstr `optstr' 'ccp_alpha':`ccp_alpha',
	}
	** alpha
	if `alpha'>=0 {
		local optstr `optstr' 'alpha':`alpha',
	}
	** init
	local optstr `optstr' 'init':None,
	** validation fraction
	if `validation_fraction'>=0 {
		local optstr `optstr' 'validation_fraction':`validation_fraction',
	}
	** n iter no change
	if `n_iter_no_change'>=0 {
		local optstr `optstr' 'n_iter_no_change':`n_iter_no_change',
	} 
	else {
		local optstr `optstr' 'n_iter_no_change':None,
	}
	** tolerance
	if `tol'>=0 {
		local optstr `optstr' 'tol':`tol',
	}
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end 

program define parse_rfClass, rclass
	syntax [anything] , [ ///
					n_estimators(int 100) ///
					criterion(string) ///
					max_depth(int -1) ///  
					min_samples_split(int 2) /// only int supported
					min_samples_leaf(int 1) ///
					min_weight_fraction_leaf(real 0) ///
					///max_features(string) ///
					max_leaf_nodes(int -1) ///
					min_impurity_decrease(real 0) ///
					NOBOOTStrap  ///
					oob_score  ///
					n_jobs(int 1) ///
					random_state(int -1) ///
					warm_start ///
					ccp_alpha(real 0) ///
					max_samples(int -1) ///
					]

	local optstr 
	** number of trees in the forest
	if `n_estimators'>0 {
		local optstr `optstr' 'n_estimators':`n_estimators',
	}
	** criterion
	if "`criterion'"=="mae"|"`criterion'"=="mse" {
		local optstr `optstr' 'criterion':'`criterion'',		
	}
	** max depth
	if `max_depth'>0 {
		local optstr `optstr' 'max_depth':`max_depth',
	} 
	else {
		local optstr `optstr' 'max_depth':None,
	}
	** min sample split
	if `min_samples_split'>0 {
		local optstr `optstr' 'min_samples_split':`min_samples_split',
	}
	** min samples leaf
	if `min_samples_leaf'>0 {
		local optstr `optstr' 'min_samples_leaf':`min_samples_leaf',
	}
	** min weight fraction leaf
	if `min_weight_fraction_leaf'>=0 {
		local optstr `optstr' 'min_weight_fraction_leaf':`min_weight_fraction_leaf',
	}
	** max features
	if "`max_features'"=="auto"|"`max_features'"=="sqrt"|"`max_features'"=="log2" {
		local optstr `optstr' 'max_features':'`max_features'',
	} 
	else {
		local optstr `optstr' 'max_features':'auto',
	}
	** max leaf nodes
	if `max_leaf_nodes'>0 {
		local optstr `optstr' 'max_leaf_nodes':`max_leaf_nodes',
	}	
	else {
		local optstr `optstr' 'max_leaf_nodes':None,
	}
	** min impurity decrease
	if `min_impurity_decrease'>=0 {
		local optstr `optstr' 'min_impurity_decrease':`min_impurity_decrease',
	}
	** bootstrap
	if "`nobootstrap'"!="" {
		local optstr `optstr' 'bootstrap':False,
	}
	else {
		local optstr `optstr' 'bootstrap':True,
	}
	** oob score
	if "`oob_score'"!="" {
		local optstr `optstr' 'oob_score':True,
	}
	else {
		local optstr `optstr' 'oob_score':False,
	}
	** n jobs
	local optstr `optstr' 'n_jobs':1,
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** warm start 
	if "`warm_start'"!="" {
		local optstr `optstr' 'warm_start':True,
	}
	else {
		local optstr `optstr' 'warm_start':False,
	}
	** ccp alpha
	if `ccp_alpha'>=0 {
		local optstr `optstr' 'ccp_alpha':`ccp_alpha',
	}
	** max samples
	if `max_samples'>0 {
		local optstr `optstr' 'max_samples':`max_samples',
	}
	else {
		local optstr `optstr' 'max_samples':None,
	}
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end 

program define parse_gradboostClass, rclass
	syntax [anything] , [ ///
					loss(string) ///
					learning_rate(real 0.1) ///
					n_estimators(int 100) ///
					subsample(real 1) ///
					criterion(string) ///
					min_samples_split(int 2) /// only int supported
					min_samples_leaf(int 1) ///
					min_weight_fraction_leaf(real 0) ///
					max_depth(int 3) ///  
					min_impurity_decrease(real 0) ///
					///init()
					///max_features(string) ///
					alpha(real 0.9) ///
					max_leaf_nodes(int -1) ///
					oob_score  ///
					n_jobs(int 1) ///
					validation_fraction(real 0.1) ///
					n_iter_no_change(int -1) ///
					tol(real 1e-4) ///
					random_state(int -1) ///
					warm_start ///
					ccp_alpha(real 0) ///
					]

	local optstr 
	** loss
	if "`loss'"=="ls"|"`loss'"=="lad"|"`loss'"=="huber"|"`loss'"=="huber" {
		local optstr `optstr' 'loss':'`loss'',
	} 
	else if "`loss'"=="" {
		local optstr `optstr' 'loss':'ls',
	}
	** learning rate
	if `learning_rate'>=0 {
		local optstr `optstr' 'learning_rate':`learning_rate',
	}
	** number of trees in the forest
	if `n_estimators'>0 {
		local optstr `optstr' 'n_estimators':`n_estimators',
	}
	** subsample
	if `subsample'>=0 {
		local optstr `optstr' 'subsample':`subsample',
	}
	** criterion
	if "`criterion'"=="mae"|"`criterion'"=="mse"|"`criterion'"=="friedman_mse" {
		local optstr `optstr' 'criterion':'`criterion'',		
	}
	else if "`criterion'"=="" {
		local optstr `optstr' 'criterion':'friedman_mse',	
	}
	** max depth
	if `max_depth'>0 {
		local optstr `optstr' 'max_depth':`max_depth',
	} 
	else {
		local optstr `optstr' 'max_depth':None,
	}
	** min sample split
	if `min_samples_split'>0 {
		local optstr `optstr' 'min_samples_split':`min_samples_split',
	}
	** min samples leaf
	if `min_samples_leaf'>0 {
		local optstr `optstr' 'min_samples_leaf':`min_samples_leaf',
	}
	** min weight fraction leaf
	if `min_weight_fraction_leaf'>=0 {
		local optstr `optstr' 'min_weight_fraction_leaf':`min_weight_fraction_leaf',
	}
	** max features
	if "`max_features'"=="auto"|"`max_features'"=="sqrt"|"`max_features'"=="log2" {
		local optstr `optstr' 'max_features':'`max_features'',
	} 
	else {
		local optstr `optstr' 'max_features':'auto',
	}
	** max leaf nodes
	if `max_leaf_nodes'>0 {
		local optstr `optstr' 'max_leaf_nodes':`max_leaf_nodes',
	}	
	else {
		local optstr `optstr' 'max_leaf_nodes':None,
	}
	** min impurity decrease
	if `min_impurity_decrease'>=0 {
		local optstr `optstr' 'min_impurity_decrease':`min_impurity_decrease',
	}
	** warm start 
	if "`warm_start'"!="" {
		local optstr `optstr' 'warm_start':True,
	}
	else {
		local optstr `optstr' 'warm_start':False,
	}
	** ccp alpha
	if `ccp_alpha'>=0 {
		local optstr `optstr' 'ccp_alpha':`ccp_alpha',
	}
	** alpha
	if `alpha'>=0 {
		local optstr `optstr' 'alpha':`alpha',
	}
	** init
	local optstr `optstr' 'init':None,
	** validation fraction
	if `validation_fraction'>=0 {
		local optstr `optstr' 'validation_fraction':`validation_fraction',
	}
	** n iter no change
	if `n_iter_no_change'>=0 {
		local optstr `optstr' 'n_iter_no_change':`n_iter_no_change',
	} 
	else {
		local optstr `optstr' 'n_iter_no_change':None,
	}
	** tolerance
	if `tol'>=0 {
		local optstr `optstr' 'tol':`tol',
	}
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end 

program define parse_SVR, rclass
	syntax [anything] , [ ///
					KERnel(string) ///
					degree(int 3) ///
					GAMma(string) ///
					coef0(real 0) ///
					tol(real 1e-3) ///
					C(real 1) ///
					epsilon(real 0.1) ///
					NOSHRinking ///
					cache_size(real 200) ///
					max_iter(int -1) ///
					]

	local optstr 
	** kernel
	if "`kernel'"=="linear"|"`kernel'"=="poly"|"`kernel'"=="rbf"|"`kernel'"=="sigmoid"|"`kernel'"=="precomputed" {
		local optstr `optstr' 'kernel':`kernel',
	}
	else {
		local optstr `optstr' 'kernel':`rbf',
	}
	** degree
	if `degree'>=1 {
		local optstr `optstr' 'degree':`degree',
	}
	** gamma
	if "`gamma'"=="scale"|"`gamma'"=="auto" {
		local optstr `optstr' 'gamma':`gamma',
	}
	else {
		local optstr `optstr' 'gamma':`scale',
	}
	** coef0
	if `coef0'>0 {
		local optstr `optstr' 'coef0':`coef0',
	}
	** epsilon
	if `epsilon'>0 {
		local optstr `optstr' 'epsilon':`epsilon',
	}
	** cache size
	if `cache_size'>0 {
		local optstr `optstr' 'cache_size':`cache_size',
	}
	** max iter
	if `max_iter'==-1|`max_iter'>0 {
		local optstr `optstr' 'max_iter':`max_iter',
	}
	** tol 
	local optstr `optstr' 'tol':`tol'
	** C 
	local optstr `optstr' 'C':`C'
	** shrinking
	if "`noshrinking'"!="" {
		local optstr `optstr' 'shrinking':False,
	}
	else {
		local optstr `optstr' 'shrinking':True,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end

program define parse_LinearSVR, rclass
	syntax [anything] , [ ///
					tol(real 1e-4) ///
					C(real 1) ///
					epsilon(real 0) ///
					NOSHRinking ///
					cache_size(real 200) ///
					max_iter(int 1000) ///
					loss(string) ///
					NOCONStant /// fit_intercept
					intercept_scaling(real 1) ///
					primal /// dual
					random_state(int -1) ///
					]

	local optstr 
	** loss
	if "`loss'"=="epsilon_insensitive"|"`loss'"=="squared_epsilon_insensitive" {
		local optstr `optstr' 'loss':`loss',
	}
	else {
		local optstr `optstr' 'loss':`epsilon_insensitive',
	}
	** epsilon
	if `epsilon'>0 {
		local optstr `optstr' 'epsilon':`epsilon',
	}
	** cache size
	if `cache_size'>0 {
		local optstr `optstr' 'cache_size':`cache_size',
	}
	** max iter
	if `max_iter'>0 {
		local optstr `optstr' 'max_iter':`max_iter',
	}
	** intercept 
	if "`noconstant'"!="" {
		local optstr `optstr' 'fit_intercept':False,
	}
	else {
		local optstr `optstr' 'fit_intercept':True,
	}
	** dual/primal 
	if "`primal'"!="" {
		local optstr `optstr' 'dual':False,
	}
	else {
		local optstr `optstr' 'dual':True,
	}
	** tol 
	local optstr `optstr' 'tol':`tol'
	** intercept scaling
	local optstr `optstr' 'intercept_scaling':`intercept_scaling'
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end

program define parse_SVC, rclass
	syntax [anything] , [ ///
					KERnel(string) ///
					degree(int 3) ///
					GAMma(string) ///
					coef0(real 0) ///
					tol(real 1e-3) ///
					C(real 1) ///
					epsilon(real 0.1) ///
					NOSHRinking ///
					cache_size(real 200) ///
					max_iter(int -1) ///
					]

	local optstr 
	** kernel
	if "`kernel'"=="linear"|"`kernel'"=="poly"|"`kernel'"=="rbf"|"`kernel'"=="sigmoid"|"`kernel'"=="precomputed" {
		local optstr `optstr' 'kernel':`kernel',
	}
	else {
		local optstr `optstr' 'kernel':`rbf',
	}
	** degree
	if `degree'>=1 {
		local optstr `optstr' 'degree':`degree',
	}
	** gamma
	if "`gamma'"=="scale"|"`gamma'"=="auto" {
		local optstr `optstr' 'gamma':`gamma',
	}
	else {
		local optstr `optstr' 'gamma':`scale',
	}
	** coef0
	if `coef0'>0 {
		local optstr `optstr' 'coef0':`coef0',
	}
	** epsilon
	if `epsilon'>0 {
		local optstr `optstr' 'epsilon':`epsilon',
	}
	** cache size
	if `cache_size'>0 {
		local optstr `optstr' 'cache_size':`cache_size',
	}
	** max iter
	if `max_iter'==-1|`max_iter'>0 {
		local optstr `optstr' 'max_iter':`max_iter',
	}
	** tol 
	local optstr `optstr' 'tol':`tol'
	** C 
	local optstr `optstr' 'C':`C'
	** shrinking
	if "`noshrinking'"!="" {
		local optstr `optstr' 'shrinking':False,
	}
	else {
		local optstr `optstr' 'shrinking':True,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end

program define parse_LinearSVC, rclass
	syntax [anything] , [ ///
					tol(real 1e-4) ///
					C(real 1) ///
					epsilon(real 0) ///
					NOSHRinking ///
					cache_size(real 200) ///
					max_iter(int 1000) ///
					loss(string) ///
					NOCONStant /// fit_intercept
					intercept_scaling(real 1) ///
					primal /// dual
					random_state(int -1) ///
					]

	local optstr 
	** loss
	if "`loss'"=="epsilon_insensitive"|"`loss'"=="squared_epsilon_insensitive" {
		local optstr `optstr' 'loss':`loss',
	}
	else {
		local optstr `optstr' 'loss':`epsilon_insensitive',
	}
	** epsilon
	if `epsilon'>0 {
		local optstr `optstr' 'epsilon':`epsilon',
	}
	** cache size
	if `cache_size'>0 {
		local optstr `optstr' 'cache_size':`cache_size',
	}
	** max iter
	if `max_iter'>0 {
		local optstr `optstr' 'max_iter':`max_iter',
	}
	** intercept 
	if "`noconstant'"!="" {
		local optstr `optstr' 'fit_intercept':False,
	}
	else {
		local optstr `optstr' 'fit_intercept':True,
	}
	** dual/primal 
	if "`primal'"!="" {
		local optstr `optstr' 'dual':False,
	}
	else {
		local optstr `optstr' 'dual':True,
	}
	** tol 
	local optstr `optstr' 'tol':`tol'
	** intercept scaling
	local optstr `optstr' 'intercept_scaling':`intercept_scaling'
	** random state
	if `random_state'>0 {
		local optstr `optstr' 'random_state':`random_state',
	}
	else {
		local optstr `optstr' 'random_state':None,
	}
	** return
	local optstr {`optstr'}
	local optstr = subinstr("`optstr'",",}","}",.)
	local optstr = subinstr("`optstr'"," ","",.)
	return local optstr `optstr'
end
