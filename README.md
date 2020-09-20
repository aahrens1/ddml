## A simple Stata package for Double Debiased Machine Learning

This package implements the DDML estimator for the partial linear model. 
It's an early incomplete version; so please use with care. 

```
use https://statalasso.github.io/dta/AJR.dta

cap ado uninstall
net install ddml, from(https://raw.githubusercontent.com/aahrens1/ddml/master/)

set seed 123
  
*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init partial

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq lassoy: (lasso2 logpgp95 lat_abst edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres)
ddml yeq rfy: (randomforest logpgp95 lat_abst edes1975 avelf temp* humid* steplow-oilres, type(reg))
ddml yeq pfory: (pyforest logpgp95 lat_abst edes1975 avelf temp* humid* steplow-oilres, type(regress))

* d-equation:
ddml deq lassod: (lasso2 avexpr lat_abst edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres)
ddml deq rigd: (rlasso avexpr lat_abst edes1975 avelf temp* humid* steplow-oilres)

*** cross-fitting and display mean-squared prediction error
ddml crossfit, kfolds(5)  

*** estimation of parameter of interest
ddml estimate, show(all)
```
