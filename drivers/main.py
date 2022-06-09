#!/usr/bin/env python
import os,sys
#--set lhapdf data path
version = int(sys.version[0])
os.environ["LHAPDF_DATA_PATH"] = '/work/JAM/ccocuzza/lhapdf/python%s/sets'%version
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import kmeanconf as kc

#--from corelib
from analysis.corelib import core, inspect, predict, classifier, optpriors, jar, mlsamples, summary

#--from qpdlib
from analysis.qpdlib import g1T

#--from obslib
#from analysis.obslib import wormgear

#--from parlib
from analysis.parlib  import params, corr

#--primary working directory
try: wdir=sys.argv[1]
except: wdir = None

Q2 = 4

#g1T.gen_xf(wdir,Q2)
#g1T.plot_xf(wdir,Q2,mode=0)                
#g1T.plot_xf(wdir,Q2,mode=1)                
#sys.exit()

######################
##--Initial Processing
######################

inspect.get_msr_inspected(wdir,limit=None)
predict.get_predictions(wdir,force=False)
classifier.gen_labels(wdir,kc)
jar.gen_jar_file(wdir,kc)
summary.print_summary(wdir,kc)

#classifier.plot_chi2_dist_per_exp(wdir,kc,'dy',20002)
#classifier.print_chi2_per_exp(wdir,kc)

###################
##--Optimize priors
###################

#optpriors.gen_priors(wdir,kc,10)

###################
#--Plot g1T
###################

#--generate and plot data
g1T.gen_xf(wdir,Q2)
g1T.plot_xf(wdir,Q2,mode=0)                
g1T.plot_xf(wdir,Q2,mode=1)                



####################
##--Observable plots
####################

wormgear. plot_obs (wdir,kc)


##---------------------------------------------------------------
##--Parameter distributions
##---------------------------------------------------------------
hist=False

params.plot_params(wdir,'g1T',False)
params.plot_params(wdir,'g1T',True)











