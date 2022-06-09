import os
conf={}

#--fitting setups
conf['bootstrap']=True
conf['flat par']=True

#--setups for DGLAP
conf['dglap mode']='truncated'
conf['alphaSmode']='backward'
conf['order'] = 'NLO'
conf['Q20']   = 1.27**2

#--datasets

conf['datasets']={}

#--lepton-hadron reactions

qT_over_Q_cut=0.5

##--IDIS
conf['datasets']['wormgear']={}
conf['datasets']['wormgear']['filters']=[]
conf['datasets']['wormgear']['filters'].append("qT/Q>%f"%qT_over_Q_cut)
conf['datasets']['wormgear']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wormgear']['xlsx'][10000]='wormgear/expdata/10000.xlsx' # HERMES   | proton | pi0
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wormgear']['norm']={}


#--parameters
conf['params']={}

#--pdf parameters
conf['params']['g1T']={}

conf['params']['g1T']['u1 N']   ={'value':    3.47549e-01   , 'min':   0.0, 'max':    10, 'fixed': False}
conf['params']['g1T']['u1 a']   ={'value':   -1.21835956e-01, 'min':  -0.9, 'max':     1, 'fixed': False}
conf['params']['g1T']['u1 b']   ={'value':    3.20766744e+00, 'min':     0, 'max':    10, 'fixed': False}
#conf['params']['g1T']['u1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
#conf['params']['g1T']['u1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['g1T']['d1 N']   ={'value':    3.47549e-01   , 'min':   -10, 'max':     0, 'fixed': False}
conf['params']['g1T']['d1 a']   ={'value':   -1.21835956e-01, 'min':  -0.9, 'max':     1, 'fixed': False}
conf['params']['g1T']['d1 b']   ={'value':    3.20766744e+00, 'min':     0, 'max':    10, 'fixed': False}
#conf['params']['g1T']['d1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
#conf['params']['g1T']['d1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}


#--steps
conf['steps']={}

istep=1
#--start: 
conf['ftol']=1e-6
conf['steps'][istep]={}
conf['steps'][istep]['dep']=[]
conf['steps'][istep]['active distributions']=['g1T']
conf['steps'][istep]['passive distributions']=[]
#------------------------------------------------------------------------------------------------------------------
conf['steps'][istep]['datasets']={}
conf['steps'][istep]['datasets']['wormgear']=[]
conf['steps'][istep]['datasets']['wormgear'].append(10000) # HERMES   | proton | pi0
