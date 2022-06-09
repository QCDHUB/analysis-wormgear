import sys,os
import numpy as np
import copy

#--matplotlib
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import pylab as py

#--from scipy stack 
from scipy.integrate import quad, cumtrapz

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier
import lhapdf

import kmeanconf as kc

import warnings
warnings.filterwarnings("ignore")

FLAV=[]
FLAV.append('u')
FLAV.append('d')

cmap = matplotlib.cm.get_cmap('plasma')

def gen_xf(wdir,Q2):
    
    print('\ngenerating g1T at Q2 = %s from %s'%(Q2,wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    if 'g1T' not in conf['steps'][istep]['active distributions']:
        if 'g1T' not in conf['steps'][istep]['passive distributions']:
                print('g1T is not an active or passive distribution')
                return 

    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman
    parman.order=replicas[0]['order'][istep]

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']
    

    g1T=conf['g1T']
    #--setup kinematics
    X=10**np.linspace(-6,-1,200)
    X=np.append(X,np.linspace(0.101,0.99,200))

    f1 = get_f1(X,Q2)
    #--compute XF for all replicas        
    XF={}
    cnt=0
    for par in replicas:
        core.mod_conf(istep,core.get_replicas(wdir)[cnt])   
        cnt+=1
        lprint('%d/%d'%(cnt,len(replicas)))

        parman.set_new_params(par,initial=True)

        for flav in FLAV:
            if flav not in XF:  XF[flav]=[]
            func=lambda i: g1T.get_xF(X[i],Q2,flav,f1[flav][i]) 

            XF[flav].append(np.array([func(i) for i in range(len(X))]))

    print() 
    checkdir('%s/data'%wdir)
    filename='%s/data/g1T-Q2=%3.5f.dat'%(wdir,Q2)

    save({'X':X,'Q2':Q2,'XF':XF},filename)
    print ('Saving data to %s'%filename)

def get_f1(X,Q2):

    f1 = {}

    if 'lhapdf_pdf' in conf: name = conf['lhapdf_pdf']
    else:                    name = 'JAM22-PDF_proton_nlo'

    os.environ['LHAPDF_DATA_PATH'] = '%s/obslib/wormgear/lhapdf'%(os.environ['FITPACK'])
    #--get PDFs, 0 for mean value
    PDF   = lhapdf.mkPDF(name,0)

    f1['u'] = np.array([PDF.xfxQ2(2,X[i],Q2)/X[i] for i in range(len(X))])
    f1['d'] = np.array([PDF.xfxQ2(1,X[i],Q2)/X[i] for i in range(len(X))])

    return f1

def plot_xf(wdir,Q2,mode=0):
    #--mode 0: plot each replica
    #--mode 1: plot average and standard deviation of replicas 

    nrows,ncols=1,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*8,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    hand = {}
    thy  = {}
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    if Q2==None: Q2 = 1.27**2

    scale = classifier.get_scale(wdir)

    filename='%s/data/g1T-Q2=%3.5f.dat'%(wdir,Q2)
    #--load data if it exists
    try:
        data=load(filename)
    #--generate data and then load it if it does not exist
    except:
        gen_xf(wdir,Q2)
        data=load(filename)
        
    replicas=core.get_replicas(wdir)
    cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
    best_cluster=cluster_order[0]

    X=data['X']

    for flav in data['XF']:
        mean = np.mean(data['XF'][flav],axis=0)
        std  = np.std (data['XF'][flav],axis=0)

        if flav=='u':     ax = ax11
        elif flav=='d':   ax = ax12
        else: continue

        #--plot each replica
        if mode==0:
            for i in range(len(data['XF'][flav])):

                thy ,= ax.plot(X,np.array(data['XF'][flav][i]),color=cmap(scale[i]),alpha=0.3)
 
        #--plot average and standard deviation
        if mode==1:
            thy  = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=5)


    for ax in [ax11,ax12]:
        ax.set_xlim(0,1.0)

        ax.tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
        ax.tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
        ax.tick_params(axis='x',    which='major', pad = 8)
        ax.set_xticks([0.2,0.4,0.6,0.8])
        minorLocator = MultipleLocator(0.1)
        ax.xaxis.set_minor_locator(minorLocator)

    ax11.set_ylim(-0.01,0.09)
    ax12.set_ylim(-0.09,0.01)

    ax11.set_yticks([0,0.02,0.04,0.06,0.08])
    ax12.set_yticks([-0.08,-0.06,-0.04,-0.02,0])

    ax11.axhline(0.0,ls='--',color='black',alpha=0.5)
    ax12.axhline(0.0,ls='--',color='black',alpha=0.5)

    ax11.set_xlabel(r'\boldmath$x$',size=40)
    ax12.set_xlabel(r'\boldmath$x$',size=40)   
    ax11.xaxis.set_label_coords(0.92,0.00)
    ax12.xaxis.set_label_coords(0.92,0.00)

    ax11.text(0.80 ,0.60  ,r'\boldmath{$u$}'            , transform=ax11.transAxes,size=40)
    ax12.text(0.80 ,0.20  ,r'\boldmath{$d$}'            , transform=ax12.transAxes,size=40)

    if Q2 == 1.27**2: ax12.text(0.05,0.08,r'$Q^2 = m_c^2$'                                  , transform=ax12.transAxes,size=30)
    else:             ax12.text(0.05,0.08,r'$Q^2 = %s$'%Q2 + ' ' + r'\textrm{GeV}' + r'$^2$', transform=ax12.transAxes,size=30)

    minorLocator = MultipleLocator(0.01)
    ax11.yaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.01)
    ax12.yaxis.set_minor_locator(minorLocator)

    if mode==0:
        sm   = py.cm.ScalarMappable(cmap=cmap)
        sm.set_array([])
        cax = fig.add_axes([0.77,0.95,0.20,0.02])
        cax.tick_params(axis='both',which='both',labelsize=20,direction='in')
        cax.xaxis.set_label_coords(0.65,-1.4)
        cbar = py.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0.2,0.4,0.6,0.8])
        cbar.set_label(r'\boldmath${\rm scaled}~\chi^2_{\rm red}$',size=30)

    ax11.set_ylabel(r'\boldmath$xg_{1T}^{(1)(x)}$',size=30)

    py.tight_layout()
    py.subplots_adjust(hspace = 0, wspace = 0.20)

    filename = '%s/gallery/g1T-Q2=%3.5f'%(wdir,Q2)
    if mode==1: filename += '-bands'

    filename+='.png'

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    py.clf()
    print ('Saving figure to %s'%filename)







