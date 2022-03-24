####################################################
# imports

import glob
import pandas as pd
import numpy as np
import scipy.ndimage as nd

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=10)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=1.0)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')

import matplotlib.pyplot as plt

####################################################
# load data

evasion_type = ['evasion','no_evasion','evasion_pe']
sigma=4

results_k=[]
for (i,s) in zip(range(3),evasion_type):

    if(i==0 or i==1):
        results2=pd.read_csv('../c/output/ss0_%s.csv'%s)
    else:
        results2=pd.read_csv('../c/output/ss0_%s.csv'%evasion_type[i-2])
        
    if(i==0):
        csv_files = [f for f in glob.glob("../c/output/ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion.csv")]
    elif(i==1):
        csv_files = [f for f in glob.glob("../c/output/ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion.csv")]
    elif(i==2):
        csv_files = [f for f in glob.glob("../c/output/ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_pe.csv")]
        
    for f in csv_files:
        tmp=pd.read_csv(f)
        if(np.abs(tmp.tauk[0]-results2.tauk[0])>0.0051):
            results2 = results2.append(tmp,ignore_index=True)

    results2.KtaxRev = results2.KtaxRev*results2.Y/results2.tauk
    results2['tauk_elast'] = (np.log(results2.KtaxRev)-np.log(results2.KtaxRev[0]))/(np.log(1-results2.tauk)-np.log(1-results2.tauk[0]))

    results_k.append(results2)

results_a=[]
for (i,s) in zip(range(3),evasion_type):

    if(i==0 or i==1):
        results2=pd.read_csv('../c/output/ss0_%s.csv'%s)
    else:
        results2=pd.read_csv('../c/output/ss0_%s.csv'%evasion_type[i-2])

    if(i==0):
        csv_files = [f for f in glob.glob("../c/output/ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion.csv")]
    elif(i==1):
        csv_files = [f for f in glob.glob("../c/output/ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion.csv")]
'    elif(i==2):
        csv_files = [f for f in glob.glob("../c/output/ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_pe.csv")]

    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True)

    results2.A_rep = results2.A_rep*results2.Y
    results2['taua_elast'] = (np.log(results2.A_rep)-np.log(results2.A_rep[0]))/(np.log(1-results2.taua)-np.log(1-results2.taua[0]))
        
    results_a.append(results2)

for r in results_k:
    r.tauk=r.tauk-r.tauk[0]
    r.drop(0,inplace=True)
    r.sort_values(by='tauk',ascending=True,inplace=True)
    r.reset_index(drop=True,inplace=True)
    r.tauk=100*r.tauk

for r in results_a:
    r.drop(0,inplace=True)
    r.sort_values(by='taua',ascending=True,inplace=True)
    r.reset_index(drop=True,inplace=True)
    r.taua=100*r.taua


########################################################
colors = ['#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33']
dashes = [(None,None),(6,1),(2,1),(1,0.5)]

fig,axes=plt.subplots(1,2,figsize=(7.5,4),sharex=False,sharey=False)

axes[0].plot(results_k[0]['tauk'],results_k[0]['tauk_elast'],color=colors[0],dashes=dashes[0],label='Baseline (long run)',linewidth=1.5,alpha=0.8)
axes[0].plot(results_k[1]['tauk'],results_k[1]['tauk_elast'],color=colors[1],dashes=dashes[1],label='No evasion (long run)',linewidth=1.5,alpha=0.8)
axes[0].plot(results_k[2]['tauk'],results_k[2]['tauk_elast'],color=colors[2],dashes=dashes[2],label='Baseline (short run)',linewidth=1.5,alpha=0.8)
#axes[0].plot(results_k[3]['tauk'],results_k[3]['tauk_elast'],color=colors[3],dashes=dashes[3],label='No evasion (short run)',linewidth=1.5,alpha=0.8)

axes[0].set_xlabel('Capital income tax (p.p. chg.)')
axes[0].set_title('(a) Reported capital income',y=1.0,size=10)
#axes[0].set_xlim(-15,45)

axes[1].plot(results_a[0]['taua'],results_a[0]['taua_elast'],color=colors[0],dashes=dashes[0],label='Baseline (long run)',linewidth=1.5,alpha=0.8)
axes[1].plot(results_a[1]['taua'],results_a[1]['taua_elast'],color=colors[1],dashes=dashes[1],label='No evasion (long run)',linewidth=1.5,alpha=0.8)
axes[1].plot(results_a[2]['taua'],results_a[2]['taua_elast'],color=colors[2],dashes=dashes[2],label='Baseline (short run)',linewidth=1.5,alpha=0.8)
#axes[1].plot(results_a[3]['taua'],results_a[3]['taua_elast'],color=colors[3],dashes=dashes[3],label='No evasion (short run)',linewidth=1.5,alpha=0.8)
axes[1].set_xlabel('Wealth tax (\%)')
axes[1].set_title('(b) Reported wealth',y=1.0,size=10)
#axes[1].set_xlim(0,6)

axes[1].legend(loc='upper right',fontsize=8)
fig.subplots_adjust(hspace=0.2,wspace=0.3)
plt.savefig('output/fig/evasion_elast.pdf',bbox='tight')

plt.close('all')
