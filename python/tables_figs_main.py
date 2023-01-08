####################################################
# imports

import glob
import pandas as pd
import numpy as np
import scipy.ndimage as nd

import matplotlib as mpl
mpl.use('Agg')
mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':8})
mpl.rc('font',size=8)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=2.0)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='pdf')


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

####################################################
# load data

print('Loading aggregate results')

evasion_type = ['evasion','no_evasion','evasion_pe','no_evasion_pe']
gama=[0.43172371688765909, 0.42913860042472707,0.43172371688765909,0.42913860042472707]
sigma=4
policy = ['warren','sanders']

inpath = '../c/output/'

results_tauk=[]
imaxR_tauk = [0,0]
imaxW_tauk = [0,0]
tauk0=[0,0,0,0]

for (i,s,g) in zip(range(4),evasion_type,gama):

    results2=None
    if(i==0 or i==2):
        results2=pd.read_csv(inpath+'ss0_evasion.csv')
    elif(i==1 or i==3):
        results2=pd.read_csv(inpath+'ss0_no_evasion.csv')

    results2.approval0 = np.nan
    tauk0[i]=results2.tauk.values[0]
    
    csv_files=[]
    if(i==0):
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion.csv")]
    elif i==1:
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion.csv")]
    elif i==2:
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_pe.csv")]
    elif i==3:
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion_pe.csv")]

        
    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True,sort=False)
        
    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2['Kinc'] = results2.KtaxRev/results2.tauk

    results2['tauk_elast'] = (np.log(results2.Kinc)-np.log(results2.Kinc[0]))/(np.log(1-results2.tauk)-np.log(1-results2.tauk[0]))
    

    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost
    results2.Rev_lost = (results2.Rev_lost/results2.Rev_lost[0])
    #results2.Rev_lost = 100*(results2.Rev_lost/results2.Rev_lost[0]-1.0)
    #results2.Rev_lost  = 100*results2.Rev_lost/results2.Y
    #results2.Rev_lost = results2.Rev_lost - results2.Rev_lost[0]

    results2.tauk=100*(results2.tauk-tauk0[i])
    results2.taua=100*(results2.taua)
    results2.A = results2.A*results2.Y
    results2.A_rep = results2.A_rep*results2.Y
    results2.A_hid = results2.A_hid*results2.Y
    results2.Y = 100*(results2.Y/results2.Y[0]-1.0)
    results2.Q = 100*(results2.Q/results2.Q[0]-1.0)
    results2.K = 100*(results2.K/results2.K[0]-1.0)
    results2.L = 100*(results2.L/results2.L[0]-1.0)
    results2.W = 100*(results2.W/results2.W[0]-1.0)
    results2.A = 100*(results2.A/results2.A[0]-1.0)
    results2.A_rep = 100*(results2.A_rep/results2.A_rep[0]-1.0)
    #results2.A_hid = 100*(results2.A_hid/results2.A_hid[0]-1.0)
    results2.A_hid = results2.A_hid/results2.A_hid[0]
    results2.P = 100*(results2.P/results2.P[0]-1.0)
    results2.r = 100*(results2.r-results2.r[0])
    #results2.lump_sum = 100*(results2.lump_sum/results2.ylbar)
    results2.welfare = 100*((results2.welfare/results2.welfare[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.welfare0 = 100*((results2.welfare0/results2.welfare0[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.approval0 = 100*results2.approval0
    results2.approval = 100*results2.approval

    results2.p999 = 100*(1.0-results2.p999)
    results2.p999 = results2.p999 - results2.p999[0]
    results2.p90 = results2.p90 - results2.p90[0]
    results2.p999_total = 100*(1.0-results2.p999_total)
    results2.p999_total = results2.p999_total - results2.p999_total[0]
    results2.p90_total = results2.p90_total - results2.p90_total[0]


    
    results2.sort_values(by='tauk',ascending=True,inplace=True)
    results2.reset_index(drop=True,inplace=True)

    if(i<2):
        imaxW_tauk[i] = results2.welfare0.idxmax()
        imaxR_tauk[i] = results2.lump_sum.idxmax()
    
    results_tauk.append(results2)


########################################################
results_taua=[]
imaxR_taua = [0,0]
imaxW_taua = [0,0]

for (i,s,g) in zip(range(4),evasion_type,gama):

    results2=None
    if(i==0 or i==2):
        results2=pd.read_csv(inpath+'ss0_evasion.csv')
    elif(i==1 or i==3):
        results2=pd.read_csv(inpath+'ss0_no_evasion.csv')

    results2.approval0 = np.nan    
    csv_files=[]
    if(i==0):
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion.csv")]
    elif i==1:
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion.csv")]
    elif i==2:
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_pe.csv")]
    elif i==3:
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion_pe.csv")]
        
    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True,sort=False)



        
    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost
    results2.Rev_lost = (results2.Rev_lost/results2.Rev_lost[0])
    #results2.Rev_lost = 100*(results2.Rev_lost/results2.Rev_lost[0]-1.0)
    #results2.Rev_lost  = 100*results2.Rev_lost/results2.Y
    #results2.Rev_lost = results2.Rev_lost - results2.Rev_lost[0]

   
    results2.tauk=100*(results2.tauk-results2.tauk[0])
    results2.A = results2.A*results2.Y
    results2.A_rep = results2.A_rep*results2.Y
    results2['taua_elast'] = (np.log(results2.A_rep)-np.log(results2.A_rep[0]))/(np.log(1-results2.taua)-np.log(1-results2.taua[0]))
    results2.taua=100*(results2.taua)
    results2.A_hid = results2.A_hid*results2.Y
    results2.Y = 100*(results2.Y/results2.Y[0]-1.0)
    results2.Q = 100*(results2.Q/results2.Q[0]-1.0)
    results2.K = 100*(results2.K/results2.K[0]-1.0)
    results2.L = 100*(results2.L/results2.L[0]-1.0)
    results2.W = 100*(results2.W/results2.W[0]-1.0)
    results2.A = 100*(results2.A/results2.A[0]-1.0)
    results2.A_rep = 100*(results2.A_rep/results2.A_rep[0]-1.0)
    #results2.A_hid = 100*(results2.A_hid/results2.A_hid[0]-1.0)
    results2.A_hid = results2.A_hid/results2.A_hid[0]
    results2.P = 100*(results2.P/results2.P[0]-1.0)
    results2.r = 100*(results2.r-results2.r[0])
    #results2.lump_sum = 100*(results2.lump_sum/results2.ylbar)
    #results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])

    results2.welfare = 100*((results2.welfare/results2.welfare[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.welfare0 = 100*((results2.welfare0/results2.welfare0[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.approval0 = 100*results2.approval0
    results2.approval = 100*results2.approval
    
    results2.p999 = 100*(1.0-results2.p999)
    results2.p999 = results2.p999 - results2.p999[0]
    results2.p90 = results2.p90 - results2.p90[0]
    results2.p999_total = 100*(1.0-results2.p999_total)
    results2.p999_total = results2.p999_total - results2.p999_total[0]
    results2.p90_total = results2.p90_total - results2.p90_total[0]


    results2.sort_values(by='taua',ascending=True,inplace=True)
    results2.reset_index(drop=True,inplace=True)

    if(i<2):
        imaxW_taua[i] = results2.welfare0.idxmax()
        imaxR_taua[i] = results2.lump_sum.idxmax()

    results_taua.append(results2)

########################################################
# progressive wealth taxes

evasion_type = ['evasion','no_evasion']
policy = ['warren','sanders']
gama=[0.43172371688765909,0.43622648670727987]
sigma=4

results_prog=[]
for (i,e,g) in zip(range(2),evasion_type,gama):
    
    results2=pd.read_csv(inpath+'ss0_%s.csv'%e)

    for p in policy:
        results2 = results2.append(pd.read_csv(inpath+'ss1_%s_%s.csv'%(p,e)),ignore_index=True)

    csv_files = [f for f in glob.glob((inpath+"optpol/ss_opt_0.*[0-9]_*[0-9]_%s.csv"%e))]
    results3 = None
    for f in csv_files:
        tmp = pd.read_csv(f)
        if results3 is None:
            results3 = tmp
        else:
            results3 = results3.append(tmp,ignore_index=True,sort=False)

    if(results3 is not None):
        results3 = results3[(results3.approval0>=0.45) & (results3.taua>0.001)]
            
    if (results3 is not None and len(results3)>0):
        imax=results3.welfare0.idxmax()
        results3=results3.loc[imax,:]
        results2 = results2.append(results3,ignore_index=True,sort=False)
    else:
        results2 = results2.append(pd.read_csv(inpath+'ss0_%s.csv'%e),ignore_index=True,sort=False)

    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)
    #results2.KtaxRev_lost = results2.KtaxRev_lost*results2.Y
    #results2.WtaxRev_lost = results2.WtaxRev_lost*results2.Y

    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost

    results2.Rev_lost = (results2.Rev_lost/results2.Rev_lost[0])
    #results2.Rev_lost  = 100*results2.Rev_lost/results2.Y
    #results2.Rev_lost = results2.Rev_lost - results2.Rev_lost[0]

    results2.tauk=100*(results2.tauk-tauk0[i])
    results2.taua=100*(results2.taua)
    results2.A = results2.A*results2.Y
    results2.A_rep = results2.A_rep*results2.Y
    results2.A_hid = results2.A_hid*results2.Y
    results2.Y = 100*(results2.Y/results2.Y[0]-1.0)
    results2.Q = 100*(results2.Q/results2.Q[0]-1.0)
    results2.K = 100*(results2.K/results2.K[0]-1.0)
    results2.L = 100*(results2.L/results2.L[0]-1.0)
    results2.W = 100*(results2.W/results2.W[0]-1.0)
    results2.A = 100*(results2.A/results2.A[0]-1.0)
    results2.A_rep = 100*(results2.A_rep/results2.A_rep[0]-1.0)
    #results2.A_hid = 100*(results2.A_hid/results2.A_hid[0]-1.0)
    results2.A_hid = results2.A_hid/results2.A_hid[0]
    results2.P = 100*(results2.P/results2.P[0]-1.0)
    results2.r = 100*(results2.r-results2.r[0])
    #results2.lump_sum = 100*(results2.lump_sum/results2.ylbar)
    results2.welfare = 100*((results2.welfare/results2.welfare[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.welfare0 = 100*((results2.welfare0/results2.welfare0[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.approval0 = 100*results2.approval0
    results2.approval = 100*results2.approval
    
    results2.p999 = 100*(1.0-results2.p999)
    results2.p999 = results2.p999 - results2.p999[0]
    results2.p90 = results2.p90 - results2.p90[0]
    results2.p999_total = 100*(1.0-results2.p999_total)
    results2.p999_total = results2.p999_total - results2.p999_total[0]
    results2.p90_total = results2.p90_total - results2.p90_total[0]
    
    results2.reset_index(drop=True,inplace=True)
        
    results_prog.append(results2)


########################################################
print('Making figures')

colors = ['#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33']
dashes = [(None,None),(6,1),(2,1)]

fig = plt.subplots(figsize=(6.5,3.85),sharex=True,sharey=False)
gs = gridspec.GridSpec(2,6)

# --------------------------
# figure 1a: tauk/tax evasion
#axes=AX[0][0]
axes= plt.subplot(gs[0,0:2])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
#axes.set_xlabel('Change in tax rate (p.p.)')

c1='A_hid'
c2 = 'Rev_lost'
ln1=axes.plot(results_tauk[0]['tauk'],results_tauk[0][c1],color=colors[0],dashes=dashes[0],label='Hidden wealth',linewidth=2,alpha=0.8)
#ax = axes.twinx()
ln2=axes.plot(results_tauk[0]['tauk'],results_tauk[0][c2],color=colors[1],dashes=dashes[1],label='Lost revenues',linewidth=2,alpha=0.8)
#axes.set_yscale('log',basey=10)
lns=ln1+ln2
labs=[l.get_label() for l in lns]
axes.legend(lns,labs,loc='upper left',prop={'size':6})
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(-15,45)
axes.set_title('(a) Tax evasion (benchmark = 1)',size=8)


# --------------------------
# figure 1b: tauk/laffer
#axes=AX[0][1]
axes= plt.subplot(gs[0,2:4])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
#axes.set_xlabel('Change in tax rate (p.p.)')

c = 'lump_sum'
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_tauk[1]['tauk'],results_tauk[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.legend(loc='lower right',prop={'size':6})
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(-15,45)
axes.set_title('(b) Transfer (\% avg. labor inc.)',size=8)


# --------------------------
# figure 1c: tauk/ineq
#axes=AX[0][2]
axes= plt.subplot(gs[0,4:6])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
#axes.set_xlabel('Change in tax rate (p.p.)')

c1 = 'p999'
c2 = 'p999_total'
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c1],color=colors[0],dashes=dashes[0],label='Evasion: reported',linewidth=2,alpha=0.8)
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c2],color=colors[2],dashes=dashes[2],label='Evasion: actual',linewidth=2,alpha=0.8)
axes.plot(results_tauk[1]['tauk'],results_tauk[1][c1],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.legend(loc='upper left',prop={'size':6})
axes.set_title('(c) Top 0.1\% share (p.p. chg.)',size=8)
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(-15,45)


# --------------------------
# figure 1d: tauk/GDP
#axes=AX[1][0]
axes= plt.subplot(gs[1,1:3])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes.set_xlabel('Change in tax rate (p.p.)')

c = 'Y'
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_tauk[1]['tauk'],results_tauk[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.legend(loc='lower left',prop={'size':6})
axes.set_title('(d) GDP (\% chg.)',size=8)
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(-15,45)

# --------------------------
# figure 1e: tauk/welfare
#axes=AX[1][1]
axes= plt.subplot(gs[1,3:5])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes.set_xlabel('Change in tax rate (p.p.)')

c = 'welfare0'
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_tauk[1]['tauk'],results_tauk[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.legend(loc='lower left',prop={'size':6})
axes.set_title('(e) Welfare (\% chg.)',size=8)
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(-15,45)

plt.tight_layout(pad=0,h_pad=1,w_pad=0)
#fig.subplots_adjust(hspace=0.26,wspace=0.3)
plt.savefig('output/fig/fig1_tauk.pdf',bbox='tight',dpi=300)

plt.close('all')

# --------------------------
# --------------------------
# --------------------------

# --------------------------
# figure 2a: taua/tax evasion

fig = plt.subplots(figsize=(6.5,3.85),sharex=True,sharey=False)
gs = gridspec.GridSpec(2,6)


#axes=AX[0][0]
axes = plt.subplot(gs[0,0:2])
axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)

c1='A_hid'
c2 = 'Rev_lost'
ln1=axes.plot(results_taua[0]['taua'],results_taua[0][c1],color=colors[0],dashes=dashes[0],label='Hidden wealth',linewidth=2,alpha=0.8)
#ax = axes.twinx()
ln2=axes.plot(results_taua[0]['taua'],results_taua[0][c2],color=colors[1],dashes=dashes[1],label='Lost revenues',linewidth=2,alpha=0.8)
#axes.set_yscale('log',ba)
axes.set_title('(a) Tax evasion (benchmark = 1)',size=8)

#axes[0].set_title('(a) Tax evasion',y=1.025)
lns=ln1+ln2
labs=[l.get_label() for l in lns]
axes.legend(lns,labs,loc='upper left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_xticks(range(9))

# --------------------------
# figure 2b: taua/laffer

#axes=AX[0][1]
axes = plt.subplot(gs[0,2:4])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)

c = 'lump_sum'
axes.plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.set_ylim(-2,4)
axes.set_title(r'(b) Transfer (\% avg. labor inc.)',size=8)
axes.legend(loc='lower left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_xticks(range(9))

# --------------------------
# figure 2c: taua/ineq
#axes=AX[0][2]
axes = plt.subplot(gs[0,4:6])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)

c1 = 'p999'
c2 = 'p999_total'
axes.plot(results_taua[0]['taua'],results_taua[0][c1],color=colors[0],dashes=dashes[0],label='Evasion: reported',linewidth=2,alpha=0.8)
axes.plot(results_taua[0]['taua'],results_taua[0][c2],color=colors[2],dashes=dashes[2],label='Evasion: actual',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c1],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.set_title(r'(c) Top 0.1\% share (p.p. chg.)',size=8)
axes.legend(loc='upper left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_xticks(range(9))

# -------------------------
# figure 2d: taua/GDP
#axes=AX[1][0]
axes = plt.subplot(gs[1,1:3])
  
axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes.set_xlabel('Tax rate (p.p.)')

c = 'Y'
axes.plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.set_title(r'(d) GDP (\% chg.)',size=8)
axes.legend(loc='lower left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_xticks(range(9))

# --------------------------
# figure 2e: taua/welfare
#axes=AX[1][1]
axes = plt.subplot(gs[1,3:5])
  
axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes.set_xlabel('Tax rate (p.p.)')

c = 'welfare0'
axes.plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.set_title(r'(e) Welfare (\% chg.)',size=8)
axes.legend(loc='lower left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_xticks(range(9))


plt.tight_layout(pad=0,h_pad=1,w_pad=0)
#fig.subplots_adjust(hspace=0.26,wspace=0.3)
plt.savefig('output/fig/fig2_taua.pdf',bbox='tight',dpi=300)


plt.close('all')

######################################################

print('Making latex tables')

for i in range(2):
    results_prog[i].abar = 50*results_prog[i].abar/results_prog[i].abar[1]
    results_prog[i].abar[2] = 32

def fmt_res(file,tmp):
    if(np.isnan(tmp) or abs(tmp)<1.0e-9):
        file.write('&--')
    else:
        file.write('&%0.2f'%tmp)

vars = ['tau','abar','lump_sum','A_hid','Rev_lost','Y','W','p999','p999_total','welfare0','approval0']
labs = ['Tax rate (p.p. chg.)','Tax threshold (\$M)','Transfer (\% avg. wage)','Concealed wealth (benchmark = 1)',
        'Lost revenues (benchmark = 1)','GDP (\% chg.)','Wage (\% chg.)','Top 0.1\% share, reported (p.p. chg.)','Top 0.1\% share, actual (p.p. chg.)',
        'Welfare (\% chg.)','Approval rate (\%)']
flags = [1,1,1,0,0,1,1,0,1,1,1]
breaks = [0,0,1,0,1,0,1,0,1,0,0]

panels = ['(a) Baseline model','(b) No-evasion counterfactual']


file=open('output/tex/table4_prog_taua.tex','w')

#file.write('\\begin{table}[p]\n')
file.write('\\footnotesize\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
#file.write('\\begin{center}\n')
#file.write('\\caption{Effects taxing capital income and wealth}\n')
#file.write('\\label{tab:agg_results}\n')

file.write('\\begin{tabular}{lccc}\n')
file.write('\\toprule\n')
#file.write('&\\multicolumn{2}{c}{Revenue-maximizing} & \\multicolumn{3}{c}{Progressive wealth taxes}\\\\\n')
#file.write('\\cmidrule(rl){2-3}\\cmidrule(rl){4-6}\n')
#file.write('Outcome & Capital income tax & Flat wealth tax & Warren & Sanders & Optimal\\\\\n')
file.write('Outcome & Warren & Sanders & Optimal\\\\\n')

for i in range(len(panels)):

    file.write('\\midrule\n')
    file.write('\\multicolumn{4}{l}{\\textit{%s}}\\\\\n'%panels[i])
    
    for j in range(len(vars)):

        v=vars[j]
        vk=v
        va=v
        if(j==0):
            vk=v+'k'
            va=v+'a'
            
        l=labs[j]
        if(i==1 and v == 'p999_actual'):
            l='Top 0.1\% wealth share'

        if(i==0 or (i==1 and flags[j]==1)):
            file.write(l)
            
            #tmp = results_tauk[i][vk].values[imaxR_tauk[i]]
            #fmt_res(file,tmp)

            #tmp = results_tauk[i][vk].values[imaxW_tauk[i]]
            #fmt_res(file,tmp)

            #tmp = results_taua[i][va].values[imaxR_taua[i]]
            #fmt_res(file,tmp)

            #tmp = results_taua[i][va].values[imaxW_taua[i]]
            #if(i==0 and v=='Rev_lost'):
            #    tmp=0.0
            #fmt_res(file,tmp)

            tmp = results_prog[i][va].values[1]
            fmt_res(file,tmp)

            tmp = results_prog[i][va].values[2]
            fmt_res(file,tmp)

            if(v=='Rev_lost' or v=='A_hid'):
                file.write('&--')
            else:
                tmp = results_prog[i][va].values[3]
                fmt_res(file,tmp)

            #            if(i==0 and j==(len(vars)-1)):
            #                file.write('\\\\\n[1.5ex]')
            if(breaks[j]==1):
                file.write('\\\\\n[1.0ex]')
            else:
                file.write('\\\\\n')
                

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
#file.write('\\end{center}\n')
file.write('\\normalsize\n')
#file.write('\\end{table}\n')
file.close()

    
######################################################

results_taua[3] = results_taua[3].drop(4,axis=0).reset_index(drop=True)

print('Evasion elasticities:')
print('\tCapital income:')
print('\t\tEvasion SR: [%0.3f,%0.3f]' % (results_tauk[2].tauk_elast.min(),results_tauk[2].tauk_elast.max()))
print('\t\tEvasion LR: [%0.3f,%0.3f]' % (results_tauk[0].tauk_elast.min(),results_tauk[0].tauk_elast.max()))
#print('\t\tNo evasion SR: [%0.3f,%0.3f]' % (results_tauk[3].tauk_elast.min(),results_tauk[3].tauk_elast.max()))
print('\t\tNo evasion LR: [%0.3f,%0.3f]' % (results_tauk[1].tauk_elast.min(),results_tauk[1].tauk_elast.max()))
print('\tWealth:')
print('\t\tEvasion SR total: [%0.3f,%0.3f]' % ( (results_taua[2].taua_elast).min(),(results_taua[2].taua_elast).max()))
print('\t\tEvasion SR evasion: [%0.3f,%0.3f]' % ( (results_taua[2].taua_elast-results_taua[3].taua_elast).min(),(results_taua[2].taua_elast-results_taua[3].taua_elast).max()))
print('\t\tEvasion LR: [%0.3f,%0.3f]' % (results_taua[0].taua_elast.min(),results_taua[0].taua_elast.max()))
print('\t\tNo evasion SR: [%0.3f,%0.3f]' % (results_taua[3].taua_elast.min(),results_taua[3].taua_elast.max()))
print('\t\tNo evasion LR: [%0.3f,%0.3f]' % (results_taua[1].taua_elast.min(),results_taua[1].taua_elast.max()))


# -----------------------------------
# -----------------------------------
# -----------------------------------
# extras

########################################################

colors = ['#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33']
dashes = [(None,None),(6,1),(2,1),(1,0.5)]

fig,axes=plt.subplots(1,2,figsize=(6.5,3),sharex=False,sharey=False)

axes[0].plot(results_tauk[0]['tauk'],results_tauk[0]['tauk_elast'],color=colors[0],dashes=dashes[0],label='Baseline (long run)',linewidth=1.5,alpha=0.8)
axes[0].plot(results_tauk[1]['tauk'],results_tauk[1]['tauk_elast'],color=colors[1],dashes=dashes[1],label='No evasion (long run)',linewidth=1.5,alpha=0.8)
axes[0].plot(results_tauk[2]['tauk'],results_tauk[2]['tauk_elast'],color=colors[2],dashes=dashes[2],label='Baseline (short run)',linewidth=1.5,alpha=0.8)

axes[0].set_xlabel('Capital income tax (p.p. chg.)')
axes[0].set_title('(a) Reported capital income',y=1.0,size=8)
axes[0].set_xlim(-15,45)

axes[1].plot(results_taua[0]['taua'],results_taua[0]['taua_elast'],color=colors[0],dashes=dashes[0],label='Baseline (long run)',linewidth=1.5,alpha=0.8)
axes[1].plot(results_taua[1]['taua'],results_taua[1]['taua_elast'],color=colors[1],dashes=dashes[1],label='No evasion (long run)',linewidth=1.5,alpha=0.8)

x=results_taua[2]['taua']
y1=results_taua[2]['taua_elast']
y2=results_taua[3]['taua_elast']
axes[1].plot(x,y1-y2,color=colors[2],dashes=dashes[2],label='Baseline (short run)',linewidth=1.5,alpha=0.8)
#axes[1].plot(results_taua[3]['taua'],results_taua[3]['taua_elast'],color=colors[3],dashes=dashes[3],label='No evasion (short run)',linewidth=1.5,alpha=0.8)
axes[1].set_xlabel('Wealth tax (\%)')
axes[1].set_title('(b) Reported wealth',y=1.0,size=8)
axes[1].set_xlim(results_taua[3].taua[1],6)

axes[1].legend(loc='upper right',prop={'size':6})
fig.subplots_adjust(hspace=0.2,wspace=0.15)
plt.savefig('output/fig/evasion_elast.pdf',bbox='tight')

plt.close('all')

########################################################
c='lump_sum'

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':10})
mpl.rc('font',size=18)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=3.0)
mpl.rc('savefig',bbox='tight')
mpl.rc('savefig',format='png')

fig,axes=plt.subplots(1,2,figsize=(11,5),sharex=False,sharey=True)
        
axes[0].axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes[0].axhline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)      
axes[0].plot(results_tauk[0]['tauk'],results_tauk[0][c],color=colors[0],dashes=dashes[0],label='With evasion',linewidth=3,alpha=0.8)
axes[0].plot(results_tauk[1]['tauk'],results_tauk[1][c],color=colors[1],dashes=dashes[1],label='Ignoring evasion',linewidth=3,alpha=0.8)
axes[0].set_xlabel(r'Change in tax rate (p.p.)')
axes[0].set_ylabel(r'Change in tax revenues (pct. GDP)')
axes[0].set_title(r'(a) Capital income tax reform',y=1.025)
axes[0].set_xticks([-15,0,15,30,45])
axes[0].set_xlim(-15,45)

axes[1].axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
axes[1].axhline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)      
axes[1].plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='With evasion',linewidth=3,alpha=0.8)
axes[1].plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='Ignoring evasion',linewidth=3,alpha=0.8)
axes[1].set_xlabel(r'Tax rate (p.p. chg.)')
axes[1].set_title(r'(b) Introducing a wealth tax',y=1.025)
#axes[0].set_ylabel(r'Lump-sum transfer (pct. avg. labor income)')
axes[1].set_xlim(0,8)
axes[1].set_xticks(range(9))

axes[0].legend(loc='lower right',fontsize=14)

fig.subplots_adjust(hspace=0.2,wspace=0.15)
plt.savefig('output/fig/webfig_evasion.png',bbox='tight',dpi=300)

plt.close('all')
