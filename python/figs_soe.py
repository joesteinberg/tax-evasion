####################################################
# imports

import glob
import pandas as pd
import numpy as np
import scipy.ndimage as nd
import weighted
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

evasion_type = ['evasion','no_evasion']
gama=[0.43172371688765909, 0.42913860042472707]
sigma=4
policy = ['warren','sanders']

inpath = '../c/output/'

results_tauk=[]
imaxR_tauk = [0,0]
imaxW_tauk = [0,0]
tauk0=[0,0]

for (i,s,g) in zip(range(2),evasion_type,gama):

    results2=pd.read_csv(inpath+'ss0_%s.csv'%s)
    results2.approval0 = np.nan
    tauk0[i]=results2.tauk.values[0]
    
    csv_files=[]
    if(i==0):
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_soe.csv")]
    elif i==1:
        csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion_soe.csv")]
        
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

for (i,s,g) in zip(range(3),evasion_type,gama):

    results2=None
    if(i<2):
        results2=pd.read_csv(inpath+'ss0_%s.csv'%s)
    else:
        results2=pd.read_csv(inpath+'ss0_evasion.csv')

    results2.approval0 = np.nan    
    csv_files=[]
    if(i==0):
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_evasion_soe.csv")]
    elif i==1:
        csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_no_evasion_soe.csv")]
        
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
        try:
            results2 = results2.append(pd.read_csv(inpath+'ss1_%s_%s_soe.csv'%(p,e)),ignore_index=True)
        except:
            results2 = results2.append(pd.read_csv(inpath+'ss0_%s.csv'%e),ignore_index=True)

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
plt.savefig('output/fig/fig_app_tauk_soe.pdf',bbox='tight',dpi=300)

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
axes.set_xlim(0,12)

# --------------------------
# figure 2b: taua/laffer

#axes=AX[0][1]
axes = plt.subplot(gs[0,2:4])

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)

c = 'lump_sum'
axes.plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)

axes.set_title(r'(b) Transfer (\% avg. labor inc.)',size=8)
axes.legend(loc='lower left',prop={'size':6})
axes.set_xlim(0,12)

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
axes.set_xlim(0,12)

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
axes.set_xlim(0,12)

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
axes.set_xlim(0,12)


plt.tight_layout(pad=0,h_pad=1,w_pad=0)
#fig.subplots_adjust(hspace=0.26,wspace=0.3)
plt.savefig('output/fig/fig_app_taua_soe.pdf',bbox='tight',dpi=300)


plt.close('all')
