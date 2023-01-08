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

evasion_type = ['evasion','no_evasion','evasion_hi_ec1','evasion_hi_ec2','evasion_hi_ec3','2type_0.350','evasion_nodrev','evasion_chi0']
gama=[0.43172371688765909, 0.42913860042472707,0.43172371688765909,0.43172371688765909,0.43172371688765909,
      0.43172371688765909,0.43172371688765909,0.43172371688765909]
sigma=4
policy = ['warren']

inpath = '../c/output/'

results_tauk=[]
imaxR_tauk = np.zeros(len(evasion_type))
tauk0=np.zeros(len(evasion_type))

for (i,s,g) in zip(range(len(evasion_type)),evasion_type,gama):

    results2=pd.read_csv(inpath+'ss0_%s.csv'%s)

    results2.approval0 = np.nan
    tauk0[i]=results2.tauk.values[0]
    
    csv_files = [f for f in glob.glob(inpath+"ss1_tauk_0.[0-9][0-9][0-9][0-9][0-9][0-9]_%s.csv"%s)]
        
    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True,sort=False)
        
    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2['Kinc'] = results2.KtaxRev/results2.tauk

    results2['tauk_elast'] = (np.log(results2.Kinc)-np.log(results2.Kinc[0]))/(np.log(1-results2.tauk)-np.log(1-results2.tauk[0]))
    

    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

   
    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost


    if(i<5):
        results2['Rev_lost_level'] = results2['Rev_lost'] - results2['detection_revenue']
    else:
        results2['Rev_lost_level'] = results2['Rev_lost']
        
    results2['A_hid_level'] = 100*results2['A_hid']/results2['A']
    
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

    if(i<5):
        results2.p999 = 100*(1.0-results2.p999)
        results2.p999 = results2.p999 - results2.p999[0]
        results2.p90 = results2.p90 - results2.p90[0]
        results2.p999_total = 100*(1.0-results2.p999_total)
        results2.p999_total = results2.p999_total - results2.p999_total[0]
        results2.p90_total = results2.p90_total - results2.p90_total[0]

    
    results2.sort_values(by='tauk',ascending=True,inplace=True)
    results2.reset_index(drop=True,inplace=True)

    imaxR_tauk[i] = results2.lump_sum.idxmax()
    
    results_tauk.append(results2)


########################################################
results_taua=[]
imaxR_taua = np.zeros(len(evasion_type))

for (i,s,g) in zip(range(len(evasion_type)),evasion_type,gama):

    results2=pd.read_csv(inpath+'ss0_%s.csv'%s)

    results2.approval0 = np.nan    
    csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_%s.csv"%s)]

    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True,sort=False)



        
    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost

    if(i<5):
        results2['Rev_lost_level'] = results2['Rev_lost'] - results2['detection_revenue']
    else:
        results2['Rev_lost_level'] = results2['Rev_lost']
        
    results2['A_hid_level'] = 100*results2['A_hid']/results2['A']


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

    if(i<5):
        results2.p999 = 100*(1.0-results2.p999)
        results2.p999 = results2.p999 - results2.p999[0]
        results2.p90 = results2.p90 - results2.p90[0]
        results2.p999_total = 100*(1.0-results2.p999_total)
        results2.p999_total = results2.p999_total - results2.p999_total[0]
        results2.p90_total = results2.p90_total - results2.p90_total[0]


    results2.sort_values(by='taua',ascending=True,inplace=True)
    results2.reset_index(drop=True,inplace=True)

    imaxR_taua[i] = results2.lump_sum.idxmax()

    results_taua.append(results2)

########################################################
# progressive wealth taxes

results_prog=[]
for (i,e,g) in zip(range(len(evasion_type)),evasion_type,gama):
    
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

    if(i<5):
        results2['Rev_lost_level'] = results2['Rev_lost'] - results2['detection_revenue']
    else:
        results2['Rev_lost_level'] = results2['Rev_lost']
        
    results2['A_hid_level'] = 100*results2['A_hid']/results2['A']

    
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

    if(i<5):
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

colors = ['#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628']
dashes = [(None,None),(12,4),(8,3),(4,2),(2,1),(1,0.5),(0.5,0.25)]

fig, AX = plt.subplots(1,2,figsize=(6.5,3),sharex=False,sharey=False)

# --------------------------
# figure 1a: tauk/laffer

axes = AX[0]

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
#axes.set_xlabel('Change in tax rate (p.p.)')

c = 'lump_sum'
axes.plot(results_tauk[0]['tauk'],results_tauk[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_tauk[1]['tauk'],results_tauk[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)
axes.plot(results_tauk[2]['tauk'],results_tauk[2][c],color=colors[2],dashes=dashes[2],label='High transfer cost',linewidth=2,alpha=0.8)
axes.plot(results_tauk[3]['tauk'],results_tauk[3][c],color=colors[3],dashes=dashes[3],label='High detection rate ',linewidth=2,alpha=0.8)
axes.plot(results_tauk[4]['tauk'],results_tauk[4][c],color=colors[4],dashes=dashes[4],label='High penalty',linewidth=2,alpha=0.8)
axes.plot(results_tauk[5]['tauk'],results_tauk[5][c],color=colors[5],dashes=dashes[5],label='Low rep. wealth elast.',linewidth=2,alpha=0.8)
axes.plot(results_tauk[6]['tauk'],results_tauk[6][c],color=colors[6],dashes=dashes[6],label='No detection revenues',linewidth=2,alpha=0.8)
#axes.plot(results_tauk[6]['tauk'],results_tauk[6][c],color=colors[6],dashes=dashes[6],label='No concealed collateral',linewidth=2,alpha=0.8)

axes.legend(loc='lower right',ncol=2,prop={'size':6})
axes.set_xticks([-15,0,15,30,45])
axes.set_xlim(0,45)
axes.set_title('(a) Capital income tax reform',size=8)
axes.set_xlabel('Change in tax rate (p.p.)')


# --------------------------
# figure 1a: tauk/laffer
axes=AX[1]

axes.axvline(0,color='black',linewidth=1,linestyle='-',alpha=0.6)
#axes.set_xlabel('Change in tax rate (p.p.)')

c = 'lump_sum'
axes.plot(results_taua[0]['taua'],results_taua[0][c],color=colors[0],dashes=dashes[0],label='Evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[1]['taua'],results_taua[1][c],color=colors[1],dashes=dashes[1],label='No evasion',linewidth=2,alpha=0.8)
axes.plot(results_taua[2]['taua'],results_taua[2][c],color=colors[2],dashes=dashes[2],label='High transfer cost',linewidth=2,alpha=0.8)
axes.plot(results_taua[3]['taua'],results_taua[3][c],color=colors[3],dashes=dashes[3],label='High detection rate ',linewidth=2,alpha=0.8)
axes.plot(results_taua[4]['taua'],results_taua[4][c],color=colors[4],dashes=dashes[4],label='High penalty',linewidth=2,alpha=0.8)
axes.plot(results_taua[5]['taua'],results_taua[5][c],color=colors[5],dashes=dashes[5],label='Lower rep. wealth elast.',linewidth=2,alpha=0.8)
axes.plot(results_taua[6]['taua'],results_taua[6][c],color=colors[6],dashes=dashes[6],label='No detection revenues',linewidth=2,alpha=0.8)

#axes.legend(loc='lower left',prop={'size':6})
axes.set_xlim(0,8)
axes.set_title('(b) Wealth tax',size=8)
axes.set_xlabel('Tax rate (\%)')

plt.tight_layout(pad=0,h_pad=1,w_pad=0)
#fig.subplots_adjust(hspace=0.26,wspace=0.3)
plt.savefig('output/fig/fig_app_sens_laffer.pdf',bbox='tight',dpi=300)

plt.close('all')

######################################################
# load SR PE results

results_taua_pe=[]

for (i,s,g) in zip(range(len(evasion_type)),evasion_type,gama):

    results2=pd.read_csv(inpath+'ss0_%s.csv'%s)

    results2.approval0 = np.nan    
    csv_files = [f for f in glob.glob(inpath+"ss1_taua_0.[0-9][0-9][0-9][0-9][0-9][0-9]_%s_pe.csv"%s)]

    for f in csv_files:
        results2 = results2.append(pd.read_csv(f),ignore_index=True,sort=False)

    results2.A_rep = results2.A_rep*results2.Y
    results2['taua_elast'] = (np.log(results2.A_rep)-np.log(results2.A_rep[0]))/(np.log(1-results2.taua)-np.log(1-results2.taua[0]))
    results2.taua=100*(results2.taua)

    results2.sort_values(by='taua',ascending=True,inplace=True)
    results2.reset_index(drop=True,inplace=True)

    results_taua_pe.append(results2)

#results_taua_pe[1].drop(4,axis=0).reset_index(drop=True)
    
#########################################

print('Making latex tables')


def fmt_res(file,tmp):
    if(np.isnan(tmp) or abs(tmp)<1.0e-9):
        file.write('&--')
    else:
        file.write('&%0.2f'%tmp)
N=6
file=open('output/tex/table6_sens.tex','w')

file.write('\\footnotesize\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
file.write('\\renewcommand{\\cellalign}{bc}\n')
file.write('\\begin{tabular}{lcccccc}\n')
file.write('\\toprule\n')
file.write('& & & \\multicolumn{3}{c}{Higher evasion costs}\\\\\n')
file.write('\\cmidrule(rl){4-6}\n')
file.write('Outcome & Baseline & \makecell{No\\\\evasion} & \\makecell{High\\\\transfer cost} & \\makecell{High\\\\detection rate} & \\makecell{High\\\\penalty} & \\makecell{Lower rep.\\\\wealth elast.} \\\\\n')
file.write('\\midrule')


file.write('\\multicolumn{7}{l}{{\\textit{(a) Benchmark equilibrium}}}\\\\\n')
file.write('Concealed wealth (\% total)')
for i in range(N):
    tmp = results_taua[i].A_hid_level[0]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Lost revenues (\% GDP)')
for i in range(N):
    tmp = results_taua[i].Rev_lost_level[0]
    fmt_res(file,tmp)

file.write('\\\\\n[1ex]')

file.write('\\multicolumn{7}{l}{{\\textit{(b) Revenue-maximizing capital income tax}}}\\\\\n')

file.write('Change in tax rate (p.p.)')
for i in range(N):
    j = imaxR_tauk[i]
    tmp = results_tauk[i].tauk[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Transfer (\% avg. wage)')
for i in range(N):
    j = imaxR_tauk[i]
    tmp = results_tauk[i].lump_sum[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Concealed wealth (bench. = 1)')
for i in range(N):
    j = imaxR_tauk[i]
    tmp = results_tauk[i].A_hid[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Lost revenues (bench. = 1)')
for i in range(N):
    j = imaxR_tauk[i]
    tmp = results_tauk[i].Rev_lost[j]
    fmt_res(file,tmp)

file.write('\\\\\n[1ex]')
    

file.write('\\multicolumn{7}{l}{{\\textit{(c) Revenue-maximizing wealth tax}}}\\\\\n')

file.write('Tax rate (\%)')
for i in range(N):
    j = imaxR_taua[i]
    tmp = results_taua[i].taua[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Transfer (\% avg. wage)')
for i in range(N):
    j = imaxR_taua[i]
    tmp = results_taua[i].lump_sum[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Concealed wealth (bench. = 1)')
for i in range(N):
    j = imaxR_taua[i]
    tmp = results_taua[i].A_hid[j]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Lost revenues (bench. = 1)')
for i in range(N):
    j = imaxR_taua[i]
    tmp = results_taua[i].Rev_lost[j]
    fmt_res(file,tmp)

file.write('\\\\\n[1ex]')


file.write('\\multicolumn{7}{l}{\\textit{(d) Warren wealth tax}}\\\\\n')

file.write('Transfer (\% avg. wage)')
for i in range(N):
    tmp = results_prog[i].lump_sum[1]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Concealed wealth (bench. = 1)')
for i in range(N):
    tmp = results_prog[i].A_hid[1]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Lost revenues (bench. = 1)')
for i in range(N):
    tmp = results_prog[i].Rev_lost[1]
    fmt_res(file,tmp)

file.write('\\\\\n')
    
file.write('Welfare (\% chg.)')
for i in range(N):
    tmp = results_prog[i].welfare0[1]
    fmt_res(file,tmp)

file.write('\\\\\n')

file.write('Approval (\%)')
for i in range(N):
    tmp = results_prog[i].approval0[1]
    fmt_res(file,tmp)

file.write('\\\\\n[1ex]')

file.write('\\multicolumn{7}{l}{\\textit{(e) Reported wealth elasticity to 1.5\% tax}}\\\\\n')

file.write('Short run')
for i in range(N):
    idx=3
    tmp = results_taua_pe[i].taua_elast[idx] - results_taua_pe[1].taua_elast[idx]
    fmt_res(file,tmp)

file.write('\\\\\n')

file.write('Long run')
for i in range(N):

    idxn = 3
    idx=0
    if(i==0):
        idx = 8
    elif(i==1):
        idx=3
    elif(i==2):
        idx=10
    elif(i==3):
        idx=8
    elif(i==4):
        idx=10
    elif(i==5):
        idx=3
        
    tmp = results_taua[i].taua_elast[idx]
    fmt_res(file,tmp)

file.write('\\\\\n')

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
file.write('\\normalsize\n')
file.close()


