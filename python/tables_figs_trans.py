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

gama=[0.43172371688765909,0.43622648670727987]

####################################################
# load data

def load_data(i,j):
    
    e = evasion_type[i]
    s = scenario[j]
    g = gama[i]

    results2 = None
    try:
        results2 = pd.read_csv('../c/output/trans_%s_%s.csv'%(s,e))
    except:
        return None
        
    results2.lump_sum = 100*(results2.lump_sum*60/results2.Y[0])
    results2.KtaxRev = results2.KtaxRev*results2.Y
    results2['Kinc'] = results2.KtaxRev/results2.tauk

    results2['tauk_elast'] = (np.log(results2.Kinc)-np.log(results2.Kinc[0]))/(np.log(1-results2.tauk)-np.log(1-results2.tauk[0]))
    

    results2.KtaxRev = 100*(results2.KtaxRev/results2.KtaxRev[0]-1)

    results2['Rev_lost'] = results2.KtaxRev_lost + results2.WtaxRev_lost
    results2.Rev_lost = (results2.Rev_lost/results2.Rev_lost[0])
    #results2.Rev_lost  = 100*results2.Rev_lost/results2.Y
    #results2.Rev_lost = results2.Rev_lost - results2.Rev_lost[0]
    
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
    results2.A_hid = (results2.A_hid/results2.A_hid[0])
    results2.P = 100*(results2.P/results2.P[0]-1.0)
    results2.r = 100*(results2.r-results2.r[0])
    #results2.lump_sum = 100*(results2.lump_sum/results2.ylbar)
    results2.welfare = 100*((results2.welfare/results2.welfare[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.welfare0 = 100*((results2.welfare0/results2.welfare0[0])**(1.0/(g*(1.0-4)))-1.0)
    results2.approval0 = 100*(results2.approval0)
    
    results2.p999 = 100*(1.0-results2.p999)
    results2.p999 = results2.p999 - results2.p999[0]
    results2.p90 = results2.p90 - results2.p90[0]
    results2.p999_total = 100*(1.0-results2.p999_total)
    results2.p999_total = results2.p999_total - results2.p999_total[0]
    results2.p90_total = results2.p90_total - results2.p90_total[0]


    results2['Rev_lost_f'] = results2.Rev_lost
    results2['A_hid_f'] = results2.A_hid
    results2['lump_sum_f'] = results2.lump_sum
    
    results2['Rev_lost_f'][2:] = nd.gaussian_filter1d(results2.Rev_lost[2:],0.75)
    results2['A_hid_f'][2:] = nd.gaussian_filter1d(results2.A_hid[2:],0.75)
    results2['lump_sum_f'][2:] = nd.gaussian_filter1d(results2.lump_sum[2:],0.75)
    
    results2 = results2[1:].reset_index(drop=True)

    return results2

evasion_type = ['evasion','no_evasion']
scenario = ['tauk','taua','warren','covid']
results = [[load_data(i,j) for i in range(len(evasion_type))] for j in range(len(scenario))]

####################################################
# make figures


colors = ['#377eb8','#e41a1c','#4daf4a','#984ea3','#ff7f00','#ffff33']
dashes = [(None,None),(6,1),(2,1)]

fig, AX = plt.subplots(4,4,figsize=(8.25,5.25),sharex=False,sharey=False)
#gs = gridspec.GridSpec(2,6)

for i in range(4):
    for j in range(4):
        AX[j][i].set_xlim(1,50)
        AX[j][i].set_xticks([1,10,20,30,40,50])
    AX[3][i].set_xlim(1,20)
    AX[3][i].set_xticks([1,5,10,15,20])
    AX[3][i].set_xlabel('Periods since policy change')



AX[0][0].set_title('Tax evasion (bench. = 1)',size=8)
AX[0][1].set_title('Transfer (\\% avg. labor inc.)',size=8)
AX[0][2].set_title('GDP (\\% chg.)',size=8)
AX[0][3].set_title('Welfare (\\% chg.)',size=8)

AX[0][0].set_ylabel('10p.p. increase in tauk',size=8)
AX[1][0].set_ylabel('2\% wealth tax',size=8)
AX[2][0].set_ylabel('Warren (perm)',size=8)
AX[3][0].set_ylabel('Warren (temp)',size=8)

# --------------------------
# figure 1a: evasion

labels = ['Evasion','No evasion']

for j in range(4):

    axes=AX[j][0]
    if(results[j][0] is not None):
        c1='A_hid'
        c2 = 'Rev_lost'
        if(j==3):
            c1 = 'A_hid_f'
            c2 = 'Rev_lost_f'
            
        ln1=axes.plot(results[j][0]['t'],results[j][0][c1],color=colors[0],dashes=dashes[0],label='Hidden wealth',linewidth=2,alpha=0.8)
        ln2=axes.plot(results[j][0]['t'],results[j][0][c2],color=colors[1],dashes=dashes[1],label='Lost revenues',linewidth=2,alpha=0.8)
        lns=ln1+ln2
        labs=[l.get_label() for l in lns]

    if(j==3):
        axes.legend(lns,labs,loc='lower right',prop={'size':6})

    axes=AX[j][1]
    for i in range(2):
        if(results[j][i] is not None):
            c='lump_sum'
            axes.plot(results[j][i]['t'],results[j][i][c],color=colors[i],dashes=dashes[i],label=labels[i],linewidth=2,alpha=0.8)

    if(j==0):
        axes.legend(loc='upper right',prop={'size':6})

    axes=AX[j][2]
    for i in range(2):
        if(results[j][i] is not None):
            c='Y'
            axes.plot(results[j][i]['t'],results[j][i][c],color=colors[i],dashes=dashes[i],label=labels[i],linewidth=2,alpha=0.8)

    if(j==0):
        axes.legend(loc='upper right',prop={'size':6})

    axes=AX[j][3]
    for i in range(2):
        if(results[j][i] is not None):
            c='welfare0'
            axes.plot(results[j][i]['t'],results[j][i][c],color=colors[i],dashes=dashes[i],label=labels[i],linewidth=2,alpha=0.8)

    if(j==0):
        axes.legend(loc='upper right',prop={'size':6})

plt.tight_layout(pad=0,h_pad=1,w_pad=0)
fig.subplots_adjust(hspace=0.2,wspace=0.25)
plt.savefig('output/fig/fig3_trans.pdf',bbox='tight',dpi=300)

##############################################################################
######################################################

print('Making latex tables')

T = [0,4,9,-1]
T2 = [0,1,4,9]

def fmt_res(file,tmp):
    if(np.isnan(tmp) or abs(tmp)<1.0e-9):
        file.write('&--')
    else:
        file.write('&%0.2f'%tmp)

vars = ['lump_sum','A_hid','Rev_lost','Y','welfare0','approval0']
labs = ['Transfer (\% avg. wage)','Concealed wealth (benchmark = 1)','Lost revenues (benchmark = 1)',
        'GDP (\% chg.)','Welfare (\% chg.)','Approval rate (\%)']
flags = [1,0,0,1,1,1]

panels = ['(a) Baseline model','(b) No-evasion counterfactual']

file=open('output/tex/table7_trans.tex','w')

file.write('\\footnotesize\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')

file.write('\\begin{tabular}{lrcccccccccccccccc}\n')
file.write('\\toprule\n')
file.write('& &\\multicolumn{4}{c}{\\makecell{10p.p. increase in\\\\capital income tax}} & \\multicolumn{4}{c}{\\makecell{2\\% wealth tax}} & \\multicolumn{4}{c}{\\makecell{Warren wealth tax\\\\(permanent)}} & \\multicolumn{4}{c}{\\makecell{Warren wealth tax\\ (temporary)}}\\\\\n')
file.write('\\cmidrule(rl){3-6}\\cmidrule(rl){7-10}\\cmidrule(rl){11-14}\\cmidrule(rl){15-18}\n')
file.write('Outcome & $t$ = ')

for i in range(3):
    file.write('& 1 & 5 & 10 & $\infty$')
file.write('& 1 & 2 & 5 & 10')
file.write('\\\\\n')

for i in range(len(panels)):

    file.write('\\midrule\n')
    file.write('\\textit{%s}&&&&&&&&&&&&&&&&&\\\\\n'%panels[i])

    for v,l,f in zip(vars,labs,flags):

        if(i==0 or f):
            
            file.write('%s & '%l)

            v2=v
            if(v=='Rev_lost' or v=='A_hid'):
                v2 = v+'_f'
        
            for t in T:
                tmp = 0
                if(results[0][i] is not None):
                    tmp = results[0][i][v].values[t]
                fmt_res(file,tmp)

            for t in T:
                tmp = 0
                if(results[1][i] is not None):
                    tmp = results[1][i][v].values[t]
                fmt_res(file,tmp)

            for t in T:
                tmp = 0
                if(results[2][i] is not None):
                    tmp = results[2][i][v].values[t]
                fmt_res(file,tmp)

            for t in T2:
                tmp = 0
                if(results[3][i] is not None):
                    tmp = results[3][i][v2].values[t]
                fmt_res(file,tmp)


            file.write('\\\\\n')

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
#file.write('\\end{center}\n')
file.write('\\normalsize\n')
#file.write('\\end{table}\n')
file.close()


plt.close('all')

