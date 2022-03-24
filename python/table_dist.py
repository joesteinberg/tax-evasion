import glob
import pandas as pd
import numpy as np
import weighted

########################################################
print('Loading distributional results')

J = 61
NE=5
NZ=7
NI=2
NV=12
ND=2
NA=200
sigma=4.0

evasion_type = ['evasion','no_evasion']
gama=[0.43172371688765909, 0.42913860042472707]
inpath = '../c/output/'


V0 = [None,None]
Psi0 = [None,None]
V_k = [None,None]
Psi_k = [None,None]
V_a = [None,None]
Psi_a = [None,None]
V_ap = [None,None]
Psi_ap = [None,None]

for i in range(2):
        
    e = evasion_type[i]

    s0 = inpath + "ss0_%s" % e
    
    if(i==1):
        V20 = np.fromfile(s0 + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
        Psi20 = np.fromfile(s0 + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
    else:
        V20 = np.fromfile(s0 + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))
        Psi20 = np.fromfile(s0 + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))

    V0[i] = V20
    Psi0[i] = Psi20

        
    s = inpath + "ss1_tauk_0.344737"
    s = s+'_'+e

    if(i==1):
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
    else:
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))

    V_k[i]= V2
    Psi_k[i] = Psi2

    s = inpath + "ss1_taua_0.021053"
    s = s+'_'+e

    if(i==1):
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
    else:
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))

    V_a[i] = V2
    Psi_a[i] = Psi2

                
    s = inpath + "ss1_warren_%s" % e
        
    if(e=='no_evasion'):
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NA))
    else:
        V2 = np.fromfile(s + "_value_function.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))
        Psi2 = np.fromfile(s + "_dist.bin",dtype=float,count=-1).reshape((J,NE,NZ,NI,NV,ND,NA))

    V_ap[i] = V2
    Psi_ap[i] = Psi2


########################################################
print('Calculating welfare by wealth group')


cutoffs = [20,40,60,80,90,95,99,99.9,99.99]

indices0=[None,None]
indices_k=[None,None]
indices_a=[None,None]
indices_ap=[None,None]
      
for e in range(len(evasion_type)):
          
    total_mass = Psi0[e].sum()
    cdf=np.zeros(NA)
    if(e==1):
        cdf[0]=Psi0[e][:,:,:,:,0].sum()/total_mass
    else:
        cdf[0]=Psi0[e][:,:,:,:,:,:,0].sum()/total_mass
        
    for ia in range(1,NA):
        if e==1:
            cdf[ia]=cdf[ia-1]+Psi0[e][:,:,:,:,ia].sum()/total_mass
        else:
            cdf[ia]=cdf[ia-1]+Psi0[e][:,:,:,:,:,:,ia].sum()/total_mass

    indices2=[np.where(cdf<=c/100.0)[0].max()+1 for c in cutoffs]
    indices0[e] = indices2

    total_mass = Psi0[e].sum()
    cdf_k=np.zeros(NA)
    cdf_a=np.zeros(NA)
    cdf_ap=np.zeros(NA)

    if(e==1):
        cdf_k[0]=Psi_k[e][:,:,:,:,0].sum()/total_mass
        cdf_a[0]=Psi_a[e][:,:,:,:,0].sum()/total_mass
        cdf_ap[0]=Psi_ap[e][:,:,:,:,0].sum()/total_mass
    else:
        cdf_k[0]=Psi_k[e][:,:,:,:,:,:,0].sum()/total_mass
        cdf_a[0]=Psi_a[e][:,:,:,:,:,:,0].sum()/total_mass
        cdf_ap[0]=Psi_ap[e][:,:,:,:,:,:,0].sum()/total_mass
                
    for ia in range(1,NA):
        if e==1:
            cdf_k[ia]=cdf_k[ia-1]+Psi_k[e][:,:,:,:,ia].sum()/total_mass
            cdf_a[ia]=cdf_a[ia-1]+Psi_a[e][:,:,:,:,ia].sum()/total_mass
            cdf_ap[ia]=cdf_ap[ia-1]+Psi_ap[e][:,:,:,:,ia].sum()/total_mass
        else:
            cdf_k[ia]=cdf_k[ia-1]+Psi_k[e][:,:,:,:,:,:,ia].sum()/total_mass
            cdf_a[ia]=cdf_a[ia-1]+Psi_a[e][:,:,:,:,:,:,ia].sum()/total_mass
            cdf_ap[ia]=cdf_ap[ia-1]+Psi_ap[e][:,:,:,:,:,:,ia].sum()/total_mass

        indices2=[ np.append(np.where(cdf_k<=c/100.0)[0],0).max()+1 for c in cutoffs]
        indices_k[e] = indices2

        indices2=[ np.append(np.where(cdf_a<=c/100.0)[0],0).max()+1 for c in cutoffs]
        indices_a[e] = indices2

        indices2=[ np.append(np.where(cdf_ap<=c/100.0)[0],0).max()+1 for c in cutoffs]
        indices_ap[e] = indices2

def calc_CE(V,Psi,V0,Psi0,g):
    tmp = (V/V0)**(1.0/(g*(1.0-sigma)))-1.0
    tmp[Psi0<1e-11] = 0.0
    tmp2 = weighted.median(tmp.flatten(),Psi0.flatten())
    #tmp2 = np.multiply(tmp,Psi0).sum()/Psi0.sum()
    return tmp2

            
W_by_a_k=[None,None]
W_by_a_a=[None,None]
W_by_a_ap=[None,None]

for e in range(2):

    W_by_a3_k=np.zeros((len(cutoffs)+1))
    
    for i in range(len(cutoffs)+1):
        lb0=0
        ub0=0
        lb=0
        ub=0
            
        if i==0:
            lb0=0
            ub0=indices0[e][i]+1
            lb=0
            ub=indices_k[e][i]+1
                
        elif i<len(cutoffs):
            lb0=indices0[e][i-1]+1
            ub0=indices0[e][i]+1
            lb=indices_k[e][i-1]+1
            ub=indices_k[e][i]+1
                    
        else:
            lb0=indices0[e][i-1]+1
            ub0=NA
            lb=indices_k[e][i-1]+1
            ub=NA
                
        if e==1:
            W_by_a3_k[i] = 100*calc_CE(V_k[e][:,:,:,:,lb:ub],
                                       Psi_k[e][:,:,:,:,lb:ub],
                                       V0[e][:,:,:,:,lb:ub],
                                       Psi0[e][:,:,:,:,lb:ub],
                                       gama[e])
                

        else:
            W_by_a3_k[i] = 100*calc_CE(V_k[e][:,:,:,:,:,:,lb:ub],
                                       Psi_k[e][:,:,:,:,:,:,lb:ub],
                                       V0[e][:,:,:,:,:,:,lb:ub],
                                       Psi0[e][:,:,:,:,:,:,lb:ub],
                                       gama[e])
            
    W_by_a3_a=np.zeros((len(cutoffs)+1))
        
    for i in range(len(cutoffs)+1):
        lb0=0
        ub0=0
        lb=0
        ub=0
        
        if i==0:
            lb0=0
            ub0=indices0[e][i]+1
            lb=0
            ub=indices_a[e][i]+1
                
        elif i<len(cutoffs):
            lb0=indices0[e][i-1]+1
            ub0=indices0[e][i]+1
            lb=indices_a[e][i-1]+1
            ub=indices_a[e][i]+1
                    
        else:
            lb0=indices0[e][i-1]+1
            ub0=NA
            lb=indices_a[e][i-1]+1
            ub=NA
                
        if e==1:
            W_by_a3_a[i] = 100*calc_CE(V_a[e][:,:,:,:,lb:ub],
                                       Psi_a[e][:,:,:,:,lb:ub],
                                       V0[e][:,:,:,:,lb:ub],
                                       Psi0[e][:,:,:,:,lb:ub],
                                       gama[e])
                

        else:
            W_by_a3_a[i] = 100*calc_CE(V_a[e][:,:,:,:,:,:,lb:ub],
                                       Psi_a[e][:,:,:,:,:,:,lb:ub],
                                       V0[e][:,:,:,:,:,:,lb:ub],
                                       Psi0[e][:,:,:,:,:,:,lb:ub],
                                       gama[e])


    W_by_a3_ap=np.zeros((len(cutoffs)+1))
        
    for i in range(len(cutoffs)+1):
        lb0=0
        ub0=0
        lb=0
        ub=0
            
        if i==0:
            lb0=0
            ub0=indices0[e][i]+1
            lb=0
            ub=indices_ap[e][i]+1
                
        elif i<len(cutoffs):
            lb0=indices0[e][i-1]+1
            ub0=indices0[e][i]+1
            lb=indices_ap[e][i-1]+1
            ub=indices_ap[e][i]+1
                    
        else:
            lb0=indices0[e][i-1]+1
            ub0=NA
            lb=indices_ap[e][i-1]+1
            ub=NA
                
        if e==1:
            W_by_a3_ap[i] = 100*calc_CE(V_ap[e][:,:,:,:,lb:ub],
                                        Psi_ap[e][:,:,:,:,lb:ub],
                                        V0[e][:,:,:,:,lb:ub],
                                        Psi0[e][:,:,:,:,lb:ub],
                                        gama[e])
                

        else:
            W_by_a3_ap[i] = 100*calc_CE(V_ap[e][:,:,:,:,:,:,lb:ub],
                                        Psi_ap[e][:,:,:,:,:,:,lb:ub],
                                        V0[e][:,:,:,:,:,:,lb:ub],
                                        Psi0[e][:,:,:,:,:,:,lb:ub],
                                        gama[e])



    W_by_a_k[e] = W_by_a3_k
    W_by_a_a[e] = W_by_a3_a
    W_by_a_ap[e] = W_by_a3_ap



def fmt_w_res(file,x):
    s=''
    if(x<-1e-4):
        s= '&\\textcolor{Red}{%+0.3f}'%x
    elif(x>1e-4):
        s= '&\\textcolor{Green}{%+0.3f}'%x
    elif(abs(x)<1.0e-10):
        s='&--'
    else:
        s= '&%0.3f'%x

    file.write(s)
    
panels = ['(a) Baseline model','(b) No-evasion counterfactual']


file=open('output/tex/table3_dist_results.tex','w')

file.write('\\footnotesize\n')
file.write('\\renewcommand{\\arraystretch}{1.2}\n')
file.write('\\renewcommand{\\cellalign}{bc}\n')
file.write('\\begin{tabular}{lccc}\n')
file.write('\\toprule\n')
#file.write('&\\multicolumn{2}{c}{2p.p. increase in}\\\\\n')
#file.write('\\cmidrule(rl){2-3}\n')
file.write('Percentile & \\makecell{10p.p. increase in\\\\capital income tax} & \\makecell{2\\% flat wealth tax} & Warren wealth tax\\\\\n')

for i in range(len(panels)):

    file.write('\\midrule\n')
    file.write('\\multicolumn{4}{l}{\\textit{%s}}\\\\\n'%panels[i])

    for j in range(len(cutoffs)+1):
        lab=''
        if j==0:
            lab='0--%d' % cutoffs[j]
        elif j<len(cutoffs)-2:
            lab='%d--%d' % (cutoffs[j-1],cutoffs[j])
        elif j==len(cutoffs)-2:
            lab='%d--%0.1f' % (cutoffs[j-1],cutoffs[j])
        elif j==len(cutoffs)-1:
            lab='%0.1f--%0.2f' % (cutoffs[j-1],cutoffs[j])
        else:
            lab='%0.2f--100' % cutoffs[j-1]

        file.write(lab)

                        
        tmp = W_by_a_k[i][j]
        fmt_w_res(file,tmp)
        
        tmp = W_by_a_a[i][j]
        fmt_w_res(file,tmp)

        tmp = W_by_a_ap[i][j]
        fmt_w_res(file,tmp)

        if(i==0 and j==len(cutoffs)):
            file.write('\\\\\n[1.0ex]')
        else:
            file.write('\\\\\n')
                

file.write('\\bottomrule\n')
file.write('\\end{tabular}\n')
file.write('\\normalsize\n')
file.close()
