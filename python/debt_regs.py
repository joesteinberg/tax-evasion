import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
import weighted

df = pd.read_csv('../c/output/regdata.csv')

df.debt = df.capital - df.assets
df.loc[df.debt<0,'debt']=0
df['leverage'] = df.capital/(df.capital - df.debt)
df['output'] = df.z*df.capital
df['constrained'] = np.abs(df.kconst - df.capital)<1.0e-5

df=df[(df.debt>0) & (df.capital>0) & (np.isfinite(df.leverage))]

tmp1 = weighted.median(df.output/df.capital,df.weight)
tmp2 = np.average(df.output/df.capital,weights=df.weight)
df2 = df[(df.constrained==True) & (df.output/df.capital>tmp2)]

###############################

x1 = np.log(df2.output.values/df2.capital.values)

w = df2.weight.values
y = np.log(df2.leverage)

#X = sm.add_constant(x1)
X=x1

mod_wls = sm.WLS(y, X, weights=w)
res_wls = mod_wls.fit()
print(res_wls.summary())

###############################

x1 = df.profit.values/df.capital.values

w = df.weight.values
y = df.debt.values/df.capital.values

#X = np.column_stack((x1,x2))
X = sm.add_constant(x1)

mod_wls = sm.WLS(y, X, weights=w)
res_wls = mod_wls.fit()
print(res_wls.summary())
