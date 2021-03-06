import os
os.chdir('STARmap_AllenVISp/')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import scipy.stats as st
import pickle

with open ('data/SpaGE_pkl/Starmap.pkl', 'rb') as f:
    datadict = pickle.load(f)

Starmap_data = datadict['Starmap_data']
del datadict

test_set = pd.read_csv('Results/10SpaGE_New_genes.csv',header=0,index_col=0,sep=',').columns
Starmap_data = Starmap_data[test_set]

### SpaGE
#10
SpaGE_imputed_10 = pd.read_csv('Results/10SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_10 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_10[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_10[i])[0]
  
#30
SpaGE_imputed_30 = pd.read_csv('Results/30SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_30 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_30[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_30[i])[0]

#50    
SpaGE_imputed_50 = pd.read_csv('Results/50SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_50 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_50[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_50[i])[0]

#100    
SpaGE_imputed_100 = pd.read_csv('Results/100SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_100 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_100[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_100[i])[0]

#200    
SpaGE_imputed_200 = pd.read_csv('Results/200SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_200 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_200[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_200[i])[0]
    
#500    
SpaGE_imputed_500 = pd.read_csv('Results/500SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_500 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_500[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_500[i])[0]

#944    
SpaGE_imputed_944 = pd.read_csv('Results/944SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_944 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_944[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed_944[i])[0]

### Comparison plots
plt.style.use('ggplot')
plt.figure(figsize=(9, 3))
plt.boxplot([SpaGE_Corr_10, SpaGE_Corr_30, SpaGE_Corr_50,
             SpaGE_Corr_100,SpaGE_Corr_200,SpaGE_Corr_944])

y = SpaGE_Corr_10
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_30
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_50
x = np.random.normal(3, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_100
x = np.random.normal(4, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_200
x = np.random.normal(5, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_500
x = np.random.normal(6, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_944
x = np.random.normal(7, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)


plt.xticks((1,2,3,4,5,6,7),('10','30', '50','100','200','500','944'),size=10)
plt.yticks(size=8)
plt.xlabel('Number of shared genes',size=12)
plt.gca().set_ylim([-0.3,0.8])
plt.ylabel('Spearman Correlation',size=12)
plt.show()
