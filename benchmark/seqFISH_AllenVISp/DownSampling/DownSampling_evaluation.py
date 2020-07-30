import os
os.chdir('seqFISH_AllenVISp/')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import scipy.stats as st
import pickle

### Original data
with open ('data/SpaGE_pkl/seqFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

seqFISH_data = datadict['seqFISH_data']
del datadict

test_set = pd.read_csv('Results/10SpaGE_New_genes.csv',header=0,index_col=0,sep=',').columns
seqFISH_data = seqFISH_data[test_set]

### SpaGE
#10
SpaGE_imputed_10 = pd.read_csv('Results/10SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_10 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_10[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_10[i])[0]
  
#30
SpaGE_imputed_30 = pd.read_csv('Results/30SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_30 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_30[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_30[i])[0]

#50    
SpaGE_imputed_50 = pd.read_csv('Results/50SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_50 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_50[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_50[i])[0]

#100    
SpaGE_imputed_100 = pd.read_csv('Results/100SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_100 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_100[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_100[i])[0]

#200    
SpaGE_imputed_200 = pd.read_csv('Results/200SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_200 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_200[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_200[i])[0]

#500    
SpaGE_imputed_500 = pd.read_csv('Results/500SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_500 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_500[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_500[i])[0]

#1000    
SpaGE_imputed_1000 = pd.read_csv('Results/1000SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_1000 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_1000[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_1000[i])[0]

#2000    
SpaGE_imputed_2000 = pd.read_csv('Results/2000SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_2000 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_2000[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_2000[i])[0]

#5000    
SpaGE_imputed_5000 = pd.read_csv('Results/5000SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_5000 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_5000[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_5000[i])[0]

#7000    
SpaGE_imputed_7000 = pd.read_csv('Results/7000SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_7000 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_7000[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_7000[i])[0]

#9701   
SpaGE_imputed_9701 = pd.read_csv('Results/9701SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_Corr_9701 = pd.Series(index = test_set)
for i in test_set:
    SpaGE_Corr_9701[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed_9701[i])[0]

### Comparison plots
plt.style.use('ggplot')
plt.figure(figsize=(9, 3))
plt.boxplot([SpaGE_Corr_10, SpaGE_Corr_30, SpaGE_Corr_50,
             SpaGE_Corr_100,SpaGE_Corr_200,SpaGE_Corr_500,
             SpaGE_Corr_1000,SpaGE_Corr_2000,SpaGE_Corr_5000,
             SpaGE_Corr_7000,SpaGE_Corr_9701])

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
y = SpaGE_Corr_1000
x = np.random.normal(7, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_2000
x = np.random.normal(8, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_5000
x = np.random.normal(9, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_7000
x = np.random.normal(10, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_Corr_9701
x = np.random.normal(11, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)


plt.xticks((1,2,3,4,5,6,7,8,9,10,11),
           ('10','30','50','100','200','500','1000','2000','5000','7000','9701'),size=10)
plt.yticks(size=8)
plt.xlabel('Number of shared genes',size=12)
plt.gca().set_ylim([-0.3,0.8])
plt.ylabel('Spearman Correlation',size=12)
plt.show()
