import os
os.chdir('osmFISH_Ziesel/')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import scipy.stats as st
from matplotlib.lines import Line2D
import pickle

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data = datadict['osmFISH_data']
del datadict

Gene_Order = osmFISH_data.columns

### SpaGE
SpaGE_imputed_Ziesel = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed_Ziesel = SpaGE_imputed_Ziesel.loc[:,Gene_Order]

SpaGE_Corr_Ziesel = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr_Ziesel[i] = st.spearmanr(osmFISH_data[i],SpaGE_imputed_Ziesel[i])[0]

### gimVI
gimVI_imputed_Ziesel = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')

gimVI_imputed_Ziesel = gimVI_imputed_Ziesel.loc[:,[x.upper() for x in np.array(Gene_Order,dtype='str')]]

gimVI_Corr_Ziesel = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr_Ziesel[i] = st.spearmanr(osmFISH_data[i],gimVI_imputed_Ziesel[str(i).upper()])[0]
gimVI_Corr_Ziesel[np.isnan(gimVI_Corr_Ziesel)] = 0

### Seurat
Seurat_imputed_Ziesel = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Seurat_imputed_Ziesel = Seurat_imputed_Ziesel.loc[:,Gene_Order]

Seurat_Corr_Ziesel = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr_Ziesel[i] = st.spearmanr(osmFISH_data[i],Seurat_imputed_Ziesel[i])[0]

### Liger
Liger_imputed_Ziesel = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Liger_imputed_Ziesel = Liger_imputed_Ziesel.loc[:,Gene_Order]

Liger_Corr_Ziesel = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr_Ziesel[i] = st.spearmanr(osmFISH_data[i],Liger_imputed_Ziesel[i])[0]
Liger_Corr_Ziesel[np.isnan(Liger_Corr_Ziesel)] = 0


os.chdir('osmFISH_AllenVISp/')

### SpaGE
SpaGE_imputed_AllenVISp = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed_AllenVISp = SpaGE_imputed_AllenVISp.loc[:,Gene_Order]

SpaGE_Corr_AllenVISp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr_AllenVISp[i] = st.spearmanr(osmFISH_data[i],SpaGE_imputed_AllenVISp[i])[0]

### gimVI
gimVI_imputed_AllenVISp = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')

gimVI_imputed_AllenVISp = gimVI_imputed_AllenVISp.loc[:,[x.upper() for x in np.array(Gene_Order,dtype='str')]]

gimVI_Corr_AllenVISp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr_AllenVISp[i] = st.spearmanr(osmFISH_data[i],gimVI_imputed_AllenVISp[str(i).upper()])[0]
gimVI_Corr_AllenVISp[np.isnan(gimVI_Corr_AllenVISp)] = 0

### Seurat
Seurat_imputed_AllenVISp = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Seurat_imputed_AllenVISp = Seurat_imputed_AllenVISp.loc[:,Gene_Order]

Seurat_Corr_AllenVISp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr_AllenVISp[i] = st.spearmanr(osmFISH_data[i],Seurat_imputed_AllenVISp[i])[0]

### Liger
Liger_imputed_AllenVISp = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Liger_imputed_AllenVISp = Liger_imputed_AllenVISp.loc[:,Gene_Order]

Liger_Corr_AllenVISp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr_AllenVISp[i] = st.spearmanr(osmFISH_data[i],Liger_imputed_AllenVISp[i])[0]
Liger_Corr_AllenVISp[np.isnan(Liger_Corr_AllenVISp)] = 0


os.chdir('osmFISH_AllenSSp/')

### SpaGE
SpaGE_imputed_AllenSSp = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed_AllenSSp = SpaGE_imputed_AllenSSp.loc[:,Gene_Order]

SpaGE_Corr_AllenSSp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr_AllenSSp[i] = st.spearmanr(osmFISH_data[i],SpaGE_imputed_AllenSSp[i])[0]

### gimVI
gimVI_imputed_AllenSSp = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')

gimVI_imputed_AllenSSp = gimVI_imputed_AllenSSp.loc[:,[x.upper() for x in np.array(Gene_Order,dtype='str')]]

gimVI_Corr_AllenSSp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr_AllenSSp[i] = st.spearmanr(osmFISH_data[i],gimVI_imputed_AllenSSp[str(i).upper()])[0]
gimVI_Corr_AllenSSp[np.isnan(gimVI_Corr_AllenSSp)] = 0

### Seurat
Seurat_imputed_AllenSSp = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Seurat_imputed_AllenSSp = Seurat_imputed_AllenSSp.loc[:,Gene_Order]

Seurat_Corr_AllenSSp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr_AllenSSp[i] = st.spearmanr(osmFISH_data[i],Seurat_imputed_AllenSSp[i])[0]

### Liger
Liger_imputed_AllenSSp = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Liger_imputed_AllenSSp = Liger_imputed_AllenSSp.loc[:,Gene_Order]

Liger_Corr_AllenSSp = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr_AllenSSp[i] = st.spearmanr(osmFISH_data[i],Liger_imputed_AllenSSp[i])[0]
Liger_Corr_AllenSSp[np.isnan(Liger_Corr_AllenSSp)] = 0

### Comparison plots
plt.style.use('ggplot')
fig, ax = plt.subplots(figsize=(5, 3))
ax.boxplot([SpaGE_Corr_Ziesel, Seurat_Corr_Ziesel, Liger_Corr_Ziesel,gimVI_Corr_Ziesel],
            positions = [1,5,9,13],boxprops=dict(color='red'),
            capprops=dict(color='red'),whiskerprops=dict(color='red'),flierprops=dict(color='red', markeredgecolor='red'),
            medianprops=dict(color='red'))
ax.boxplot([SpaGE_Corr_AllenVISp, Seurat_Corr_AllenVISp, Liger_Corr_AllenVISp,gimVI_Corr_AllenVISp],
            positions = [2,6,10,14],boxprops=dict(color='blue'),
            capprops=dict(color='blue'),whiskerprops=dict(color='blue'),flierprops=dict(color='blue', markeredgecolor='blue'),
            medianprops=dict(color='blue'))
ax.boxplot([SpaGE_Corr_AllenSSp, Seurat_Corr_AllenSSp, Liger_Corr_AllenSSp,gimVI_Corr_AllenSSp],
            positions = [3,7,11,15],boxprops=dict(color='purple'),
            capprops=dict(color='purple'),whiskerprops=dict(color='purple'),flierprops=dict(color='purple', markeredgecolor='purple'),
            medianprops=dict(color='purple'))
y = SpaGE_Corr_Ziesel
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = SpaGE_Corr_AllenVISp
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = SpaGE_Corr_AllenSSp
x = np.random.normal(3, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Seurat_Corr_Ziesel
x = np.random.normal(5, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Seurat_Corr_AllenVISp
x = np.random.normal(6, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Seurat_Corr_AllenSSp
x = np.random.normal(7, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Liger_Corr_Ziesel
x = np.random.normal(9, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Liger_Corr_AllenVISp
x = np.random.normal(10, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Liger_Corr_AllenSSp
x = np.random.normal(11, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = gimVI_Corr_Ziesel
x = np.random.normal(13, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = gimVI_Corr_AllenVISp
x = np.random.normal(14, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = gimVI_Corr_AllenSSp
x = np.random.normal(15, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)

plt.xticks((2,6,10,14),('SpaGE','Seurat', 'Liger','gimVI'),size=12)
plt.yticks(size=8)
plt.gca().set_ylim([-0.4,0.8])
plt.ylabel('Spearman Correlation',size=12)
colors = ['red', 'blue', 'purple']
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
labels = ['Ziesel', 'AllenVISp','AllenSSp']
plt.legend(lines, labels)
plt.show()
