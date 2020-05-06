import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

### Timing
experiments = ['STARmap_AllenVISp','osmFISH_Ziesel','osmFISH_AllenVISp',
               'osmFISH_AllenSSp','MERFISH_Moffit']
# SpaGE
SpaGE_Avg_Time = pd.Series(index = experiments)
for i in experiments:
    SpaGE_Precise_Time = np.array(pd.read_csv(i + '/Results/SpaGE_PreciseTime.csv',
                                          header=0,sep=',')).reshape(-1)
    SpaGE_knn_Time = np.array(pd.read_csv(i +'/Results/SpaGE_knnTime.csv',
                                      header=0,sep=',')).reshape(-1)

    SpaGE_total_time = SpaGE_Precise_Time + SpaGE_knn_Time
    SpaGE_Avg_Time[i] = np.mean(SpaGE_total_time)
del SpaGE_Precise_Time,SpaGE_knn_Time,SpaGE_total_time

# gimVI
gimVI_Avg_Time = pd.Series(index = experiments)
for i in experiments:
    gimVI_Time = np.array(pd.read_csv(i+'/Results/gimVI_Time.csv',header=0,sep=',')).reshape(-1)
    gimVI_Avg_Time[i] = np.mean(gimVI_Time)
del gimVI_Time

# Seurat
Seurat_Avg_Time = pd.Series(index = experiments)
for i in experiments:
    Seurat_anchor_Time = np.array(pd.read_csv(i+'/Results/Seurat_anchor_time.csv',
                                              header=0,sep=',')).reshape(-1)
    Seurat_transfer_Time = np.array(pd.read_csv(i+'/Results/Seurat_transfer_time.csv'
                                                ,header=0,sep=',')).reshape(-1)
    
    Seurat_total_time = Seurat_anchor_Time + Seurat_transfer_Time
    Seurat_Avg_Time[i] = np.mean(Seurat_total_time)
del Seurat_anchor_Time,Seurat_transfer_Time,Seurat_total_time

# Liger
Liger_Avg_Time = pd.Series(index = experiments)
for i in experiments:
    Liger_NMF_Time = np.array(pd.read_csv(i+'/Results/Liger_NMF_time.csv',
                                          header=0,sep=',')).reshape(-1)
    Liger_knn_Time = np.array(pd.read_csv(i+'/Results/Liger_knn_time.csv',
                                          header=0,sep=',')).reshape(-1)
    
    Liger_total_time = Liger_NMF_Time + Liger_knn_Time
    Liger_Avg_Time[i] = np.mean(Liger_total_time)
del Liger_NMF_Time,Liger_knn_Time,Liger_total_time

plt.style.use('ggplot')
plt.figure(figsize=(9, 3))
plt.plot([1,2,3,4,5],SpaGE_Avg_Time,color='black',marker='s',linewidth=3)
plt.plot([1,2,3,4,5],Seurat_Avg_Time,color='blue',marker='s',linewidth=3)
plt.plot([1,2,3,4,5],Liger_Avg_Time,color='red',marker='s',linewidth=3)
plt.plot([1,2,3,4,5],gimVI_Avg_Time,color='purple',marker='s',linewidth=3)
#plt.yscale('log')
plt.xticks((1,2,3,4,5),('STARmap_AllenVISp\n(1549,14249)','osmFISH_Zeisel\n(3405,1691)',
            'osmFISH_AllenVISp\n(3405,14249)','osmFISH_AllenSSp\n(3405,5613)',
            'MERFSIH_Moffit\n(64373,31299)'),size=10)
plt.yticks(size=8)
plt.ylabel('Avergae computation time (seconds)',size=12)
colors = ['black','blue', 'red', 'purple']
lines = [Line2D([0], [0], color=c, linewidth=3, marker='s') for c in colors]
labels = ['SpaGE','Seurat', 'Liger','gimVI']
plt.legend(lines, labels)
plt.show()