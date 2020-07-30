import os
os.chdir('/osmFISH_AllenVISp/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

seqFISH_data = pd.read_csv('data/seqFISH/sourcedata/cortex_svz_counts.csv',header=0)
seqFISH_meta = pd.read_csv('data/seqFISH/sourcedata/cortex_svz_cellcentroids.csv',header=0)

seqFISH_data = seqFISH_data.iloc[np.where(seqFISH_meta['Region'] == 'Cortex')[0],:]
seqFISH_meta = seqFISH_meta.iloc[np.where(seqFISH_meta['Region'] == 'Cortex')[0],:]

seqFISH_data = seqFISH_data.T

cell_count = np.sum(seqFISH_data,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

seqFISH_data = seqFISH_data.apply(Log_Norm,axis=0)
seqFISH_data_scaled = pd.DataFrame(data=st.zscore(seqFISH_data.T),index = seqFISH_data.columns,columns=seqFISH_data.index)

datadict = dict()
datadict['seqFISH_data'] = seqFISH_data.T
datadict['seqFISH_data_scaled'] = seqFISH_data_scaled
datadict['seqFISH_meta'] = seqFISH_meta

with open('data/SpaGE_pkl/seqFISH_Cortex.pkl','wb') as f:
    pickle.dump(datadict, f)
