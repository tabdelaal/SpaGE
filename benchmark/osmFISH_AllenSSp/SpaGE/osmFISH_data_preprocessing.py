import os
os.chdir('osmFISH_AllenSSp/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle
import loompy

## Read loom data
ds = loompy.connect('data/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')

FISH_Genes = ds.ra['Gene']
    
colAtr = ds.ca.keys()

df = pd.DataFrame()
for i in colAtr:
    df[i] = ds.ca[i]

osmFISH_meta = df.iloc[np.where(df.Valid == 1)[0], :]
osmFISH_data = ds[:,:]
osmFISH_data = osmFISH_data[:,np.where(df.Valid == 1)[0]]
osmFISH_data = pd.DataFrame(data= osmFISH_data, index= FISH_Genes)

del ds, colAtr, i, df, FISH_Genes

Cortex_Regions = ['Layer 2-3 lateral', 'Layer 2-3 medial', 'Layer 3-4', 
                  'Layer 4','Layer 5', 'Layer 6', 'Pia Layer 1']
Cortical = np.stack(i in Cortex_Regions for i in osmFISH_meta.Region)

osmFISH_meta = osmFISH_meta.iloc[Cortical,:]
osmFISH_data = osmFISH_data.iloc[:,Cortical]

cell_count = np.sum(osmFISH_data,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

osmFISH_data = osmFISH_data.apply(Log_Norm,axis=0)
osmFIH_data_scaled = pd.DataFrame(data=st.zscore(osmFISH_data.T),index = osmFISH_data.columns,columns=osmFISH_data.index)

datadict = dict()
datadict['osmFISH_data'] = osmFISH_data.T
datadict['osmFIH_data_scaled'] = osmFIH_data_scaled
datadict['osmFISH_meta'] = osmFISH_meta

with open('data/SpaGE_pkl/osmFISH_Cortex.pkl','wb') as f:
    pickle.dump(datadict, f)