import os
os.chdir('MERFISH_Moffit/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

MERFISH = pd.read_csv('data/MERFISH/Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv')
#Select 1st replicate, Naive female
MERFISH_1 = MERFISH.loc[MERFISH['Animal_ID']==1,:]

#remove Blank-1 to 5 and Fos --> 155 genes
MERFISH_1 = MERFISH_1.loc[MERFISH_1['Cell_class']!='Ambiguous',:]
MERFISH_meta = MERFISH_1.iloc[:,0:9]
MERFISH_data = MERFISH_1.iloc[:,9:171]
MERFISH_data = MERFISH_data.drop(columns = ['Blank_1','Blank_2','Blank_3','Blank_4','Blank_5','Fos'])
del MERFISH, MERFISH_1

MERFISH_data = MERFISH_data.T

cell_count = np.sum(MERFISH_data,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

MERFISH_data = MERFISH_data.apply(Log_Norm,axis=0)
MERFISH_data_scaled = pd.DataFrame(data=st.zscore(MERFISH_data.T),index = MERFISH_data.columns,columns=MERFISH_data.index)

datadict = dict()
datadict['MERFISH_data'] = MERFISH_data.T
datadict['MERFISH_data_scaled'] = MERFISH_data_scaled
datadict['MERFISH_meta'] = MERFISH_meta

with open('data/SpaGE_pkl/MERFISH.pkl','wb') as f:
    pickle.dump(datadict, f)