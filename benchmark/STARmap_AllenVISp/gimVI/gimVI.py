import os
os.chdir('STARmap_AllenVISp/')

from scvi.dataset import CsvDataset
from scvi.models import JVAE, Classifier
from scvi.inference import JVAETrainer
import numpy as np
import pandas as pd
import copy
import torch
import time as tm

### STARmap data
Starmap_data = CsvDataset('data/gimVI_data/STARmap_data_scvi.csv', save_path = "", sep = ",", gene_by_cell = False)

### AllenVISp
RNA_data = CsvDataset('data/gimVI_data/Allen_data_scvi.csv', save_path = "", sep = ",", gene_by_cell = False)

### Leave-one-out validation
Gene_set = np.intersect1d(Starmap_data.gene_names,RNA_data.gene_names)
Starmap_data.gene_names = Gene_set
Starmap_data.X = Starmap_data.X[:,np.reshape(np.vstack(np.argwhere(i==Starmap_data.gene_names) for i in Gene_set),-1)]
Common_data = copy.deepcopy(RNA_data)
Common_data.gene_names = Gene_set
Common_data.X = Common_data.X[:,np.reshape(np.vstack(np.argwhere(i==RNA_data.gene_names) for i in Gene_set),-1)]
Imp_Genes = pd.DataFrame(columns=Gene_set)
gimVI_time = []

for i in Gene_set:
    print(i)
    # Create copy of the fish dataset with hidden genes
    data_spatial_partial = copy.deepcopy(Starmap_data)
    data_spatial_partial.filter_genes_by_attribute(np.setdiff1d(Starmap_data.gene_names,i))
    data_spatial_partial.batch_indices += Common_data.n_batches
    
    datasets = [Common_data, data_spatial_partial]
    generative_distributions = ["zinb", "nb"]
    gene_mappings = [slice(None), np.reshape(np.vstack(np.argwhere(i==Common_data.gene_names) for i in data_spatial_partial.gene_names),-1)]
    n_inputs = [d.nb_genes for d in datasets]
    total_genes = Common_data.nb_genes
    n_batches = sum([d.n_batches for d in datasets])
    
    model_library_size = [True, False]
    
    n_latent = 8
    kappa = 5
    
    start = tm.time()
    torch.manual_seed(0)
    
    model = JVAE(
        n_inputs,
        total_genes,
        gene_mappings,
        generative_distributions,
        model_library_size,
        n_layers_decoder_individual=0,
        n_layers_decoder_shared=0,
        n_layers_encoder_individual=1,
        n_layers_encoder_shared=1,
        dim_hidden_encoder=64,
        dim_hidden_decoder_shared=64,
        dropout_rate_encoder=0.2,
        dropout_rate_decoder=0.2,
        n_batch=n_batches,
        n_latent=n_latent,
    )
    discriminator = Classifier(n_latent, 32, 2, 3, logits=True)
    
    trainer = JVAETrainer(model, discriminator, datasets, 0.95, frequency=1, kappa=kappa)
    trainer.train(n_epochs=200)
    _,Imputed = trainer.get_imputed_values(normalized=True)
    Imputed = np.reshape(Imputed[:,np.argwhere(i==Common_data.gene_names)[0]],-1)
    Imp_Genes[i] = Imputed
    gimVI_time.append(tm.time()-start)
    
Imp_Genes = Imp_Genes.fillna(0)    
Imp_Genes.to_csv('Results/gimVI_LeaveOneOut.csv')
gimVI_time = pd.DataFrame(gimVI_time)
gimVI_time.to_csv('Results/gimVI_Time.csv', index = False)

### New genes
Imp_New_Genes = pd.DataFrame(columns=["TESC","PVRL3","SOX10","GRM2","TCRB"])

# Create copy of the fish dataset with hidden genes
data_spatial_partial = copy.deepcopy(Starmap_data)
data_spatial_partial.filter_genes_by_attribute(Starmap_data.gene_names)
data_spatial_partial.batch_indices += RNA_data.n_batches

datasets = [RNA_data, data_spatial_partial]
generative_distributions = ["zinb", "nb"]
gene_mappings = [slice(None), np.reshape(np.vstack(np.argwhere(i==RNA_data.gene_names) for i in data_spatial_partial.gene_names),-1)]
n_inputs = [d.nb_genes for d in datasets]
total_genes = RNA_data.nb_genes
n_batches = sum([d.n_batches for d in datasets])

model_library_size = [True, False]

n_latent = 8
kappa = 5

torch.manual_seed(0)

model = JVAE(
    n_inputs,
    total_genes,
    gene_mappings,
    generative_distributions,
    model_library_size,
    n_layers_decoder_individual=0,
    n_layers_decoder_shared=0,
    n_layers_encoder_individual=1,
    n_layers_encoder_shared=1,
    dim_hidden_encoder=64,
    dim_hidden_decoder_shared=64,
    dropout_rate_encoder=0.2,
    dropout_rate_decoder=0.2,
    n_batch=n_batches,
    n_latent=n_latent,
)
discriminator = Classifier(n_latent, 32, 2, 3, logits=True)

trainer = JVAETrainer(model, discriminator, datasets, 0.95, frequency=1, kappa=kappa)
trainer.train(n_epochs=200)
    
for i in ["TESC","PVRL3","SOX10","GRM2","TCRB"]:
    _,Imputed = trainer.get_imputed_values(normalized=True)
    Imputed = np.reshape(Imputed[:,np.argwhere(i==RNA_data.gene_names)[0]],-1)
    Imp_New_Genes[i] = Imputed
    
Imp_New_Genes.to_csv('Results/gimVI_New_genes.csv')
