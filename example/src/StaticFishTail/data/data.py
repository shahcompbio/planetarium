import scanpy as sc
import pandas as pd
import json
from scipy.sparse import csr_matrix, find

protein = sc.read('protein_h5ad.h5ad')

rna = sc.read('rna_h5ad.h5ad')

rna_umap = pd.DataFrame(rna.obsm['X_umap'])
rna_umap.columns = ['rna_UMAP_1', 'rna_UMAP_2']
rna_cells = rna.obs.reset_index().rename(columns={'index': 'cell_id'})
rna_umap = rna_umap.merge(rna_cells[['cell_id']], left_index=True, right_index=True)

protein_umap = pd.DataFrame(protein.obsm['X_umap'])
protein_umap.columns = ['protein_UMAP_1', 'protein_UMAP_2']
protein_cells = protein.obs.reset_index().rename(columns={'index': 'cell_id'})
protein_umap = protein_umap.merge(protein_cells[['cell_id']], left_index=True, right_index=True)

## join two tables together
df = protein_umap.merge(rna_umap)


## now we prep the matrix for DISTRIBUTION

rna_matrix = csr_matrix(rna.X)
nonzero = find(rna_matrix)
rna_matrix = pd.DataFrame(list(zip(*nonzero)), columns=['cell_idx', 'gene_idx', 'expression'])
rna_matrix = rna_matrix.merge(rna_cells[['cell_id']], left_on='cell_idx', right_index=True)

genes = rna.var.reset_index().rename(columns={'index': 'gene_id'})
rna_matrix = rna_matrix.merge(genes[['gene_id']], left_on='gene_idx', right_index=True)


protein_matrix = csr_matrix(protein.X)
nonzero = find(protein_matrix)
protein_matrix = pd.DataFrame(list(zip(*nonzero)), columns=['cell_idx', 'protein_idx', 'expression'])
protein_matrix = protein_matrix.merge(protein_cells[['cell_id']], left_on='cell_idx', right_index=True)

proteins = protein.var.reset_index().rename(columns={'index': 'protein_id'})
protein_matrix = protein_matrix.merge(proteins[['protein_id']], left_on='protein_idx', right_index=True)


## Now format into something usable

cells = df.to_dict(orient='records')

records = []
for cell in cells:
    rna_records = rna_matrix[rna_matrix['cell_id'] == cell['cell_id']]
    rna_obj = dict(zip(rna_records['gene_id'], rna_records['expression']))

    protein_records = protein_matrix[protein_matrix['cell_id'] == cell['cell_id']]
    protein_obj = dict(zip(protein_records['protein_id'], protein_records['expression']))

    record = {
        **cell,
        **rna_obj,
        **protein_obj
    }

    records.append(record)

with open('cells.json', 'w', encoding='utf-8') as file:
    json.dump(records, file, ensure_ascii=False, indent=4)

with open('genes.json', 'w', encoding='utf-8') as file:
    gene_list = list(genes['gene_id'])
    gene_list.sort()
    json.dump(gene_list, file, ensure_ascii=False, indent=4)


with open('proteins.json', 'w', encoding='utf-8') as file:
    protein_list = list(proteins['protein_id'])
    protein_list.sort()
    json.dump(protein_list, file, ensure_ascii=False, indent=4)
