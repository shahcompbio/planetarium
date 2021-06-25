from jinja2 import Template
import subprocess
import pandas as pd
import sys
import os
import json

import scanpy as sc
from scipy.sparse import csr_matrix, find
import json


def get_data(filepath):
    adata = sc.read(filepath)

    genes = get_gene_data(adata)
    metadata = get_metadata(adata)


def get_gene_data(adata):
    genes = adata.var
    genes = genes.reset_index(drop=True)
    genes = genes.reset_index().rename(columns={'index': 'gene_idx'})

    ens = pd.read_csv('ensemble_ids.tsv', sep='\t')
    ens = ens.rename(columns={'Unnamed: 0': 'gene'})

    genes = genes.merge(ens, how='left')

    cells = adata.obs
    cells = cells.reset_index().rename(columns={'index': 'cell_id'})
    cells = cells.reset_index().rename(columns={'index': 'cell_idx'})

    normalized_matrix = csr_matrix(adata.X)
    indices = normalized_matrix < 0.0
    normalized_matrix[indices] = 0.0
    nonzero = find(normalized_matrix)

    df = pd.DataFrame(list(zip(*nonzero)), columns=['cell_idx', 'gene_idx', 'expression'])

    df = df.merge(genes, how='left')
    df = df.merge(cells, how='left')

    df = df[df['gene'].notna()]
    df = df[['cell_id', 'gene', 'expression']]

    return df

def get_metadata(adata):
    umap = pd.DataFrame(adata.obsm['X_umap'])
    umap.columns = ['UMAP_1', 'UMAP_2']
    df = adata.obs[['clone', 'timepoint', 'therapy']]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'index': 'cell_id'})
    df = df.replace(to_replace="nan", value="None")
    return df


data_dir = sys.argv[1]
metadata = pd.read_csv(os.path.join(data_dir, "vdjTimeSeries.tsv"), sep="\t")
metadata = metadata.to_dict('records')
with open(os.path.join(data_dir, "genes.json")) as f:
  genes = json.load(f)
print(len(genes))
#genes = pd.read_csv(os.path.join(data_dir, "genes.json"), sep="")
#genes = genes.to_dict("records")
data = {
    "metadata": metadata,
    "genes": genes,
}
data = json.dumps(data, indent=4)
app_dir = os.path.dirname(os.path.abspath(__file__))
app_dir = os.path.abspath(os.path.join(app_dir, "../..", "build"))
index_template=os.path.join(app_dir, "index.html")
template = Template(open(index_template,"r").read())
html = template.render(data=data)
output_html=os.path.join(app_dir,"Timeseries.html")
output = open(output_html,"w")
output.write(html)
output.close()
# datalake_build=os.path.join(sys.argv[1], "build")
# shutil.copytree(build_folder, datalake_build)