from jinja2 import Template
from elasticsearch import Elasticsearch
from elasticsearch import helpers
import pandas as pd
import numpy as np
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
    umap.columns = ['umap_1', 'umap_2']
    df = adata.obs[['clone', 'timepoint', 'therapy']]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'index': 'cell_id'})
    df = df.replace(to_replace="nan", value="None")
    return df

MAPPING = {
    'mappings': {
        "dynamic_templates": [
            {
                "string_values": {
                    "match": "*",
                    "match_mapping_type": "string",
                    "mapping": {
                        "type": "keyword"
                    }
                }
            }
        ]
    }
}

def load_df(df, index_name, host='localhost', port=9200, batch_size=int(1e5)):
    """Batch load dataframe"""
    total_records = df.shape[0]
    num_records = 0

    es = Elasticsearch(hosts=[{'host': host,'port': port }])

    for batch_start_idx in range(0, df.shape[0], batch_size):
        batch_end_idx = min(batch_start_idx + batch_size, df.shape[0])
        batch_data = df.loc[df.index[batch_start_idx:batch_end_idx]]

        clean_fields(batch_data)

        records = []
        for record in batch_data.to_dict(orient='records'):
            clean_nans(record)
            records.append(record)

        load_records(es, records, index_name)
        num_records += batch_data.shape[0]
        print(
            f"Loading {len(records)} records. Total: {num_records} / {total_records} ({(num_records * 100 / total_records): .1f}%)")
    if total_records != num_records:
        raise ValueError(
            'mismatch in {num_records} records loaded to {total_records} total records')

def load_records(es, records, index):
    """Load batch of records"""
    if not es.indices.exists(index):
        es.indices.create(index=index, body=MAPPING)

    for success, info in helpers.parallel_bulk(es, records, index=index):
        if not success:
            print('Doc failed in parallel loading', fg="red")
            print(info)

def clean_fields(df):
    """Remove invalid characters from column names"""
    invalid_chars = ['.']
    invalid_cols = [col for col in df.columns if any(
        [char in col for char in invalid_chars])]
    for col in invalid_cols:
        found_chars = [char for char in invalid_chars if char in col]

        for char in found_chars:
            new_col = col.replace(char, '_')
            df.rename(columns={col: new_col}, inplace=True)


def clean_nans(record):
    """Delete any fields that contain nan values"""
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            del record[field]


## inject data into built template

# data_dir = sys.argv[1]
# metadata = pd.read_csv(os.path.join(data_dir, "vdjTimeSeries.tsv"), sep="\t")
# metadata = metadata.to_dict('records')
# with open(os.path.join(data_dir, "genes.json")) as f:
#   genes = json.load(f)
# print(len(genes))
# #genes = pd.read_csv(os.path.join(data_dir, "genes.json"), sep="")
# #genes = genes.to_dict("records")
# data = {
#     "metadata": metadata,
#     "genes": genes,
# }
# data = json.dumps(data, indent=4)
# app_dir = os.path.dirname(os.path.abspath(__file__))
# app_dir = os.path.abspath(os.path.join(app_dir, "../..", "build"))
# index_template=os.path.join(app_dir, "index.html")
# template = Template(open(index_template,"r").read())
# html = template.render(data=data)
# output_html=os.path.join(app_dir,"Timeseries.html")
# output = open(output_html,"w")
# output.write(html)
# output.close()
# datalake_build=os.path.join(sys.argv[1], "build")
# shutil.copytree(build_folder, datalake_build)

if __name__ == "__main__":

    title = sys.argv[2]
    filename = sys.argv[1]

    print(filename)
    print(os.listdir())

    adata = sc.read(filename)

    genes = get_gene_data(adata)
    genes = genes.reset_index(drop=True)
    load_df(genes, title.lower())

    metadata = get_metadata(adata).to_dict(orient='records')
    data = {
        "metadata": metadata
    }
    data = json.dumps(data, indent=4)

    # app_dir = os.path.abspath(os.path.join(app_dir, "../..", "build"))
    index_template=os.path.join("build", "index.html")
    template = Template(open(index_template,"r").read())


    html = template.render(data=data, dashboard_id=json.dumps(title), api_url=json.dumps('http://localhost:9200'))
    output_html=os.path.join("build", "timeseries.html")
    output = open(output_html,"w")

    js_txt = open(os.path.join("build", "main.js"), 'r').read()
    css_txt = open(os.path.join("build", "main.css"), 'r').read()

    html = html.replace('<script src="./main.js"></script>', f"<script>{js_txt}</script>")
    html = html.replace('<link href="./main.css" rel="stylesheet">', f"<style>{css_txt}</style>")

    output.write(html)
    output.close()