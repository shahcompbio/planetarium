from elasticsearch import Elasticsearch
from elasticsearch import helpers
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, find
import sys

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


if __name__ == "__main__":
    title = sys.argv[2]
    filename = sys.argv[1]
    print(title)
    print(filename)

    adata = sc.read(filename)

    genes = get_gene_data(adata)
    genes = genes.reset_index(drop=True)
    load_df(genes, title.lower())