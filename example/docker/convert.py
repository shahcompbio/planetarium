import scanpy as sc
import numpy as np
import pandas as pd
import math
import argparse

parser = argparse.ArgumentParser(description='Convert h5ad file to tsv')
parser.add_argument('--file', type=str, help='Name of h5ad file to convert')
args = parser.parse_args()

file=args.file
#import os
#ile = os.environ['FILE']
print(file)
adata = sc.read(file)
umap1, umap2 = list(zip(*adata.obsm['X_umap']))
adata.obs = adata.obs.assign(UMAP_1=umap1, UMAP_2=umap2 )
final_frame = adata.obs[['UMAP_1', 'UMAP_2', 'IR_VDJ_1_cdr3', 'IR_VDJ_1_v_gene', 'subtype']]
final_frame.rename(columns={'IR_VDJ_1_v_gene':'v_gene', 'IR_VDJ_1_cdr3': 'cdr3s_aa'}, inplace=True)
final_frame.replace(to_replace="nan", value=np.nan, inplace=True)
final_frame.to_csv("metadata.tsv", index=True, index_label="cell_id", sep="\t", na_rep="None")
genes = list(zip(*adata.uns['rank_genes_groups']["names"].tolist()))
subtypes = adata.uns['rank_genes_groups']["names"]
subtypes = subtypes.dtype.names
adjpvals = list(zip(*adata.uns['rank_genes_groups']["pvals_adj"].tolist()))
logfc = list(zip(*adata.uns['rank_genes_groups']["logfoldchanges"].tolist()))
def subset_significant(genes, pvalues, fcs, pthresh=0.01, fcthresh=0.01, top_n=20, sort_by=2):
    subset = []
    for gene, pvalue, fc in zip(genes,pvalues,fcs):
        if pvalue < pthresh and fc > fcthresh:
            subset.append((gene, pvalue, fc))
    subset.sort(key=lambda x:-1.0 *x[sort_by])
    return subset[:top_n]
df_genes   = []
df_pvals   = []
df_fcs     = []
df_subtype = []
for subtype, gene, pval, fc in zip(subtypes, genes, adjpvals, logfc):
    subset = subset_significant(gene, pval, fc)
    df_genes += [x[0] for x in subset]
    df_pvals += [x[1] for x in subset]
    df_fcs += [x[2] for x in subset]
    df_subtype += [subtype for _ in subset]
data = {"gene":df_genes, "adj_pval":df_pvals, "log_fc": df_fcs, "subtype":df_subtype}
df = pd.DataFrame.from_dict(data)
df.to_csv("degs.tsv", sep="\t", index=False)
pgen = adata.obs
pgen = pgen[pgen["IR_VDJ_1_cdr3"]!="None"]
pgen.index = pgen["IR_VDJ_1_cdr3"]
pgen = pgen[["pgen","subtype"]]
pgen = pgen[pgen["pgen"]!="None"]
pgen["pgen"] = [math.log10(float(prob)) for prob in pgen["pgen"].tolist()]
pgen.to_csv("probabilities.tsv", sep="\t", index=True)
