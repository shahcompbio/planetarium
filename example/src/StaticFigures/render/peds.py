import scanpy as sc
# from scipy.sparse import csr_matrix, find
import pandas as pd
# import tqdm
import json
# import statistics
import numpy as np
# from scipy.stats import ttest_ind
# import statsmodels.stats.multitest as smt

# t cell
#adata = sc.read("./data/pediatrics_nonblasts.h5ad")
# protien version

PATIENTS = ['AE', 'BM', 'MK', 'NP']
ORDER = ["Post", "Pre"]
GENES = ["HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "HLA-DRA", "HLA-DQB1", "HLA-DMA", "CD74", "CIITA", "HLA-A",
         "HLA-DQA2", "HLA-C", "HLA-DQA1", "HLA-B", "HLA-E", "HLA-DOA", "HLA-DMB", "HLA-G", "HLA-DQB1-AS1", "HLA-F", "B2M"]

adata = sc.read("hacohen_viz.h5ad")
# protien version
#adata = sc.read("./data/protein_blasts.h5ad")

adata = adata[adata.obs["timepoint"] != "HD"]
#adata = adata[adata.obs["CellType"]=="T"]

# protein

f = dict.fromkeys(GENES, {})

final = {}

for patient in PATIENTS:
    final[patient] = {**f}
    adata_patient = adata[adata.obs["patient"] == patient]
    sc.tl.rank_genes_groups(adata_patient, "timepoint")
    logfc = list(
        zip(*adata_patient.uns['rank_genes_groups']["logfoldchanges"].tolist()))
    genes = list(
        zip(*adata_patient.uns['rank_genes_groups']["names"].tolist()))
    adjpvals = list(
        zip(*adata_patient.uns['rank_genes_groups']["pvals_adj"].tolist()))

    pre = adata_patient[adata_patient.obs["timepoint"] == "Pre"]
    post = adata_patient[adata_patient.obs["timepoint"] == "Post"]

    for cond, fcs, gene, pvalues in zip(ORDER, logfc, genes, adjpvals):
        for fc, g, pval in zip(fcs, gene, pvalues):
            if g in GENES:
                premean = np.mean(pre.X[:, pre.var.index.tolist().index(g)])
                postmean = np.mean(post.X[:, post.var.index.tolist().index(g)])
                if str(fc) == "nan":
                    if premean < postmean:
                        fc = "*+"
                    if premean > postmean:
                        fc = "*-"
                final[patient][g] = {"Pre": str(premean), "Post": str(
                    postmean), "fc": str(fc), "p": str(pval)}
                print(cond, fc, g, pval, premean, postmean)

#accepted = ["PDCD1","CD226","CTLA4","TIGIT","CD28","ICOS","ICAM1","TNFRSF4","LAG3","TIM3", '41BB', "CD44", "CD43", "CD40", "BLIMP1", "CD127", "CD62L", "CD69"]
#accepted = ["HLA-DRA",'LGALS9', "CEACAM1","HMGB1","LGALS3","CLEC4G","PVR","NECTIN2", "CD80", "CD86","CD274","PDCD1LG2","ICOSLG","CD276","TNFRSF14","CD137L","TNFSF4","CD70","CD40","TNFSF9","VSIR"]
# timepoint = sc.tl.rank_genes_groups(adata, "timepoint")
# f = dict.fromkeys(accepted, {})
# print(f)
# final = {"All": {}}
# final["All"] = f
# ​
# for cond, fcs, gene, pvalues in zip(order, logfc, genes, adjpvals):
#     for fc, g, pval in zip(fcs, gene, pvalues):
#         if g in accepted:
#             pre = adata[adata.obs["timepoint"] == "Pre"]
#             post = adata[adata.obs["timepoint"] == "Post"]
#             premean = np.mean(pre.X[:, pre.var.index.tolist().index(g)])
#             postmean = np.mean(post.X[:, post.var.index.tolist().index(g)])
#             if str(fc) == "nan":
#                 if premean < postmean:
#                     fc = "*+"
#                 if premean > postmean:
#                     fc = "*-"
#             final["All"][g] = {"Pre": str(premean), "Post": str(
#                 postmean), "fc": str(fc), "p": str(pval)}
#             print(cond, fc, g, pval, premean, postmean)

with open('test.json', 'w') as json_file:
    json.dump(str(final), json_file)

print("fin")
