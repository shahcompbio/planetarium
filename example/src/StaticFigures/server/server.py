from flask import Flask, session, jsonify, request, make_response
from flask_session import Session
from datetime import timedelta
from flask_cors import CORS, cross_origin
import anndata
from flask_cors import CORS
import redis
import scanpy as sc
from scipy.sparse import csr_matrix, find
import pandas as pd
import tqdm
import json
import statistics
import numpy as np
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smt
from sklearn.preprocessing import MinMaxScaler

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_PERMANENT'] = True
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_COOKIE_NAME'] = "permanent_cookie_1"
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=False)
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=5)

app.secret_key = "render"
server_session = Session(app)

@app.route('/render/<path:file>/',methods=['GET', 'POST'])
def render(file=None):
    if file:
        adata = sc.read("./"+file)
        time = []
        patientsList = []
        #print(adata.obs.patient)

        #for d in adata.obs.index.tolist():
        #for d in adata.obs["sample"].tolist():
        #    if  "Pre" in d:
        #        time.append("Pre")
        #    else:
        #        time.append("Post")

        #    if "AE" in d:
        #        patientsList.append("AE")
        #    elif "NP" in d:
        #        patientsList.append("NP")
        #    elif "BM" in d:
        #        patientsList.append("BM")
        #    elif "MK" in d:
        #        patientsList.append("MK")

        #adata.obs["patient"] = patientsList
        #adata.obs["timepoint"] = time
        print('hello')
        print(adata.obs["timepoint"])
        timepoint = sc.tl.rank_genes_groups(adata, "timepoint")

        sc.tl.rank_genes_groups(adata,"timepoint")

        logfc = list(zip(*adata.uns['rank_genes_groups']["logfoldchanges"].tolist()))
        genes = list(zip(*adata.uns['rank_genes_groups']["names"].tolist()))
        adjpvals = list(zip(*adata.uns['rank_genes_groups']["pvals_adj"].tolist()))
        order = ["Pre","Post"]

        accepted = json.loads(request.data)
        accepted = accepted["data"]
        PATIENTS = list(set(adata.obs["patient"]))
        timepoint = sc.tl.rank_genes_groups(adata,"timepoint")

        final = {"All":{}, 'AE':{}, 'BM':{}, 'MK':{}, 'NP':{}}

        prescaled = []
        acceptedlist = accepted.split(",")

        #for g in acceptedlist:
        #    pre = adata[adata.obs["timepoint"]=="Pre"]
        #    post = adata[adata.obs["timepoint"]=="Post"]
        #    premean = np.mean(pre.X[:,pre.var.index.tolist().index(g)])
        #    postmean = np.mean(post.X[:,post.var.index.tolist().index(g)])
        #    prescaled.append([premean,postmean])
        #scaler = MinMaxScaler()
        #scaled = scaler.fit_transform(prescaled)


        for cond,fcs,gene,pvalues in zip(order,logfc,genes,adjpvals):
            for fc, g, pval in zip(fcs,gene,pvalues):
                if g in acceptedlist:
                    pre = adata[adata.obs["timepoint"]=="Pre"]
                    post = adata[adata.obs["timepoint"]=="Post"]
                    premean = np.mean(pre.X[:,pre.var.index.tolist().index(g)])
                    postmean = np.mean(post.X[:,post.var.index.tolist().index(g)])
                    if str(fc) == "nan":
                        if premean < postmean:
                            fc = "*+"
                        if premean > postmean:
                            fc = "*-"
                    index = acceptedlist.index(g)
                    #mean = scaled[index]
                    final["All"][g] = {"Pre":str(premean),"Post":str(postmean),"fc":str(fc),"p":str(pval)}

        #f = dict.fromkeys(accepted, {})
        #print(f)

        for patient in PATIENTS:
            final[patient] = {}
            print(adata.obs["patient"])
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

            prescaled = []
            indexTrack = {}
            i = 0

            for idx, g in enumerate(acceptedlist):
                if g in pre.var.index.tolist() and g in post.var.index.tolist():
                    premean = np.mean(pre.X[:,pre.var.index.tolist().index(g)])
                    postmean = np.mean(post.X[:,post.var.index.tolist().index(g)])
                    prescaled.append([premean,postmean])
                    indexTrack[g] = i
                    i = i + 1
            print(indexTrack)
            scaler = MinMaxScaler()
            scaled = scaler.fit_transform(prescaled)
            print(scaled)

            for cond,fcs,gene,pvalues in zip(order,logfc,genes,adjpvals):
                for fc, g, pval in zip(fcs,gene,pvalues):
                    if g in acceptedlist and g in indexTrack:
                        #premean = np.mean(pre.X[:,pre.var.index.tolist().index(g)])
                        #postmean = np.mean(post.X[:,post.var.index.tolist().index(g)])
                        if str(fc) == "nan":
                            if premean < postmean:
                                fc = "*+"
                            if premean > postmean:
                                fc = "*-"
                        index = indexTrack[g]
                        #index = acceptedlist.index(g)
                        mean = scaled[index]
                        final[patient][g] = {"Pre":str(mean[0]),"Post":str(mean[1]),"fc":str(fc),"p":str(pval)}
                        #final[patient][g] = {"Pre":str(premean),"Post":str(postmean),"fc":str(fc),"p":str(pval)}


        with open('../data/test.json', 'w') as json_file:
           json.dump([final], json_file)

        print("fin")
        return "hello"
