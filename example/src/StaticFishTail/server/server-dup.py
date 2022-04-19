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
#app.config['SESSION_TYPE'] = 'redis'
#app.config['SESSION_PERMANENT'] = True
#app.config['SESSION_USE_SIGNER'] = True
#app.config['SESSION_COOKIE_NAME'] = "permanent_cookie_1"
#app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
#app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=False)
#app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=5)

app.secret_key = "render1"
server_session = Session(app)

config = {"steve" :{"timepoint":"sample", "patient":"patient", "clone_id":"clone_phy"},
"normal" :{"timepoint":"timepoint", "patient":"patient", "clone_id":"clone_id",},
"peds":{"timepoint":"tp", "patient":"patient", "clone_id":"clone"},
};
type = "steve";

@app.route('/render/<path:file>/',methods=['GET', 'POST'])
def render(file=None):
    if file:
        print(file)
        adata = sc.read("./"+file)
        time = []
        patientsList = []

        print(adata.obsm)
        #print(adata.obs.columns)
        #print(config[type]["timepoint"])
        patient = adata.obs[[config[type]["timepoint"],config[type]["patient"],config[type]["clone_id"]]]

        final = []

        for index, row in patient.iterrows():
            clone = row[config[type]["clone_id"]]
            timepoint = row[config[type]["timepoint"]]
            patient = row[config[type]["patient"]]

            b = {"clone":clone, "timepoint":timepoint}
            copy = b.copy()
            final.append(copy)

        with open('../data/test.json', 'w') as json_file:
           json.dump(final, json_file)

        print("fin")
        return "hello"
#gene = "MHCII"
gene = "ERBB2"
@app.route('/render2/<path:file>/',methods=['GET', 'POST'])
def render2(file=None):
    if file:
        print(file)
        adata = sc.read("./"+file)
        time = []
        patientsList = []
        cloneList = []
        cloneNameList = []
        print(adata.var_keys)
        print(adata.obs[config[type]["clone_id"]])
        #for clone in set(adata.obs["clone"].tolist()):
             #if clone.startswith("AE"):
        #        new_clone_pre = adata[adata.obs["clone"] == clone]
        #        new_clone_pre = new_clone_pre[new_clone_pre.obs["tp"] == "Diagnosis"]

        #        new_clone_post = adata[adata.obs["clone"] == clone]
        #        new_clone_post =  new_clone_post[new_clone_post.obs["tp"] != "Diagnosis"]

        #        cloneList.append(new_clone_pre)
        #        cloneList.append(new_clone_post)
        #        cloneNameList.append(clone+"-pre")
        #        cloneNameList.append(clone+"-post")

        #print(adata.obs["clone_phy"])
        #geneIndex = 0
        for g in adata.var.index.tolist():
           if g == gene:
               for clone in set(adata.obs[config[type]["clone_id"]].tolist()):
                   #print(adata.obs[config[type]["timepoint"]])
                   for timepoint in set(adata.obs[config[type]["timepoint"]].tolist()):
                       new_clone_pre = adata[adata.obs["clone"] == clone]
                       new_clone_pre = new_clone_pre[new_clone_pre.obs[config[type]["timepoint"]] == timepoint]
                       cloneList.append(new_clone_pre)
                       cloneNameList.append(clone+"-"+timepoint)
                       clonemean = np.mean(new_clone_pre.X[:,new_clone_pre.var.index.tolist().index(g)])
                       print(clonemean)
                    #        new_clone_post = adata[adata.obs["clone"] == clone]
                    #        new_clone_post =  new_clone_post[new_clone_post.obs["tp"] != "Diagnosis"]

                    #        cloneList.append(new_clone_pre)
                    #        cloneList.append(new_clone_post)
                    #        cloneNameList.append(clone+"-pre")
                    #        cloneNameList.append(clone+"-post")
                    #print(clone)
                    #clone1mean = np.mean(clone1.X[:,clone1.var.index.tolist().index(g)])
                    #clone2mean = np.mean(clone2.X[:,clone2.var.index.tolist().index(g)])
        clones= {}
        print(cloneNameList)
        for idx, clone in enumerate(cloneList):
            clone_mean = np.mean(clone.X[:,clone.var.index.tolist().index(gene)])
            clone_name = cloneNameList[idx]
            clones[clone_name]=clone_mean
            print(clone_name,clone_mean)
        #for idx, clone in enumerate(cloneList):
        #    clone_mean = np.mean(clone.obs["MHCII"])
        #    clone_name = cloneNameList[idx]
        #    clones[clone_name]=clone_mean
        #    print(clone_name,clone_mean)



        #print(adata.obs.columns)
        #print(config[type]["timepoint"])
        patient = adata.obs[[config[type]["timepoint"],config[type]["patient"],config[type]["clone_id"]]]

        final = []

        for index, row in patient.iterrows():
            clone = row[config[type]["clone_id"]]
            timepoint = row[config[type]["timepoint"]]
            patient = row[config[type]["patient"]]

            b = {"clone":clone, "timepoint":timepoint}
            copy = b.copy()
            final.append(copy)

        with open('../data/test1.json', 'w') as json_file:
            finalfinal = {"data":final, "clones":clones}
            json.dump(finalfinal, json_file)

        print("fin")
        return "hello"
