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
        print(file)
        adata = sc.read("./"+file)
        time = []
        patientsList = []
        print(adata.obs.columns)
        patient = adata.obs[["timepoint","patient","clone_id"]]

        final = []

        for index, row in patient.iterrows():
            clone = row["clone_id"]
            timepoint = row["timepoint"]
            patient = row["patient"]

            b = {"clone":clone, "timepoint":timepoint}
            copy = b.copy()
            final.append(copy)

        with open('../data/test.json', 'w') as json_file:
           json.dump(final, json_file)

        print("fin")
        return "hello"
