from flask_session import Session
from flask import Flask, session, jsonify, request
from flask_cors import CORS, cross_origin
import redis
import tqdm
import json
import sys
from ete3 import Tree
import random

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})

app.secret_key = "renderTree"
server_session = Session(app)

@app.route('/render/<path:file>/',methods=['GET', 'POST'])
def render(file=None):
    if file:
        #print(file)
        t = Tree("./"+file)
        #print(t)
        jsonTree = get_json(t)
        with open('../data/test.json', 'w') as json_file:
            json.dump([jsonTree], json_file)

    return "hello"


def get_json(node):
  # Read ETE tag for duplication or speciation events
  if not hasattr(node, 'evoltype'):
    dup = random.sample(['N','Y'], 1)[0]
  elif node.evoltype == "S":
    dup = "N"
  elif node.evoltype == "D":
    dup = "Y"

  node.name = node.name.replace("'", '')
  print(node)
  json = { "name": node.name,
                     "duplication": dup,
                     "branch_length": str(node.dist),
                     "type": "node" if node.children else "leaf",
                     "uniprot_name": "Unknown",
                     }

  if node.children:
    json["children"] = []
    for ch in node.children:
      json["children"].append(get_json(ch))

  return json
