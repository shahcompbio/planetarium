from jinja2 import Template
import subprocess
import pandas as pd
import sys
import os
import json
import shutil


data_dir = sys.argv[1]

metadata = pd.read_csv(os.path.join(data_dir, "metadata.tsv"), sep="\t")
metadata = metadata.to_dict('records')

probabilities = pd.read_csv(os.path.join(data_dir, "probabilities.tsv"), sep="\t")
probabilities = probabilities.to_dict("records")

degs = pd.read_csv(os.path.join(data_dir, "degs.tsv"), sep="\t")
degs = degs.to_dict("records")

data = {
    "metadata": metadata,
    "probabilities": probabilities,
    "degs": degs
}

data = json.dumps(data, indent=4)

app_dir = os.path.dirname(os.path.abspath(__file__))
app_dir = os.path.abspath(os.path.join(app_dir, "../..", "build"))
index_template=os.path.join(app_dir, "index.html")
template = Template(open(index_template,"r").read())


html = template.render(data=data)
output_html=os.path.join(app_dir,"NDV.html")
output = open(output_html,"w")
output.write(html)
output.close()
# datalake_build=os.path.join(sys.argv[1], "build")
# shutil.copytree(build_folder, datalake_build)
