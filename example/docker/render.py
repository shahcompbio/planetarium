from jinja2 import Template
import subprocess
import pandas as pd
import sys
import os
import json
import shutil
#import argparse
import os

#parser = argparse.ArgumentParser(description='Convert h5ad file to tsv')
#parser.add_argument('--file', type=str, help='Name of h5ad file to convert')
#parser.add_argument('--output', type=str, help='Name of tsv output file')
#parser.add_argument('--project', type=str, help='Name of project (options: vdjVisual, vdjTimeSeries)')
#args = parser.parse_args()

data_dir = "./docker"

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
app_dir = os.path.abspath(os.path.join("build"))
index_template=os.path.join(app_dir, "index.html")
template = Template(open(index_template,"r").read())


html = template.render(data=data)
output_html=os.path.join(app_dir,"index.html")
print(output_html)
output = open(output_html,"w")
output.write(html)
output.close()
