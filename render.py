from jinja2 import Template
import subprocess
import pandas
import sys
import os
import json
import shutil


app_dir = os.path.dirname(os.path.abspath(__file__))
build_folder = os.path.join(app_dir, "build")

# tsv = sys.argv[2]
# df = pandas.read_csv(tsv,sep="\t")


# result = df.to_json(orient="records")
with open(sys.argv[1]) as file:
    parsed = json.load(file)

data = json.dumps(parsed, indent=4)

index_template=os.path.join(app_dir,"build/index.html")
template = Template(open(index_template,"r").read())

html = template.render(data=data)
output_html=os.path.join(app_dir,"build/tcells.html")
output = open(output_html,"w")
output.write(html)
output.close()
# datalake_build=os.path.join(sys.argv[1], "build")
# shutil.copytree(build_folder, datalake_build)
