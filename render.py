from jinja2 import Template
import subprocess
import pandas
import sys
import os
import json


#tsv = sys.argv[1]
df = pandas.read_csv("./vdj.tsv",sep="\t")

result = df.to_json(orient="records")
parsed = json.loads(result)

new_parsed = []
#for p in parsed:
#    p["sample"] = "NDVL"
#    new_parsed.append(p)*/

for _ in range(0):
    new_parsed = new_parsed + new_parsed
data = json.dumps(parsed, indent=4)

print(data)
template = Template(open("./build/index.html","r").read())

html = template.render(data=parsed)
output = open("./build/tcells.html","w")
output.write(html)
output.close()
