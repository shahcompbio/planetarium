FROM python:3.6
WORKDIR /usr/src/app

COPY requirements.txt .
RUN pip install -r requirements.txt
COPY render.py .
COPY load_es.py .
COPY data/ensemble_ids.tsv .

ENTRYPOINT ["python", "load_es.py", "ndv.h5ad"]