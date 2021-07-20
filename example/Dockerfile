FROM python:3.6
WORKDIR /usr/src/app

COPY src/Timeseries/requirements.txt .
RUN pip install -r requirements.txt
COPY src/Timeseries/render.py .
COPY src/Timeseries/data/ensemble_ids.tsv .
COPY build ./build

ENTRYPOINT ["python", "render.py", "data.h5ad"]