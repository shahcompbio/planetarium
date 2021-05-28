FROM python:3.8 AS builder

COPY requirements.txt .

# install dependencies
RUN pip install -r requirements.txt

COPY ./docker ./docker
COPY ./build ./build

CMD python ./docker/convert.py --file=./data/ndv.h5ad

#append data into static HTML
CMD python ./docker/render.py

FROM nginx:alpine
COPY --from=builder ./build /usr/share/nginx/html
CMD ["nginx", "-g", "daemon off;"]
