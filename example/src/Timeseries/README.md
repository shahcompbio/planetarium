## start docker
```
docker-compose up -d
```

## to load into ES

```
docker run -v <path to data>:/usr/src/app/ndv.h5ad --network='host' timeseries <ID of dashboard>
```

## to build HTML

```
yarn build

python render.py <ID of dashboard> <path to data file>
```