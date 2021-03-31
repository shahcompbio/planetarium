# Planetarium

A collection of visual components and dashboards for the lab.

Dashboards are outputted as HTML templates, where data is injected in through [Jinja](https://palletsprojects.com/p/jinja/).


## Installation

Install the dependencies:

```
yarn install
```


Then you can start the app in development mode.

```
yarn start
```

Open [http://localhost:3000](http://localhost:3000) to view it in the browser.


## Switching dashboards

To switch the dashboard that will be shown on development mode / build, go to `/src/App.js` and change the import to the correct dashboard

```
import App from '<!! path here>'
```

For example, for VDJ, we should change the import to the following

```
import App from './VDJ/VDJ
```

You should also ensure that you have test data where needed.


## Building

To build the app for production,

```
yarn build
```

This should output a `build` folder with the minified files. This is the HTML template and will need to be injected with data to create a standalone file.

To inject data, run the appropriate python script:

```
python3 render.py <path to data file>
```

