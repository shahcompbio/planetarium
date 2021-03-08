import React from "react";
import ReactDOM from "react-dom";
import "./index.css";
import App from "./App";
import * as serviceWorker from "./serviceWorker";

import { ApolloProvider } from "react-apollo";
import client from "./apollo.js";
//const data = process.env.NODE_ENV === "development" ? [] : window.isablData;
import data from "./data2.js";

ReactDOM.render(
  <ApolloProvider client={client}>
    <App data={data} />
  </ApolloProvider>,
  document.getElementById("root")
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister();
