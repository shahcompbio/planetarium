import React from "react";
import ReactDOM from "react-dom";
import "./index.css";
import App from "./App";
import * as serviceWorker from "./serviceWorker";

//const data = process.env.NODE_ENV === "development" ? [] : window.isablData;
import data from "./data/api.js";

ReactDOM.render(<App data={data} />, document.getElementById("root"));

serviceWorker.unregister();
