import React from "react";
import App from "./Timeseries/Timeseries";
import fetchFileData from "./Timeseries/data/api";

const DevApp = () => {
  const data = fetchFileData();

  return Object.keys(data).length === 0 ? null : <App data={data} />;
};

const ProdApp = () => <App data={window.isablData} />;

export default process.env.NODE_ENV === "development" ? DevApp : ProdApp;
