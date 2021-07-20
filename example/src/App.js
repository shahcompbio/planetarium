import React from "react";
import App from "./CITESEQ/CITESEQ";
import fetchFileData from "./CITESEQ/data/api";

const DevApp = () => {
  const data = fetchFileData();

  return Object.keys(data).length === 0 ? null : (
    <App data={data} dashboardID={"Test"} api={"http://localhost:9200"} />
  );
};

const ProdApp = () => (
  <App
    data={window.isablData}
    dashboardID={window.dashboardID}
    api={window.apiURL}
  />
);

export default process.env.NODE_ENV === "development" ? DevApp : ProdApp;
