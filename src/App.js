import React from "react";
import logo from "./logo.svg";
import Layout from "./Layout/Layout.js";
import "./App.css";

const App = ({ data }) => {
  return (
    <div className="App">
      <Layout chartName={"UMAP"} data={data} />
    </div>
  );
};

export default App;
