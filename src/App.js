import React from "react";
import logo from "./logo.svg";
import Layout from "./Layout/Layout.js";
import "./App.css";
import _ from "lodash";
import * as d3 from "d3";

import dashboardReducer, {
  initialState
} from "./PlotState/dashboardReducer.js";
import { DashboardProvider } from "./PlotState/dashboardState.js";

const App = ({ data }) => {
  //data manipulations
  const topTen = Object.entries(
    _.countBy(data.map(row => row[initialState["clonotypeParam"]]))
  )
    .sort(([, a], [, b]) => b - a)
    .slice(0, 10);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = data.filter(row =>
    sampleTen.hasOwnProperty(row[initialState["clonotypeParam"]])
  );
  const topTenNumbering = Object.keys(sampleTen).reduce((final, seq, index) => {
    const label = "NDVL" === "NDVL" ? "L" + (index + 1) : "R" + (index + 1);
    final[seq] = label;
    return final;
  }, {});
  const clonotypes = _.groupBy(sampleData, initialState["clonotypeParam"]);

  const types = Object.keys(clonotypes);

  const colourList = [
    "#674172",
    "#098dde",
    "#fa832f",
    "#0e5702",
    "#c20c1e",
    "#911eb4",
    "#fc97bc",
    "#469990",
    "#b5762a",
    "#5aebed",
    "#8f8f3f",
    "#ed1a1a"
  ];
  var colors = d3
    .scaleOrdinal()
    .domain([...types])
    .range([...colourList]);

  return (
    <div className="App">
      <DashboardProvider
        initialState={{
          ...initialState,
          sampleTen,
          sampleData,
          topTenNumbering,
          topTen,
          colors,
          clonotypes
        }}
        reducer={dashboardReducer}
      >
        <div className="App">
          <Layout chartName={"UMAP"} data={data} />
          <Layout chartName={"HEATMAP"} data={data} />
        </div>
      </DashboardProvider>
    </div>
  );
};

export default App;
