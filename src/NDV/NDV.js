import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";
import Modal from "react-bootstrap/Modal";
import OverlayTrigger from "react-bootstrap/OverlayTrigger";
import Overlay from "react-bootstrap/Overlay";
import Button from "react-bootstrap/Button";

const NDV = ({ data }) => {
  const [selectedSubtype, setSelectedSubtype] = useState(
    initialState["defaultSelectedObject"]
  );
  const [selectedClonotype, setSelectedClonotype] = useState(
    initialState["defaultSelectedObject"]
  );

  const { metadata, probabilities, degs } = data;
  const target = useRef(null);

  const filteredMetadata = metadata.filter(
    row => row[initialState["clonotypeParam"]] !== "None"
  );
  const topTen = Object.entries(
    _.countBy(filteredMetadata.map(row => row[initialState["clonotypeParam"]]))
  )
    .sort(([, a], [, b]) => b - a)
    .slice(0, 10);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = metadata.filter(
    row =>
      sampleTen.hasOwnProperty(row[initialState["clonotypeParam"]]) &&
      row[initialState["clonotypeParam"]] !== "None"
  );
  const topTenNumbering = Object.keys(sampleTen).reduce((final, seq, index) => {
    final[seq] = "SEQ" + (index + 1);
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
          sampleTen: sampleTen,
          sampleData: sampleData,
          topTenNumbering: topTenNumbering,
          topTen: topTen,
          colors: colors,
          clonotypes: clonotypes
        }}
        reducer={dashboardReducer}
      >
        <div>
          {(selectedClonotype["selected"] || selectedSubtype["selected"]) && (
            <Popup
              selected={
                selectedClonotype["selected"] || selectedSubtype["selected"]
              }
              setSelected={() => {
                setSelectedClonotype(initialState["defaultSelectedObject"]);
                setSelectedSubtype(initialState["defaultSelectedObject"]);
              }}
              type={selectedClonotype["selected"] ? "Clonotype" : "Subtype"}
            />
          )}
          <div style={{ display: "flex" }}>
            <Layout
              chartName={"UMAP"}
              data={metadata}
              selectedClonotype={selectedClonotype["selected"]}
              hoveredClonotype={selectedClonotype["hover"]}
              setSelectedClonotype={clonotype => {
                if (clonotype["selected"]) {
                  setSelectedSubtype(initialState["defaultSelectedObject"]);
                }
                setSelectedClonotype({ ...clonotype });
              }}
            />
            <Layout
              chartName={"SUBTYPEUMAP"}
              data={metadata}
              selectedSubtype={selectedSubtype["selected"]}
              hoveredSubtype={selectedSubtype["hover"]}
              setSelectedSubtype={subtype => {
                if (subtype["selected"]) {
                  setSelectedClonotype(initialState["defaultSelectedObject"]);
                }
                setSelectedSubtype({ ...subtype });
              }}
            />
          </div>
          <div style={{ display: "flex" }}>
            <Layout
              chartName={"HEATMAP"}
              dim={{
                chart: {
                  x1: 30,
                  x2: 500,
                  y1: 100,
                  y2: 500
                },
                height: 500,
                width: 750
              }}
              data={probabilities}
              selectedSubtype={
                selectedSubtype["selected"] || selectedSubtype["hover"]
              }
              selectedClonotype={
                selectedClonotype["selected"] || selectedClonotype["hover"]
              }
            />
            <Layout
              chartName={"TABLE"}
              data={degs}
              selectedSubtype={
                selectedSubtype["selected"]
                  ? selectedSubtype["selected"]
                  : selectedSubtype["hover"]
              }
              dim={{
                chart: {
                  x1: 50,
                  y1: 50,
                  x2: 600,
                  y2: 400
                },
                height: 500,
                width: 750
              }}
            />
          </div>
          <div style={{ display: "flex" }}>
            <Layout
              chartName={"BARPLOT"}
              data={probabilities}
              dim={{
                chart: {
                  x1: 30,
                  x2: 500,
                  y1: 50,
                  y2: 400
                },
                height: 475,
                width: 700
              }}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
            <Layout
              chartName={"HISTOGRAM"}
              data={probabilities}
              dim={{
                chart: {
                  x1: 100,
                  y1: 50,
                  x2: 600,
                  y2: 400
                },
                height: 500,
                width: 750
              }}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
          </div>
        </div>
      </DashboardProvider>
    </div>
  );
};
const Popup = ({ selected, setSelected, type }) => (
  <div
    class="fixed-top"
    style={{
      width: 150,
      float: "right",
      right: "10px",
      top: "10px",
      left: "auto"
    }}
  >
    <div class="card">
      <div class="card-header">Selected {type}:</div>
      <ul class="list-group list-group-flush">
        <li class="list-group-item">
          {selected}
          <button
            type="button"
            class="close"
            aria-label="Close"
            onClick={setSelected}
          >
            <span aria-hidden="true">&times;</span>
          </button>
        </li>
      </ul>
    </div>
  </div>
);
export default NDV;
