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
  const [selectedSubtype, setSelectedSubtype] = useState(null);
  const [selectedClonotype, setSelectedClonotype] = useState(null);

  const { metadata, probabilities, degs } = data;
  const target = useRef(null);

  const topTen = Object.entries(
    _.countBy(metadata.map(row => row[initialState["clonotypeParam"]]))
  )
    .sort(([, a], [, b]) => b - a)
    .slice(0, 10);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = metadata.filter(row =>
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
          <div ref={target}>Click me!</div>

          <Popup />

          <div style={{ display: "flex" }}>
            <Layout
              chartName={"UMAP"}
              data={metadata}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
            <Layout
              chartName={"SUBTYPEUMAP"}
              data={metadata}
              selectedSubtype={selectedSubtype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
          </div>
          <div style={{ display: "flex" }}>
            <Layout
              chartName={"HEATMAP"}
              dim={{
                chart: {
                  x1: 50,
                  x2: 500,
                  y1: 100,
                  y2: 500
                },
                height: 500,
                width: 750
              }}
              data={probabilities}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
            <Layout
              chartName={"TABLE"}
              data={degs}
              selectedSubtype={selectedSubtype}
              dim={{
                chart: {
                  x1: 50,
                  y1: 50,
                  x2: 600,
                  y2: 200
                },
                height: 400,
                width: 650
              }}
            />
          </div>
          <div style={{ display: "flex" }}>
            <Layout
              chartName={"HISTOGRAM"}
              data={probabilities}
              dim={{
                chart: {
                  x1: 100,
                  y1: 50,
                  x2: 600,
                  y2: 200
                },
                height: 300,
                width: 650
              }}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
              setSelectedClonotype={clonotype =>
                setSelectedClonotype(clonotype)
              }
            />
            <Layout
              chartName={"BARPLOT"}
              data={probabilities}
              dim={{
                chart: {
                  x1: 50,
                  y1: 50,
                  x2: 600,
                  y2: 200
                },
                height: 300,
                width: 600
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
const Popup = props => (
  <div class="fixed-top">
    <div class="card" style="width: 18rem;">
      <div class="card-header">Featured</div>
      <ul class="list-group list-group-flush">
        <li class="list-group-item">Cras justo odio</li>
        <li class="list-group-item">Dapibus ac facilisis in</li>
        <li class="list-group-item">Vestibulum at eros</li>
      </ul>
    </div>
  </div>
);
export default NDV;
