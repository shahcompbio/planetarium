import React, { useState, useEffect } from "react";

import Layout from "./Layout.js";
import _ from "lodash";
import * as d3 from "d3";
import dashboardReducer, {
  initialState
} from "../PlotState/dashboardReducer.js";
import {
  DashboardProvider,
  useDashboardState
} from "../PlotState/dashboardState.js";

const LayoutWrapper = ({ data }) => {
  const [{}, dispatch] = useDashboardState();
  const [selectedSubtype, setSelectedSubtype] = useState(null);
  const [selectedClonotype, setSelectedClonotype] = useState(null);
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

  useEffect(() => {
    if (data.length > 0) {
      dispatch({
        type: "OVERRIDE",
        value: {
          sampleTen: sampleTen,
          sampleData: sampleData,
          topTenNumbering: topTenNumbering,
          topTen: topTen,
          colors: colors,
          clonotypes: clonotypes
        }
      });
    }
  }, [data]);
  return (
    <div>
      {data.length > 0 ? (
        <div>
          <Layout
            chartName={"SUBTYPEUMAP"}
            data={data}
            selectedSubtype={selectedSubtype}
            selectedClonotype={selectedClonotype}
            setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
            setSelectedClonotype={clonotype => setSelectedClonotype(clonotype)}
          />
          <Layout
            chartName={"UMAP"}
            data={data}
            selectedSubtype={selectedSubtype}
            selectedClonotype={selectedClonotype}
            setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
            setSelectedClonotype={clonotype => setSelectedClonotype(clonotype)}
          />
          <Layout
            chartName={"HEATMAP"}
            data={data}
            selectedSubtype={selectedSubtype}
            selectedClonotype={selectedClonotype}
            setSelectedSubtype={subtype => setSelectedSubtype(subtype)}
            setSelectedClonotype={clonotype => setSelectedClonotype(clonotype)}
          />
        </div>
      ) : null}
    </div>
  );
};
export default LayoutWrapper;
