import React from "react";

import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@mui/material";
import { ReglFrame } from "react-regl";

import WebglUMAP from "./utils/drawWebglUMAP";

const COLOR_ARRAY = [
  "[0.369,0.31,0.635,1.0]",
  "[0.196,0.533,0.741,1.0]",
  "[0.4,0.761,0.647,1.0]",
  "[0.996,0.878,0.545,1.0]",
  "[0.957,0.427,0.263,1.0]",
  "[0.835,0.243,0.31,1.0]",
  "[0.788,0.8,0.463,1.0]",
  "[0.62,0.004,0.259,1.0]",
];

const getColorScale = ({ data, subsetParam, isCategorical }) => {
  if (isCategorical) {
    const subsetGroups = _.groupBy(data, subsetParam);
    const subsetValues = Object.keys(subsetGroups).sort();
    return d3
      .scaleOrdinal()
      .domain(subsetValues)
      .range(
        COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
      );
  } else {
    const subsetData = data
      .filter((d) => d.hasOwnProperty(subsetParam))
      .map((d) => parseFloat(d[subsetParam]));

    const subsetMax = Math.max(...subsetData);
    return d3
      .scaleSequential(d3.interpolateViridis)
      .domain([0, subsetMax])
      .nice();
  }
};

const ReglUmap = ({
  data,
  width = 500,
  height = 500,
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  onLasso = (data) => {},
  labels = (value) => value,
  yScale,
  xScale,
  canvasRef,
  isCategorical = true,
  pointSize = 2,
}) => {
  const subsetColors = getColorScale({
    data,
    subsetParam,
    isCategorical: isCategorical,
  });

  const newData = data.map((d) => {
    const x = xScale(d[xParam]);
    const y = yScale(d[yParam]);
    return { ...d, x: x, y: y, color: subsetColors(d[subsetParam]) };
  });

  return (
    <div style={{ padding: 0, position: "absolute" }}>
      <ReglFrame canvasRef={canvasRef} onFrame={(context, regl) => {}}>
        <WebglUMAP
          data={newData}
          stageWidth={width}
          stageHeight={height}
          pointSize={pointSize}
        />
      </ReglFrame>
    </div>
  );
};

export default ReglUmap;
