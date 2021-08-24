import React from "react";
import * as d3 from "d3";

import numData from "./data/cells.json";
import { useCanvas } from "../utils/useCanvas";
import { drawAxis, drawPoints } from "./UMAP";

const CanvasComponent = ({ canvasFunc, ...args }) => {
  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      canvasFunc({ context, ...args });
    },
    500,
    500,
    []
  );

  return <canvas ref={canvasRef} />;
};

export default {
  title: "Components/UMAP/Utils",
  component: CanvasComponent,
};

const AxisTemplate = (args) => (
  <CanvasComponent canvasFunc={drawAxis} {...args} />
);

export const Axis = AxisTemplate.bind({});
Axis.args = {
  xPos: 50,
  yPos: 450,
  xLabel: "UMAP 1",
  yLabel: "UMAP 2",
};

const PointsTemplate = (args) => (
  <CanvasComponent canvasFunc={drawPoints} {...args} />
);

const yData = numData.map((d) => parseFloat(d["rna_UMAP_2"]));
const xData = numData.map((d) => parseFloat(d["rna_UMAP_1"]));

const yMin = Math.min(...yData);
const yMax = Math.max(...yData);
const xMin = Math.min(...xData);
const xMax = Math.max(...xData);

const xScale = d3.scaleLinear().domain([xMin, xMax]).range([10, 490]);

const yScale = d3.scaleLinear().domain([yMax, yMin]).range([10, 490]);

const subsetData = numData
  .filter((d) => d.hasOwnProperty("VIM"))
  .map((d) => parseFloat(d["VIM"]));

const subsetMax = Math.max(...subsetData);
const colorScale = d3
  .scaleSequential(d3.interpolateViridis)
  .domain([0, subsetMax])
  .nice();

export const Points = PointsTemplate.bind({});
Points.args = {
  data: numData,
  xParam: "rna_UMAP_1",
  yParam: "rna_UMAP_2",
  xScale,
  yScale,
  subsetParam: "VIM",
  highlightIDs: numData.map((datum) => datum["cell_id"]),
  idParam: "cell_id",
  colorScale,
};

export const PointsSubset = PointsTemplate.bind({});
PointsSubset.args = {
  data: numData,
  xParam: "rna_UMAP_1",
  yParam: "rna_UMAP_2",
  xScale,
  yScale,
  subsetParam: "VIM",
  highlightIDs: numData
    .filter((datum) => datum.hasOwnProperty("VIM") && datum["VIM"] >= 2)
    .map((datum) => datum["cell_id"]),
  idParam: "cell_id",
  colorScale,
};
