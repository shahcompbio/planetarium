import React from "react";
import { useCanvas } from "../utils/useCanvas";
import * as d3 from "d3";

const HEATMAP_COLOR_SCALE = d3
  .scaleLinear()
  .range(["#ffec8b", "#d91e18"])
  .domain([0, 96]);
const HEATMAP_ROW_SPACING = 3;

const LABEL_COLOR = "#000000";
const LABEL_FONT = "normal 12px Helvetica";

const HeatmapPlot = ({ data, highlighted, width, height }) => {
  const ref = useCanvas((canvas) => {}, width, height);

  return <canvas ref={ref} />;
};

export default HeatmapPlot;
