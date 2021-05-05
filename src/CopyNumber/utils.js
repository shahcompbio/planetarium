import * as d3 from "d3";

const MAX_STATE = 11;
export const Y_AXIS_WIDTH = 25;
export const X_AXIS_HEIGHT = 12;
export const TOP_PADDING = 5;

const COPY_NUMBER_COLOURS = [
  "#2e7aab",
  "#9ECAE1",
  "#CCCCCC",
  "#FDCC8A",
  "#FC8D59",
  "#E34A33",
  "#B30000",
  "#980043",
  "#DD1C77",
  "#DF65B0",
  "#C994C7",
  "#D4B9DA",
];

export const BACKGROUND_COLORS = ["#fefefe", "#eee"];

export const getCopyNumberArray = (maxState) =>
  Array.from(Array(maxState + 1).keys());

export const getGenomeYScale = (maxState, height) =>
  d3.scaleLinear().domain([-0.5, maxState]).range([height, 0]);

const COPY_NUMBER_ARRAY = getCopyNumberArray(MAX_STATE);

export const colorScale = d3
  .scaleLinear()
  .domain(COPY_NUMBER_ARRAY)
  .range(COPY_NUMBER_COLOURS);
