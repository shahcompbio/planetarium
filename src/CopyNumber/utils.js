import * as d3 from "d3";

export const X_AXIS_HEIGHT = 12;

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

export const colorScale = d3
  .scaleOrdinal()
  .domain(Array.from(Array(COPY_NUMBER_COLOURS.length).keys()))
  .range(COPY_NUMBER_COLOURS);
