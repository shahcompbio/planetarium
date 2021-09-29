import React from "react";
import _ from "lodash";
import * as d3 from "d3";

import Legend from "./Horizontal";

const CATEGORICAL_COLORS = [
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
  "#ed1a1a",
];

const numericalColorScale = d3
  .scaleSequential(d3.interpolateViridis)
  .domain([0, 10]);

const Template = (args) => <Legend {...args} />;

export default {
  title: "Components/Legend/Horizontal",
  component: Legend,
};

export const Categorical = Template.bind({});
Categorical.args = {
  width: 150,
  title: "Test",
  ticks: [1, 2, 3, 4],
  colorScale: (value) => CATEGORICAL_COLORS[value],
};

export const CategoricalNoTitle = Template.bind({});
CategoricalNoTitle.args = {
  width: 150,
  ticks: [1, 2, 3, 4],
  colorScale: (value) => CATEGORICAL_COLORS[value],
};

export const CategoricalLabels = Template.bind({});
CategoricalLabels.args = {
  width: 150,
  ticks: [1, 2, 3, 4].map((value) => ({ value, label: `Label ${value}` })),
  colorScale: (value) => CATEGORICAL_COLORS[value],
};

export const CategoricalNoInteraction = Template.bind({});
CategoricalNoInteraction.args = {
  width: 150,
  ticks: [1, 2, 3, 4],
  colorScale: (value) => CATEGORICAL_COLORS[value],
  disable: true,
};

export const Font = Template.bind({});
Font.args = {
  width: 150,
  title: "Test",
  ticks: [1, 2, 3, 4],
  colorScale: (value) => CATEGORICAL_COLORS[value],
  fontFamily: { regular: "Noto Sans", bold: "Noto Sans" },
};

export const Numerical = Template.bind({});
Numerical.args = {
  width: 150,
  title: "Test",
  colorScale: numericalColorScale,
  ticks: 3,
};

export const NumericalTicks = Template.bind({});
NumericalTicks.args = {
  width: 150,
  ticks: 5,
  title: "Test",
  colorScale: numericalColorScale,
};

export const NumericalNoTitle = Template.bind({});
NumericalNoTitle.args = {
  width: 150,
  ticks: 2,
  colorScale: numericalColorScale,
};

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

export const CopyNumber = Template.bind({});
CopyNumber.args = {
  ticks: COPY_NUMBER_COLOURS.map((value, index) => ({
    value: index,
    label: index === COPY_NUMBER_COLOURS.length - 1 ? `≥${index}` : index,
  })),
  colorScale: (value) => COPY_NUMBER_COLOURS[value],
  title: "Copy Number",
};
