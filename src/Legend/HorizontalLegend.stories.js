import React from "react";
import _ from "lodash";

import Legend from "./HorizontalLegend";

const COLORS = [
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

const Template = (args) => <Legend {...args} />;

export default {
  title: "Components/Legend/Horizontal",
  component: Legend,
};

export const Horizontal = Template.bind({});
Horizontal.args = {
  labels: [1, 2, 3, 4].map((value) => ({
    value,
    label: `${value}`,
    color: COLORS[value],
  })),
  title: "Test Label",
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
  labels: COPY_NUMBER_COLOURS.map((value, index) => ({
    value: index,
    label: index === COPY_NUMBER_COLOURS.length - 1 ? `≥${index}` : index,
    color: value,
  })),
  title: "Copy Number",
};
