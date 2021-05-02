import React from "react";
import _ from "lodash";

import Legend from "./VerticalLegend";

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
  title: "Components/Legend/Vertical",
  component: Legend,
};

export const Vertical = Template.bind({});
Vertical.args = {
  width: 150,
  height: 400,
  labels: [1, 2, 3, 4].map((value) => ({
    value,
    label: `Label ${value}`,
    color: COLORS[value],
  })),
  setHighlighted: () => {},
};
