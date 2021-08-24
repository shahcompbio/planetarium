import React from "react";
import * as d3 from "d3";

import Legend from "./VerticalNumericalLegend";

const colorScale = d3.scaleSequential(d3.interpolateViridis).domain([0, 10]);

const Template = (args) => <Legend {...args} />;

export default {
  title: "Components/Legend/VerticalNumerical",
  component: Legend,
};

export const Default = Template.bind({});
Default.args = {
  width: 150,
  height: 500,
  title: "Test",
  colorScale,
};

export const MultiTick = Template.bind({});
MultiTick.args = {
  width: 150,
  height: 500,
  ticks: 5,
  title: "Test",
  colorScale,
};

export const Titleless = Template.bind({});
Titleless.args = {
  width: 150,
  height: 500,
  ticks: 2,
  colorScale,
};
