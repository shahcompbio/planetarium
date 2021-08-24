import React from "react";
import * as d3 from "d3";

import { useCanvas } from "../useCanvas";
import drawAxis from "./drawAxis";

const AxisComponent = (args) => {
  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawAxis({ context, ...args });
    },
    500,
    500,
    []
  );

  return <canvas ref={canvasRef} />;
};

const Template = (args) => <AxisComponent {...args} />;

const scale = d3.scaleLinear().range([50, 450]);

export default {
  title: "Utils/Canvas/Axis",
  component: AxisComponent,
};

export const Vertical = Template.bind({});
Vertical.args = {
  xScale: scale,
  yScale: scale,
  label: "Test label",
  ticks: 5,
};

export const Horizontal = Template.bind({});
Horizontal.args = {
  xScale: scale,
  yScale: scale,
  label: "Test label",
  ticks: 5,
  orientation: "horizontal",
};

export const Log = Template.bind({});
Log.args = {
  xScale: scale,
  yScale: d3.scaleLog().range([50, 450]),
  label: "Test label",
  ticks: 5,
};

export const Lineless = Template.bind({});
Lineless.args = {
  xScale: scale,
  yScale: scale,
  label: "Test label",
  ticks: 5,
  gridlines: false,
};
