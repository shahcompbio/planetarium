import React from "react";
import regl, { ReglFrame, Regl } from "react-regl";
import * as d3 from "d3";

import { useLasso } from "@shahlab/planetarium";

import _ from "lodash";

// Some constants to use
var MAX_WIDTH = 200;
var MAX_HEIGHT = 800;
var MAX_SPEED = 25;
var POINT_SIZE = 10;
var POINT_COUNT = 100000;
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
/*const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];*/
const PADDING = 10;
const AXIS_SPACE = 20;

const AXIS_LENGTH = 50;

const Circle = regl({
  frag: `
  precision mediump float;
  //attribute vec4 color;
  //uniform vec4 color;
  varying vec4 fcolor;
  void main () {
    gl_FragColor = fcolor;
  }`,

  vert: `
  precision mediump float;
  attribute vec2 position;
  attribute vec4 color;
  varying vec4 fcolor;

  // @change acquire the pointWidth uniform
  //  this is set by the uniforms section below
  uniform float pointWidth;
  uniform float stageWidth;
  uniform float stageHeight;

  vec2 normalizeCoords(vec2 position) {
  // read in the positions into x and y vars
  float x = position[0];
  float y = position[1];
  return vec2(
    2.0 * ((x / stageWidth) - 0.5),
    // invert y to treat [0,0] as bottom left in pixel space
    -(2.0 * ((y / stageHeight) - 0.5)));
  }

  void main () {
    // @change Set gl_PointSize global to
    //  configure point size
    fcolor = color;
    gl_PointSize = pointWidth;
    gl_Position = vec4(normalizeCoords(position), 0, 1);
  }`,

  attributes: {
    color: function (context, props) {
      return props.data.map((d) => {
        return JSON.parse(d["color"]);
      });
    },
    position: function (context, props, i) {
      return props.data.map((d) => {
        return [d.x, d.y];
      });
    },
  },

  uniforms: {
    stageWidth: 1000,

    stageHeight: 800,
    // @change: Add a pointWidth uniform -
    //  set by a prop
    pointWidth: (context, props) => 2,
    color: function (context, props) {
      return props.data.map((d) => {
        return JSON.parse(d["color"]);
      });
    },
  },

  count: function (context, props) {
    return props.data.length;
  },
  // @change: Set our primitive to points
  primitive: "points",
});

function onFrame(context, regl) {
  regl.clear({
    color: [1, 1, 1, 1],
  });
}

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

const Reg2 = ({ data, xParam, yParam, subsetParam, width, height }) => {
  const chartWidth = width;
  const chartHeight = height;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([PADDING, PADDING + chartWidth]);

  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([PADDING, PADDING + chartHeight]);

  const subsetColors = getColorScale({
    data,
    subsetParam,
    isCategorical: true,
  });

  const newData = data.map((d) => {
    const x = xScale(d[xParam]);
    const y = yScale(d[yParam]);
    return { ...d, x: x, y: y, color: subsetColors(d[subsetParam]) };
  });

  const [lassoData, drawLasso, addLassoHandler, resetLasso] = useLasso(
    data,
    xScale,
    yScale,
    xParam,
    yParam
  );

  return (
    <ReglFrame width={width} height={height} onFrame={onFrame}>
      <Circle data={newData} />
    </ReglFrame>
  );
};
export default Reg2;
