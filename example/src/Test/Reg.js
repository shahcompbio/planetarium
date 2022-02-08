import React, { useState, useRef, useEffect } from "react";
import PropTypes from "prop-types";

import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@mui/material";
import { useCanvas, useLasso } from "@shahlab/planetarium";
import regl, { ReglFrame, Regl } from "react-regl";
//import Circle from "./Circle.js";

const PADDING = 10;

const NUM_LEGEND_WIDTH = 70;
const CAT_LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;

const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";
const AXIS_LENGTH = 50;

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

const UMAP = regl({
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
      return props.data.map((d) => JSON.parse(d["color"]));
    },
    position: function (context, props, i) {
      return props.data.map((d) => [d.x, d.y]);
    },
  },

  uniforms: {
    stageWidth: 1000,

    stageHeight: 800,

    pointWidth: (context, props) => 2,
    color: function (context, props) {
      return props.data.map((d) => JSON.parse(d["color"]));
    },
  },

  count: function (context, props) {
    return props.data.length;
  },

  primitive: "points",
});

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

const Reg = ({
  data,
  width = 500,
  height = 500,
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  onLasso = (data) => {},
  yScale,
  xScale,
  canvasRef,
}) => {
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

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <ReglFrame canvasRef={canvasRef} onFrame={(context, regl) => {}}>
          <UMAP data={newData} />
        </ReglFrame>
      </Grid>
      <Grid item style={{ paddingLeft: "40px" }}></Grid>
    </Grid>
  );
};

export default Reg;
