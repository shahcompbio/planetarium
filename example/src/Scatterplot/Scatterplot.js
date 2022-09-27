import React, { useState, useRef, useEffect } from "react";
import PropTypes from "prop-types";

import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@mui/material";
import { Legend, useD3 } from "@shahlab/planetarium";
const PADDING = 40;

const NUM_LEGEND_WIDTH = 70;
const CAT_LEGEND_WIDTH = 180;
const AXIS_SPACE = 40;

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;

const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";
const AXIS_LENGTH = 50;

const COLOR_ARRAY = [
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
];

export const drawAxis = ({ svgRef, xPos, yPos, xScale, yScale }) => {
  console.log(yPos);
  svgRef
    .append("g")
    .attr("transform", "translate(" + xPos + ",0)") // This controls the vertical position of the Axis
    .call(d3.axisLeft(yScale));
  svgRef
    .append("g")
    .attr("transform", "translate(0," + yPos + ")") // This controls the vertical position of the Axis
    .call(d3.axisBottom(xScale));

  /*  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.lineWidth = 1;
  context.lineCap = "butt";
  context.moveTo(xPos, yPos);
  context.lineTo(xPos, yPos - AXIS_LENGTH);
  context.stroke();

  context.beginPath();
  context.moveTo(xPos, yPos);
  context.lineTo(xPos + AXIS_LENGTH, yPos);
  context.stroke();

  context.textAlign = "left";
  context.textBaseline = "middle";
  context.fillText(xLabel, xPos + AXIS_LENGTH + 2, yPos);
  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText(yLabel, -(yPos - AXIS_LENGTH - 2), xPos);
  context.restore();*/
};

export const drawPoints = ({
  svgRef,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  idParam,
  highlightIDs,
  colorScale,
  pointSize,
}) => {
  svgRef
    .selectAll("circle")
    .data(data)
    .enter()
    .append("circle")
    .attr("cx", function (d) {
      return xScale(d[xParam]);
    })
    .attr("cy", function (d) {
      return yScale(d[yParam]);
    })
    .attr("r", pointSize)
    .attr("fill", (d) => colorScale(parseFloat(d[subsetParam])))
    .attr("opacity", 0.8);
};

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

const Scatterplot = ({
  data,
  width = 600,
  height = 600,
  xParam = "Entropy",
  yParam = "Size",
  subsetParam = "Exhaustion",
  idParam = "id",
  highlightIDs = null,
  colorScale = null,
  pointSize = 5,
}) => {
  console.log(data);
  const isCategorical = false;

  const legendWidth = isCategorical ? CAT_LEGEND_WIDTH : NUM_LEGEND_WIDTH;

  const canvasWidth = width - legendWidth;
  const canvasHeight = height;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

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

  const subsetColors = getColorScale({ data, subsetParam, isCategorical });

  const prevHighlightRef = useRef();

  const [hoveredLegend, setHoveredLegend] = useState(null);
  const [clickedLegend, setClickedLegend] = useState(null);

  const legendIDs = hoveredLegend || clickedLegend;

  const subsettedIDs =
    legendIDs !== null
      ? legendIDs
      : highlightIDs !== null
      ? highlightIDs
      : data.map((datum) => datum[idParam]);

  const canvasRef = useD3(
    (svgRef) => {
      requestAnimationFrame(() => {
        drawPoints({
          svgRef,
          data,
          xScale,
          yScale,
          xParam,
          yParam,
          subsetParam,
          idParam,
          highlightIDs: subsettedIDs,
          colorScale: subsetColors,
          pointSize,
        });

        drawAxis({
          svgRef,
          xPos: AXIS_SPACE,
          yPos: canvasHeight - AXIS_SPACE * 2,
          xLabel: xParam,
          yLabel: yParam,
          xScale,
          yScale,
        });
      });
    },
    canvasWidth,
    canvasHeight,
    [data]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <svg ref={canvasRef} />
      </Grid>
      <Grid item style={{ paddingLeft: "40px" }}></Grid>
    </Grid>
  );
};
/*  <Legend
    width={legendWidth}
    height={height / 2}
    colorScale={subsetColors}
    /*  ticks={
      isCategorical
        ? subsetColors
            .domain()
            .sort()
            .map((value) => ({ value, label: labels(value) }))
        : 10
    }
    onHover={(value) => {
      const legendData = getLegendData(value);
      //  setHoveredLegend(legendData);
      //  onLegendHover(value);
    }}
    onClick={(value) => {
      const legendData = getLegendData(value);
      //  setClickedLegend(legendData);
      //  onLegendClick(value);
    }}
    //  disable={disable}
    reset={highlightIDs !== null}
  />*/
export default Scatterplot;
