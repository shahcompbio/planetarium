import React, { useState, useRef, useEffect } from "react";
import PropTypes from "prop-types";

import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@material-ui/core";
import { useCanvas } from "../utils/useCanvas";
import Legend from "../Legend/Vertical";
import useLasso from "./utils/useLasso";

const PADDING = 10;

const NUM_LEGEND_WIDTH = 70;
const CAT_LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;

const AXIS_FONT = "normal 10px Noto Sans";
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

export const drawAxis = ({ context, xPos, yPos, xLabel, yLabel }) => {
  context.beginPath();
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
  context.restore();
};

export const drawPoints = ({
  context,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  idParam,
  highlightIDs,
  colorScale,
}) => {
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;

  data.forEach((point) => {
    context.fillStyle =
      point.hasOwnProperty(subsetParam) && highlightIDs.includes(point[idParam])
        ? colorScale(point[subsetParam])
        : NULL_POINT_COLOR;

    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      POINT_RADIUS,
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });

  // this draws highlightedPoints on top
  if (highlightIDs && highlightIDs.length !== data.length) {
    data
      .filter((datum) => highlightIDs.includes(datum[idParam]))
      .forEach((point) => {
        context.fillStyle =
          point.hasOwnProperty(subsetParam) &&
          highlightIDs.includes(point[idParam])
            ? colorScale(point[subsetParam])
            : NULL_POINT_COLOR;

        context.beginPath();
        context.arc(
          xScale(point[xParam]),
          yScale(point[yParam]),
          POINT_RADIUS,
          0,
          Math.PI * 2,
          true
        );
        context.fill();
      });
  }
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

const UMAP = ({
  data,
  width = 500,
  height = 500,
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  highlightIDs = null,
  colorScale = null,
  disable = false,
  fontFamily = null,
  labels = (value) => value,
  onLasso = (data) => {},
  onLegendHover = (value) => {},
  onLegendClick = (value) => {},
}) => {
  const isCategorical =
    typeof data.filter((datum) => datum.hasOwnProperty(subsetParam))[0][
      subsetParam
    ] === "string";

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

  const subsetColors =
    colorScale || getColorScale({ data, subsetParam, isCategorical });

  const [lassoData, drawLasso, addLassoHandler, resetLasso] = useLasso(
    data,
    xScale,
    yScale,
    xParam,
    yParam
  );

  const prevHighlightRef = useRef();

  useEffect(() => {
    // resetLasso if highlightIDs is suddenly null and lassoData still exists
    if (highlightIDs === null && prevHighlightRef.current !== null) {
      if (lassoData !== null) {
        resetLasso();
      }
    }
  }, [highlightIDs, lassoData]);

  useEffect(() => {
    prevHighlightRef.current = highlightIDs;
  }, [highlightIDs]);

  const [hoveredLegend, setHoveredLegend] = useState(null);
  const [clickedLegend, setClickedLegend] = useState(null);

  const legendIDs = hoveredLegend || clickedLegend;

  const subsettedIDs =
    lassoData !== null
      ? lassoData.map((datum) => datum[idParam])
      : legendIDs !== null
      ? legendIDs
      : highlightIDs !== null
      ? highlightIDs
      : data.map((datum) => datum[idParam]);

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      drawAxis({
        context,
        xPos: AXIS_SPACE / 2,
        yPos: canvasHeight - AXIS_SPACE / 2,
        xLabel: xParam,
        yLabel: yParam,
      });

      drawPoints({
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        idParam,
        highlightIDs: subsettedIDs,
        colorScale: subsetColors,
      });

      drawLasso(context);
      const disableLasso = disable || legendIDs !== null;
      addLassoHandler(canvas, disableLasso, onLasso);
    },
    canvasWidth,
    canvasHeight,
    [data, disable, subsetParam, subsettedIDs]
  );

  const legendFilter = isCategorical
    ? (value, datum) => datum[subsetParam] === value
    : (value, datum) =>
        datum.hasOwnProperty(subsetParam) && datum[subsetParam] >= value[0];

  const getLegendData = (value) => {
    return value === null
      ? value
      : data
          .filter((datum) => legendFilter(value, datum))
          .map((datum) => datum[idParam]);
  };

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid item style={{ paddingLeft: "40px" }}>
        <Legend
          fontFamily={fontFamily}
          width={legendWidth}
          height={height / 2}
          colorScale={subsetColors}
          ticks={
            isCategorical
              ? subsetColors
                  .domain()
                  .sort()
                  .map((value) => ({ value, label: labels(value) }))
              : 10
          }
          onHover={(value) => {
            const legendData = getLegendData(value);
            setHoveredLegend(legendData);
            onLegendHover(value);
          }}
          onClick={(value) => {
            const legendData = getLegendData(value);
            setClickedLegend(legendData);
            onLegendClick(value);
          }}
          disable={disable || lassoData !== null}
          reset={highlightIDs !== null}
        />
      </Grid>
    </Grid>
  );
};

UMAP.propTypes = {
  /**
   * Data
   */
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  /**
   * width of plot
   */
  width: PropTypes.number,
  /**
   * height of plot
   */
  height: PropTypes.number,
  /**
   * name of key for x axis in UMAP
   */
  xParam: PropTypes.string.isRequired,
  /**
   * name of key for y axis in UMAP
   */
  yParam: PropTypes.string.isRequired,
  /**
   * name of key for colouring UMAP
   */
  subsetParam: PropTypes.string.isRequired,
  /**
   * name of ID key for UMAP points
   */
  idParam: PropTypes.string,
  /**
   * IDs of points to highlight in UMAP
   */
  highlightIDs: PropTypes.arrayOf(PropTypes.string),
  /**
   * Color Scale for points
   */
  colorScale: PropTypes.func,
  /**
   * Whether to disable interactions on plot
   */
  disable: PropTypes.bool,
  /**
   * Handler for lassoing UMAP
   */
  onLasso: PropTypes.func,
  /**
   * Handler for hovering on legend
   */
  onLegendHover: PropTypes.func,
  /**
   * Handler for clicking on legend
   */
  onLegendClick: PropTypes.func,
};

export default UMAP;
