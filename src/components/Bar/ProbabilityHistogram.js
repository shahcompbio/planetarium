/*

Probability distribution with kde curves

*/

import React from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useCanvas } from "../utils/useCanvas";

import Info from "../../Info/Info.js";
import infoText from "../../Info/InfoText.js";

const HIGHLIGHTED_BAR_COLOR = "#eb5067";
const HIGHLIGHTED_BAR_WIDTH = 2;

const PADDING = 10;
const TITLE_HEIGHT = 30;
const X_AXIS_HEIGHT = 50;
const Y_AXIS_WIDTH = 50;

const NUM_TICKS = 25;
const LABEL_FONT = "normal 12px Helvetica";
const TICK_FONT = "normal 10px Helvetica";

const BAR_COLOR = "#6bb9f0";
const BAR_STROKE_COLOR = "#5c97bf";

const LINE_COLOR = "steelblue";

const ProbabilityHistogram = ({
  data,
  binParam,
  lineParam,
  highlightBarParam,
  chartDim,
  highlightedBar,
  highlightedLine,
  chartName,
}) => {
  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - PADDING - PADDING - TITLE_HEIGHT;

  const chartWidth = canvasWidth - Y_AXIS_WIDTH;
  const chartHeight = canvasHeight - X_AXIS_HEIGHT;

  const allX = data.map((row) => parseFloat(row[binParam]));
  const xMax = Math.max(...allX);
  const xMin = Math.min(...allX);

  const x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([Y_AXIS_WIDTH, Y_AXIS_WIDTH + chartWidth - PADDING]);

  const bins = d3Array
    .bin()
    .value((d) => d[binParam])
    .domain(x.domain())
    .thresholds(x.ticks(NUM_TICKS))(data);

  const maxY = Math.max(...bins.map((row) => row.length));

  const y = d3
    .scaleLinear()
    .domain([0, maxY])
    .range([chartHeight, PADDING]);

  const barScale = d3
    .scaleLinear()
    .domain([0, maxY])
    .range([0, chartHeight - PADDING]);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      drawAxisLabels(
        context,
        x,
        y,
        binParam,
        chartWidth,
        chartHeight,
        xMin,
        xMax,
        data.length
      );
      drawBars(context, bins, x, y, barScale);
      drawHighlightedBar(
        context,
        data,
        highlightedBar,
        highlightBarParam,
        binParam,
        maxY,
        x,
        y,
        barScale
      );
      drawKde(context, data, x, y, binParam, lineParam, highlightedLine);
    },
    canvasWidth,
    canvasHeight,
    [highlightedBar, highlightedLine]
  );

  return (
    <Paper
      style={{
        margin: 10,
        padding: PADDING,
        height: chartDim["height"],
        width: chartDim["width"],
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="stretch"
      >
        <Grid
          item
          style={{
            textAlign: "right",
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
        <Grid item>
          <canvas ref={ref} />
        </Grid>
      </Grid>
    </Paper>
  );
};

const drawAxisLabels = (
  context,
  x,
  y,
  binParam,
  chartWidth,
  chartHeight,
  xMin,
  xMax,
  maxBin
) => {
  context.beginPath();
  context.globalAlpha = 1;
  context.fillStyle = "black";
  context.textAlign = "right";

  context.font = TICK_FONT;

  x.ticks(10).forEach((tick) => {
    context.fillText(tick, x(tick), y(0) + 15);
  });

  const format = (tick) => (tick === 0 ? "0" : d3.format(".2f")(tick / maxBin));

  y.ticks(10).forEach((tick) => {
    context.globalAlpha = 1;
    context.textBaseline = "middle";
    context.fillText(format(tick), x(xMin), y(tick));
    context.globalAlpha = 0.2;
    context.lineWidth = 0.5;
    context.beginPath();
    context.moveTo(x(xMin), y(tick));
    context.lineTo(x(xMax), y(tick));
    context.stroke();
  });

  context.globalAlpha = 1;
  context.font = LABEL_FONT;
  context.textAlign = "center";
  context.textBaseline = "hanging";
  context.fillText(binParam, chartWidth / 2, chartHeight + 20);

  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText("Density", -(chartHeight / 2), 10);
  context.restore();
};

const drawBars = (context, bins, x, y, barScale) => {
  context.globalAlpha = 1;
  context.fillStyle = BAR_COLOR;
  context.strokeStyle = BAR_STROKE_COLOR;

  bins.forEach((bin) => {
    const xPos = x(bin["x0"]) + 1;
    const yPos = y(bin.length);
    const width = x(bin["x1"]) - x(bin["x0"]) - 2;
    const height = barScale(bin.length);
    context.fillRect(xPos, yPos, width, height);
    context.strokeRect(xPos, yPos, width, height);
  });
};

const drawHighlightedBar = (
  context,
  data,
  highlightedBar,
  barParam,
  binParam,
  maxY,
  x,
  y,
  barScale
) => {
  if (highlightedBar) {
    const highlightedData = data.filter(
      (datum) => datum[barParam] === highlightedBar
    );

    if (highlightedData.length > 0) {
      const highlightedX = highlightedData[0][binParam];

      context.fillStyle = HIGHLIGHTED_BAR_COLOR;
      context.fillRect(
        x(highlightedX),
        y(maxY),
        HIGHLIGHTED_BAR_WIDTH,
        barScale(maxY)
      );
    }
  }
};

const drawKde = (context, data, x, y, binParam, lineParam, highlightedLine) => {
  const kde = (kernel, thresholds, data) =>
    thresholds.map((t) => [t, d3.mean(data, (d) => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return (x) =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }
  const densityData = highlightedLine
    ? data.filter((datum) => datum[lineParam] === highlightedLine)
    : data;

  const density = kde(
    epanechnikov(1),
    x.ticks(NUM_TICKS),
    densityData.map((row) => parseFloat(row[binParam]))
  );
  var line = d3
    .line()
    .curve(d3.curveBasis)
    .x(function(d) {
      return x(d[0]);
    })
    .y(function(d) {
      return y(d[1] * densityData.length);
    })
    .context(context);

  context.beginPath();
  line(density);
  context.lineWidth = 2;
  context.strokeStyle = LINE_COLOR;
  context.stroke();
};

export default ProbabilityHistogram;
