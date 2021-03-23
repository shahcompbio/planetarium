/*

Histogram with kde curves

*/

import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useCanvas } from "../components/utils/useCanvas";

import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";

import { canvasInit, drawAxis, changeFontSize } from "../DrawingUtils/utils.js";

const HIGHLIGHTED_BAR_COLOR = "#eb5067";
const HIGHLIGHTED_BAR_WIDTH = 2;

const PADDING = 10;
const TITLE_HEIGHT = 30;
const X_AXIS_HEIGHT = 50;
const Y_AXIS_WIDTH = 50;

const DataWrapper = ({
  chartName,
  data,
  chartDim,
  selectedSubtype,
  selectedClonotype,
}) => {
  const [{ logXParam, clonotypeParam, subtypeParam }] = useDashboardState();
  return (
    <Histogram
      data={data}
      binParam={logXParam}
      lineParam={subtypeParam}
      highlightBarParam={clonotypeParam}
      chartDim={chartDim}
      chartName={chartName}
      highlightedBar={
        selectedClonotype["hover"] || selectedClonotype["selected"]
      }
      highlightedLine={selectedSubtype["hover"] || selectedSubtype["selected"]}
    />
  );
};

const Histogram = ({
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
    .thresholds(x.ticks(20))(data);

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

  context.font = "normal 10px Helvetica";

  x.ticks(10).forEach((tick) => {
    context.fillText(tick, x(tick), y(0) + 15);
  });

  const format = (tick) => (tick === 0 ? "0" : d3.format(".2f")(tick / maxBin));

  y.ticks(10).map((tick) => {
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
  context.font = "normal 12px Helvetica";
  context.textAlign = "center";
  context.textBaseline = "hanging";
  context.fillText(binParam, chartWidth / 2, chartHeight + PADDING * 2);

  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText("Probability", -(chartHeight / 2), PADDING);
  context.restore();
};

const drawBars = (context, bins, x, y, barScale) => {
  context.globalAlpha = 1;
  context.fillStyle = "#6bb9f0";
  context.strokeStyle = "#5c97bf";

  bins.forEach((bin) => {
    const xPos = x(bin["x0"]);
    const yPos = y(bin.length);
    const width = x(bin["x1"]) - x(bin["x0"]);
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
    x.ticks(20),
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
  context.strokeStyle = "steelblue";
  context.stroke();
};

const Histogram2 = ({ chartName, data, chartDim, highlighted }) => {
  const [
    { logXParam, logYParam, fontSize, clonotypeParam },
  ] = useDashboardState();
  const [drawReady, setDrawReady] = useState(false);
  const [context, saveContext] = useState(null);

  const allX = data.map((row) => parseFloat(row[logXParam]));
  const allY = data.map((row) => parseFloat(row[logYParam]));
  const xMax = Math.max(...allX);
  const xMin = Math.min(...allX);

  const x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"] - 50]);

  const bins = d3Array
    .bin()
    .value((d) => d[logXParam])
    .domain(x.domain())
    .thresholds(x.ticks(10));

  const buckets = bins(data);
  const maxY = Math.max(...buckets.map((row) => row.length));

  // Y axis
  const y = d3
    .scaleLinear()
    .domain([0, maxY / data.length])
    .range([chartDim["chart"]["y2"], chartDim["chart"]["y1"]]);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      drawAxisLabels(context);
      drawBars(context, data);
      drawKde(context, data);
    },
    chartDim["width"],
    chartDim["width"],
    [highlighted]
  );
  function drawKde(context, data) {
    const density = kde(
      epanechnikov(1),
      x.ticks(20),
      data.map((row) => parseFloat(row["log10_probability"]))
    );
    var line = d3
      .line()
      .curve(d3.curveBasis)
      .x(function(d) {
        return x(d[0]);
      })
      .y(function(d) {
        return y(d[1]);
      })
      .context(context);

    context.beginPath();
    line(density);
    context.lineWidth = 2;
    context.strokeStyle = "steelblue";
    context.stroke();
  }
  const kde = (kernel, thresholds, data) =>
    thresholds.map((t) => [t, d3.mean(data, (d) => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return (x) =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }

  function drawBars(context, data) {
    context.fillStyle = "#6bb9f0";
    context.strokeStyle = "#5c97bf";
    buckets.map((bar) => {
      const yPos = y(bar.length / data.length);
      context.fillRect(
        x(bar["x0"]) - 5,
        yPos,
        x(bar["x1"]) - x(bar["x0"]) - 5,
        y(0) - yPos
      );
      context.strokeRect(
        x(bar["x0"]) - 5,
        yPos,
        x(bar["x1"]) - x(bar["x0"]) - 5,
        y(0) - yPos
      );
      context.fill();
    });

    if (highlighted) {
      const highlightedData = data.filter(
        (datum) => datum[clonotypeParam] === highlighted
      );

      if (highlightedData.length > 0) {
        const highlightedProbability = highlightedData[0][logXParam];

        context.fillStyle = HIGHLIGHTED_BAR_COLOR;
        context.fillRect(
          x(highlightedProbability),
          y(maxY / data.length),
          HIGHLIGHTED_BAR_WIDTH,
          y(0) - y(maxY / data.length)
        );
      }
    }
  }

  function drawAxisLabels(context) {
    const format = (tick) => (tick === 0 ? "0" : d3.format(".2f")(tick));
    context.beginPath();
    context.globalAlpha = 1;
    //context.lineWidth = 1;
    context.fillStyle = "black";
    context.textAlign = "right";
    const xMinValue = x(xMin) - 10;
    changeFontSize(context, fontSize["tickLabelFontSize"]);
    x.ticks(10).map((tick) => {
      context.fillText(tick, x(tick), y(0) + 15);
    });
    changeFontSize(context, fontSize["axisLabelFontSize"]);
    context.fillText(
      logXParam,
      (chartDim["chart"]["x2"] - chartDim["chart"]["x1"]) / 2 +
        chartDim["chart"]["x1"],
      chartDim["chart"]["y2"] + 50
    );

    context.restore();
    y.ticks(10).map((tick) => {
      context.globalAlpha = 1;
      changeFontSize(context, fontSize["tickLabelFontSize"]);
      context.fillText(format(tick), xMinValue, y(tick) + 3);
      context.fill();
      context.save();
      //  context.strokeStyle = "rgba(0, 0, 0, 0.5)";
      context.globalAlpha = 0.2;
      context.lineWidth = 0.5;
      context.beginPath();
      context.moveTo(xMinValue + 3, y(tick));
      context.lineTo(x(xMax), y(tick));
      context.stroke();
      context.restore();
    });

    context.globalAlpha = 1;
    context.save();
    context.rotate((270 * Math.PI) / 180);
    changeFontSize(context, fontSize["axisLabelFontSize"]);
    context.fillText(
      "Probability",
      -(chartDim["chart"]["y2"] - chartDim["chart"]["y1"]) / 2 - 15,
      chartDim["chart"]["x1"] - 50
    );
    context.restore();
  }

  return (
    <Paper>
      <Grid
        container
        direction="row"
        justify="flex-start"
        alignItems="flex-start"
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative",
        }}
      >
        <Grid
          item
          xs={17}
          sm={8}
          id="histogram"
          style={{
            pointerEvents: "all",
          }}
        >
          <canvas ref={ref} />
        </Grid>
        <Grid
          item
          xs={7}
          sm={4}
          style={{
            textAlign: "right",
            marginTop: 10,
            paddingRight: 15,
            pointerEvents: "all",
            curser: "pointer",
            zIndex: 100,
          }}
        >
          {infoText[chartName]["title"] + "    "}
          <Info name={chartName} direction="s" />
        </Grid>
      </Grid>
    </Paper>
  );
};
export default DataWrapper;
