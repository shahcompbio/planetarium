/*

Probability distribution with kde curves

Data is an array of objects.

x-axis is the probability variable, and can be associated with an observation variable
y-axis is the density

The histogram bars can show the makeup of
The bars can be hovered and shown the breakdown of the density plot based off the subgroup variable name


*/

import React, { useState } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import * as _ from "lodash";
import Grid from "@material-ui/core/Grid";

import { Tooltip } from "@material-ui/core";
import { useCanvas } from "../utils/useCanvas";

const HIGHLIGHTED_BAR_COLOR = "#eb5067";
const HIGHLIGHTED_BAR_WIDTH = 2;

const X_AXIS_HEIGHT = 40;
const Y_AXIS_WIDTH = 50;

const NUM_TICKS = 25;
const LABEL_FONT = "normal 12px Helvetica";
const TICK_FONT = "normal 10px Helvetica";

const BAR_COLOR = "#6bb9f0";
const BAR_STROKE_COLOR = "#5c97bf";
const HIGHLIGHTED_LINE_COLOR = "#47a647";

const LINE_COLOR = "steelblue";
const format = d3.format(".3f");

const drawAxisLabels = (
  context,
  x,
  y,
  probParam,
  chartWidth,
  chartHeight,
  xMin,
  xMax,
  xTickWidth
) => {
  context.beginPath();
  context.globalAlpha = 1;
  context.fillStyle = "black";
  context.textAlign = "right";

  context.font = TICK_FONT;

  context.textBaseline = "bottom";

  x.ticks(10).forEach((tick) => {
    context.fillText(tick, x(tick), y(0) + 15);
  });

  const format = (tick) => (tick === 0 ? "0" : d3.format(".2f")(tick));

  y.ticks(10).forEach((tick) => {
    context.globalAlpha = 1;
    context.textBaseline = "middle";
    context.fillText(format(tick), x(xMin) - xTickWidth, y(tick));
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
  context.fillText(probParam, chartWidth / 2, chartHeight + X_AXIS_HEIGHT / 2);

  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText("Density", -(chartHeight / 2), 5);
  context.restore();
};

const drawBars = (
  canvas,
  bins,
  x,
  y,
  barScale,
  total,
  setHighlightedBin,
  highlightedSubgroup,
  subgroupParam
) => {
  const context = canvas.getContext("2d");

  context.globalAlpha = 1;
  context.fillStyle = BAR_COLOR;
  context.strokeStyle = BAR_STROKE_COLOR;

  bins.forEach((bin) => {
    const xPos = x(bin["x0"]) + 1;
    const yPos = y(bin.length / total);
    const width = x(bin["x1"]) - x(bin["x0"]) - 2;
    const height = barScale(bin.length / total);
    context.fillRect(xPos, yPos, width, height);
    context.strokeRect(xPos, yPos, width, height);
  });

  if (highlightedSubgroup) {
    context.fillStyle = HIGHLIGHTED_LINE_COLOR;
    bins.forEach((bin) => {
      const content =
        bin.filter((datum) => datum[subgroupParam] === highlightedSubgroup)
          .length / total;
      const xPos = x(bin["x0"]) + 1;
      const yPos = y(content);
      const width = x(bin["x1"]) - x(bin["x0"]) - 2;
      const height = barScale(content);
      context.fillRect(xPos, yPos, width, height);
    });
  }

  d3.select(canvas)
    .on("mousemove", function () {
      const mouseX = d3.event.layerX || d3.event.offsetX;
      const [minX, maxX] = x.range();

      if (mouseX >= minX && mouseX <= maxX) {
        const bin = bins.filter(
          (bin) => x(bin["x0"]) <= mouseX && mouseX <= x(bin["x1"])
        )[0];
        if (bin.length > 0) {
          setHighlightedBin(bin);
        }
      }
    })
    .on("mouseout", () => setHighlightedBin(null));
};

const drawHighlightedBar = (
  context,
  data,
  highlightedObservation,
  barParam,
  probParam,
  maxY,
  x,
  y,
  barScale
) => {
  if (highlightedObservation) {
    const highlightedData = data.filter(
      (datum) => datum[barParam] === highlightedObservation
    );

    if (highlightedData.length > 0) {
      const highlightedX = highlightedData[0][probParam];

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

const getDensity = (data, x, probParam, ticks) => {
  const kde = (kernel, thresholds, data) =>
    thresholds.map((t) => [t, d3.mean(data, (d) => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return (x) =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }
  const density = kde(
    epanechnikov(1),
    x.ticks(ticks),
    data.map((row) => parseFloat(row[probParam]))
  );

  return density;
};

const drawKde = (
  context,
  data,
  x,
  y,
  probParam,
  subgroupParam,
  highlightedSubgroup
) => {
  const kde = (kernel, thresholds, data) =>
    thresholds.map((t) => [t, d3.mean(data, (d) => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return (x) =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }

  var line = d3
    .line()
    .curve(d3.curveBasis)
    .x(function (d) {
      return x(d[0]);
    })
    .y(function (d) {
      return y(d[1]);
    })
    .context(context);

  const allDensity = kde(
    epanechnikov(1),
    x.ticks(NUM_TICKS * 2),
    data.map((row) => parseFloat(row[probParam]))
  );
  context.beginPath();
  line(allDensity);
  context.lineWidth = 2;
  context.strokeStyle = LINE_COLOR;
  context.stroke();

  if (highlightedSubgroup) {
    const filteredDensity = kde(
      epanechnikov(1),
      x.ticks(NUM_TICKS * 2),
      data
        .filter((datum) => datum[subgroupParam] === highlightedSubgroup)
        .map((row) => parseFloat(row[probParam]))
    );
    context.beginPath();
    line(filteredDensity);
    context.lineWidth = 2;
    context.strokeStyle = HIGHLIGHTED_LINE_COLOR;
    context.stroke();
  }
};

const getTooltipText = (bin, subgroupParam) => {
  const subtypeGroups = _.groupBy(bin, subgroupParam);

  return (
    <div>
      {Object.keys(subtypeGroups)
        .sort((a, b) => a.length - b.length)
        .map((group) => (
          <p>
            {group}: {format(subtypeGroups[group].length / bin.length)}
          </p>
        ))}
    </div>
  );
};

const ProbabilityHistogram = ({
  data,
  width,
  height,
  probParam,
  subgroupParam,
  observationParam,
  highlightedObservation,
  highlightedSubgroup,
}) => {
  const [highlightedBin, setHighlightedBin] = useState(null);

  const chartWidth = width - Y_AXIS_WIDTH;
  const chartHeight = height - X_AXIS_HEIGHT;

  const allX = data.map((row) => parseFloat(row[probParam]));
  const xMax = Math.max(...allX);
  const xMin = Math.min(...allX);

  const startX = Y_AXIS_WIDTH;
  const x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([startX, startX + chartWidth]);

  const bins = d3Array
    .bin()
    .value((d) => d[probParam])
    .domain(x.domain())
    .thresholds(x.ticks(NUM_TICKS))(data);

  const maxYData = Math.max(...bins.map((row) => row.length / data.length));
  const density = getDensity(data, x, probParam, 1000);
  const maxYDensity = Math.max(...density.map((datum) => datum[1]));

  const maxY = Math.max(maxYData, maxYDensity) * 1.1;

  const y = d3.scaleLinear().domain([0, maxY]).range([chartHeight, 0]);

  const barScale = d3.scaleLinear().domain([0, maxY]).range([0, chartHeight]);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawAxisLabels(
        context,
        x,
        y,
        probParam,
        chartWidth,
        chartHeight,
        xMin,
        xMax,
        5
      );
      drawBars(
        canvas,
        bins,
        x,
        y,
        barScale,
        data.length,
        setHighlightedBin,
        highlightedSubgroup,
        subgroupParam
      );
      drawHighlightedBar(
        context,
        data,
        highlightedObservation,
        observationParam,
        probParam,
        maxY,
        x,
        y,
        barScale
      );
      drawKde(
        context,
        data,
        x,
        y,
        probParam,
        subgroupParam,
        highlightedSubgroup
      );
    },
    width,
    height,
    [highlightedObservation, highlightedSubgroup]
  );

  return (
    <Grid
      item
      style={{
        position: "relative",
      }}
    >
      <canvas ref={ref} />
      <Tooltip
        title={
          highlightedBin ? getTooltipText(highlightedBin, subgroupParam) : ""
        }
        open={highlightedBin !== null}
        arrow
        placement="top"
      >
        <div
          style={{
            position: "absolute",
            pointerEvents: "none",
            left: highlightedBin
              ? x(highlightedBin["x0"]) +
                (x(highlightedBin["x1"]) - x(highlightedBin["x0"])) / 2
              : null,
            top: highlightedBin
              ? chartHeight - barScale(highlightedBin.length / data.length)
              : null,
          }}
        />
      </Tooltip>
    </Grid>
  );
};

export default ProbabilityHistogram;
