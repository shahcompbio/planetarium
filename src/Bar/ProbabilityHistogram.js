/*

Probability distribution with kde curves

Data is an array of objects.

x-axis is the probability variable, and can be associated with an observation variable
y-axis is the density

The histogram bars can show the makeup of
The bars can be hovered and shown the breakdown of the density plot based off the subgroup variable name


*/

import React, { useState } from "react";
import PropTypes from "prop-types";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import * as _ from "lodash";
import Grid from "@material-ui/core/Grid";

import Tooltip from "../Tooltip/Tooltip";
import { useCanvas } from "../utils/useCanvas";
import drawAxis from "../utils/canvas/drawAxis";

const HIGHLIGHTED_BAR_COLOR = "#eb5067";
const HIGHLIGHTED_BAR_WIDTH = 2;

const X_AXIS_HEIGHT = 40;
const Y_AXIS_WIDTH = 50;

const NUM_TICKS = 25;

const BAR_COLOR = "#6bb9f0";
const BAR_STROKE_COLOR = "#5c97bf";
const HIGHLIGHTED_LINE_COLOR = "#47a647";

const LINE_COLOR = "steelblue";

const drawAxisLabels = (context, x, y, font) => {
  drawAxis({
    context,
    xScale: x,
    yScale: y,
    ticks: 10,
    orientation: "horizontal",
    label: "Log10pgen",
    gridlines: false,
    font: font,
  });

  const format = (tick) => (tick === 0 ? "0" : d3.format(".2f")(tick));
  drawAxis({
    context,
    xScale: x,
    yScale: y,
    ticks: 10,
    label: "Cell Density",
    format,
    font: font,
  });
};

const drawBars = (
  canvas,
  bins,
  x,
  y,
  barScale,
  total,
  setHighlightedBin,
  highlightedIDs,
  idParam
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
  if (highlightedIDs !== null) {
    context.fillStyle = HIGHLIGHTED_LINE_COLOR;
    bins.forEach((bin) => {
      const content =
        bin.filter((datum) => highlightedIDs.includes(datum[idParam])).length /
        total;
      const xPos = x(bin["x0"]) + 1;
      const yPos = y(content);
      const width = Math.max(x(bin["x1"]) - x(bin["x0"]) - 2, 0);
      const height = barScale(content);
      context.fillRect(xPos, yPos, width, height);
    });
  }

  // if (highlightedSubgroup) {
  //   context.fillStyle = HIGHLIGHTED_LINE_COLOR;
  //   bins.forEach((bin) => {
  //     const content =
  //       bin.filter((datum) => datum[subgroupParam] === highlightedSubgroup)
  //         .length / total;
  //     const xPos = x(bin["x0"]) + 1;
  //     const yPos = y(content);
  //     const width = Math.max(x(bin["x1"]) - x(bin["x0"]) - 2, 0);
  //     const height = barScale(content);
  //     context.fillRect(xPos, yPos, width, height);
  //   });
  // }

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

const drawKde = (context, data, x, y, probParam, idParam, highlightedIDs) => {
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

  if (highlightedIDs !== null) {
    const filteredDensity = kde(
      epanechnikov(1),
      x.ticks(NUM_TICKS * 2),
      data
        .filter((datum) => highlightedIDs.includes(datum[idParam]))
        .map((row) => parseFloat(row[probParam]))
    );
    context.beginPath();
    line(filteredDensity);
    context.lineWidth = 2;
    context.strokeStyle = HIGHLIGHTED_LINE_COLOR;
    context.stroke();
  }
};

const ProbabilityHistogram = ({
  data,
  width = 400,
  height = 400,
  probParam,
  observationParam,
  highlightedObservation = null,
  idParam = "id",
  highlightedIDs = null,
  getTooltipText = (bin) => `Count: ${bin.length}`,
  font = "MyFontLight",
}) => {
  const [highlightedBin, setHighlightedBin] = useState(null);

  const chartWidth = width - Y_AXIS_WIDTH;
  const chartHeight = height - X_AXIS_HEIGHT;

  const allX = data.map((row) => parseFloat(row[probParam]));
  const xDataMin = Math.min(...allX);
  const xDataMax = Math.max(...allX);
  const xDataWidth = Math.abs(xDataMax - xDataMin);
  const xMin = xDataMin - Math.abs(xDataWidth * 0.05);
  const xMax = xDataMax + Math.abs(xDataWidth * 0.05);

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
      drawAxisLabels(context, x, y, font);
      drawBars(
        canvas,
        bins,
        x,
        y,
        barScale,
        data.length,
        setHighlightedBin,
        highlightedIDs,
        idParam
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
      drawKde(context, data, x, y, probParam, idParam, highlightedIDs);
    },
    width,
    height,
    [highlightedObservation, highlightedIDs]
  );

  return (
    <Grid
      item
      style={{
        position: "relative",
        paddingLeft: 10,
      }}
    >
      <Tooltip
        getText={getTooltipText}
        getX={(bin) => x(bin["x0"]) + (x(bin["x1"]) - x(bin["x0"])) / 2}
        getY={(bin) => chartHeight - barScale(bin.length / data.length)}
        data={highlightedBin}
      />
      <canvas ref={ref} />
    </Grid>
  );
};

ProbabilityHistogram.propTypes = {
  /**
   * list of data points to bin
   */
  data: PropTypes.arrayOf(PropTypes.object).isRequired,
  /**
   * width of plot
   */
  width: PropTypes.number.isRequired,
  /**
   * height of plot
   */
  height: PropTypes.number.isRequired,
  /**
   * Key used to bin for histogram
   */
  probParam: PropTypes.string.isRequired,
  /**
   * Key used for observation
   */
  observationParam: PropTypes.string.isRequired,
  /**
   * Value associated with observationParam to highlight in plot
   */
  highlightedObservation: PropTypes.string,
  /**
   * Key used as ID
   */
  idParam: PropTypes.string,
  /**
   * List of IDs to highlight in plot
   */
  highlightedIDs: PropTypes.arrayOf(PropTypes.string),
  /**
   * Returns text to display when hover on bin
   */
  getTooltipText: PropTypes.func,
};

export default ProbabilityHistogram;
