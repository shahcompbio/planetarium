/*

Probability distribution with kde curves

*/

import React from "react";
import * as d3 from "d3";
import d3Tip from "d3-tip";
import * as d3Array from "d3-array";
import * as _ from "lodash";
import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useCanvas } from "../utils/useCanvas";

import Info from "../../Info/Info.js";
import infoText from "../../Info/InfoText.js";

const HIGHLIGHTED_BAR_COLOR = "#eb5067";
const HIGHLIGHTED_BAR_WIDTH = 2;

const PADDING = 10;
const TITLE_HEIGHT = 30;
const X_AXIS_HEIGHT = 20;
const Y_AXIS_WIDTH = 20;

const NUM_TICKS = 25;
const LABEL_FONT = "normal 12px Helvetica";
const TICK_FONT = "normal 10px Helvetica";

const BAR_COLOR = "#6bb9f0";
const BAR_STROKE_COLOR = "#5c97bf";

const LINE_COLOR = "steelblue";
const format = d3.format(".3f");

const ProbabilityHistogram = ({
  data,
  binParam,
  lineParam,
  highlightBarParam,
  chartDim,
  highlightedBar,
  highlightedLine,
  chartName
}) => {
  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - PADDING - PADDING - TITLE_HEIGHT;

  const chartWidth = canvasWidth - Y_AXIS_WIDTH;
  const chartHeight = canvasHeight - X_AXIS_HEIGHT;

  const allX = data.map(row => parseFloat(row[binParam]));
  const xMax = Math.max(...allX);
  const xMin = Math.min(...allX);
  const xTickSize = (xMax - xMin) / NUM_TICKS;

  const x = d3
    .scaleLinear()
    .domain([xMin - xTickSize, xMax + xTickSize])
    .range([
      Y_AXIS_WIDTH + PADDING * 2,
      Y_AXIS_WIDTH + chartWidth - PADDING - PADDING
    ]);
  const xTickWidth = (x.range()[1] - x.range()[0]) / NUM_TICKS;
  const bins = d3Array
    .bin()
    .value(d => d[binParam])
    .domain(x.domain())
    .thresholds(x.ticks(NUM_TICKS))(data);

  const maxY = Math.max(...bins.map(row => row.length));

  const y = d3
    .scaleLinear()
    .domain([0, maxY])
    .range([chartHeight, PADDING * 3]);

  const barScale = d3
    .scaleLinear()
    .domain([0, maxY])
    .range([0, chartHeight - PADDING * 3]);

  const ref = useCanvas(
    canvas => {
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
        xTickWidth,
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
        width: chartDim["width"]
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
            position: "absolute",
            zIndex: 100
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
        <Grid item style={{ position: "absolute" }}>
          <canvas ref={ref} id="probabilityCanvas" />
          <div id="probabilityTooltip" />
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
  xTickWidth,
  maxBin
) => {
  context.beginPath();
  context.globalAlpha = 1;
  context.fillStyle = "black";
  context.textAlign = "right";

  context.font = TICK_FONT;

  x.ticks(10).forEach(tick => {
    context.fillText(tick, x(tick), y(0) + 15);
  });

  const format = tick => (tick === 0 ? "0" : d3.format(".2f")(tick / maxBin));

  y.ticks(10).forEach(tick => {
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
  context.fillText(binParam, chartWidth / 2, chartHeight + 30);

  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText("Density", -(chartHeight / 2), 3);
  context.restore();
};

const drawBars = (context, bins, x, y, barScale) => {
  context.globalAlpha = 1;
  context.fillStyle = BAR_COLOR;
  context.strokeStyle = BAR_STROKE_COLOR;

  bins.forEach(bin => {
    const xPos = x(bin["x0"]) + 1;
    const yPos = y(bin.length);
    const width = x(bin["x1"]) - x(bin["x0"]) - 2;
    const height = barScale(bin.length);
    context.fillRect(xPos, yPos, width, height);
    context.strokeRect(xPos, yPos, width, height);
  });

  d3.select("#probabilityCanvas").on("mousemove", function() {
    var mouseX = d3.event.layerX || d3.event.offsetX;
    var mouseY = d3.event.layerY || d3.event.offsety;
    if (mouseX >= x.range()[0] && mouseX <= x.range()[1]) {
      const tickSize = (x.range()[1] - x.range()[0]) / NUM_TICKS;
      x.ticks(NUM_TICKS).map((tick, index) => {
        if (mouseX >= x(tick) && mouseX <= x(tick) + tickSize) {
          //show tip
          const bin = bins.filter(bin => bin["x0"] === tick)[0];
          if (bin.length > 0) {
            const yPos = y(bins[index].length);

            const subtypeGroups = _.groupBy(bin, "subtype");

            const tooltipWidth =
              Math.max(
                ...Object.keys(subtypeGroups).map(group => group.length)
              ) > 15
                ? 250
                : 150;

            d3.select("#probabilityTooltip")
              .attr("width", tooltipWidth + "px")
              .style("opacity", 0.8)
              .style(
                "left",
                x(tick) -
                  tooltipWidth / 2 +
                  tickSize / 2 +
                  PADDING / 2 -
                  (tooltipWidth === 250 ? 14 : 0) +
                  "px"
              )
              .style(
                "top",
                y(bin.length) -
                  Object.keys(subtypeGroups).length * 20 -
                  55 +
                  "px"
              )
              .html(function(d) {
                return (
                  "<p><ul style='width:" +
                  tooltipWidth +
                  "px'>" +
                  Object.keys(subtypeGroups)
                    .sort((a, b) => a.length - b.length)
                    .map(
                      group =>
                        "<li>" +
                        group +
                        " : " +
                        format(subtypeGroups[group].length / bin.length) +
                        "</li>"
                    )
                    .join("") +
                  "</ul></p>"
                );
              });
          } else {
            hideTooltip();
          }
        }
      });
    } else {
      hideTooltip();
    }
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
      datum => datum[barParam] === highlightedBar
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
    thresholds.map(t => [t, d3.mean(data, d => kernel(t - d))]);

  function epanechnikov(bandwidth) {
    return x =>
      Math.abs((x /= bandwidth)) <= 1 ? (0.75 * (1 - x * x)) / bandwidth : 0;
  }
  const densityData = highlightedLine
    ? data.filter(datum => datum[lineParam] === highlightedLine)
    : data;

  const density = kde(
    epanechnikov(1),
    x.ticks(NUM_TICKS * 2),
    densityData.map(row => parseFloat(row[binParam]))
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
const hideTooltip = () => d3.select("#probabilityTooltip").style("opacity", 0);
export default ProbabilityHistogram;
