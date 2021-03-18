import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis, changeFontSize } from "../DrawingUtils/utils.js";

const Histogram = ({ data, chartDim }) => {
  const [{ logXParam, logYParam, fontSize }] = useDashboardState();
  const [drawReady, setDrawReady] = useState(false);
  const [context, saveContext] = useState(null);

  const allX = data.map(row => parseFloat(row[logXParam]));
  const allY = data.map(row => parseFloat(row[logYParam]));
  const xMax = Math.max(...allX);
  const xMin = Math.min(...allX);

  const x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"] - 50]);

  const bins = d3Array
    .bin()
    .value(d => d[logXParam])
    .domain(x.domain())
    .thresholds(x.ticks(10));

  const buckets = bins(data);
  const maxY = Math.max(...buckets.map(row => row.length));

  // Y axis
  const y = d3
    .scaleLinear()
    .domain([0, maxY / data.length])
    .range([chartDim["chart"]["y2"], chartDim["chart"]["y1"]]);

  function drawBars(context, data) {
    context.fillStyle = "#6bb9f0";
    context.strokeStyle = "#5c97bf";
    buckets.map(bar => {
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
  }
  useEffect(() => {
    if (data.length > 0) {
      init(data, chartDim);
    }
  }, [data]);
  useEffect(() => {
    if (drawReady) {
      drawAxisLabels(context);
      drawBars(context, data);
    }
  }, [drawReady]);

  function drawAxisLabels(context) {
    const format = tick => (tick === 0 ? "0" : d3.format(".2f")(tick));
    context.beginPath();
    context.globalAlpha = 1;
    //context.lineWidth = 1;
    context.fillStyle = "black";
    context.textAlign = "right";
    const xMinValue = x(xMin) - 10;
    changeFontSize(context, fontSize["tickLabelFontSize"]);
    x.ticks(10).map(tick => {
      context.fillText(tick, x(tick), y(0) + 15);
    });
    changeFontSize(context, fontSize["axisLabelFontSize"]);
    context.fillText(
      logXParam,
      (chartDim["chart"]["x2"] - chartDim["chart"]["x1"]) / 2 +
        chartDim["chart"]["x1"],
      chartDim["chart"]["y2"] + 30
    );

    context.restore();
    y.ticks(10).map(tick => {
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
      chartDim["chart"]["x1"] - 35
    );
    context.restore();
  }
  function init(data, chartDim) {
    var canvas = d3.select("#histogramCanvas");

    var currContext = canvasInit(canvas, chartDim.width, chartDim.height);

    currContext.fillStyle = "white";
    currContext.fillRect(0, 0, chartDim.width, chartDim.height);
    saveContext(currContext);
    setDrawReady(true);
  }

  function drawLegend() {}
  return (
    <div class="card" style={{ margin: 10 }}>
      <div
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative"
        }}
      >
        <div
          id="histogram"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="histogramCanvas" />
        </div>
      </div>
    </div>
  );
};
export default Histogram;
