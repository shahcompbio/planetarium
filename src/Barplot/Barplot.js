import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";

const Barplot = ({ data, chartDim }) => {
  const [
    { xParam, yParam, cellIdParam, clonotypeParam, topTen, subtypeParam }
  ] = useDashboardState();
  const [drawReady, setDrawReady] = useState(false);
  const [context, saveContext] = useState(null);

  const groupedData = _.groupBy(data, subtypeParam);
  const subtypes = Object.keys(groupedData);
  const stackedBarData = subtypes.reduce((final, subtype) => {
    const groupedClonotypes = _.groupBy(groupedData[subtype], clonotypeParam);
    var counter = {};
    const hitList = Object.entries(groupedClonotypes).map(
      cellHits => cellHits[1].length
    );
    hitList.forEach(x => (counter[x] = (counter[x] || 0) + 1));
    final[subtype] = {
      total: hitList.length,
      counts: counter
    };
    return final;
  }, []);

  // X axis
  const x = d3
    .scaleBand()
    .domain(subtypes)
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"]]);
  // Y axis
  const y = d3
    .scaleLinear()
    .domain([100, 0])
    .range([chartDim["chart"]["y1"], chartDim["chart"]["y2"]]);
  var colors = d3
    .scaleOrdinal()
    .domain([...Array.from(Array(10).keys())])
    .range([
      "#5E4FA2",
      "#3288BD",
      "#66C2A5",
      "#ABDDA4",
      "#E6F598",
      "#FFFFBF",
      "#FEE08B",
      "#FDAE61",
      "#F46D43",
      "#D53E4F",
      "#9E0142"
    ]);
  function drawBars(context) {
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    subtypes.map(subtype => {
      var currentHeight = 0;
      [...Array.from(Array(10).keys())].map((key, index) => {
        const { counts, total } = stackedBarData[subtype];
        var height;
        if (index === 9) {
          const allOther = Object.entries(counts)
            .filter(row => row[0] > 9)
            .map(row => row[1]);
          height =
            allOther.length > 0
              ? (allOther.reduce((a, b) => a + b) / total) * 100
              : 0;
        } else {
          height = counts[key + 1] ? (counts[key + 1] / total) * 100 : 0;
        }
        context.fillStyle = colors(key);
        context.fillRect(
          x(subtype),
          y(height + currentHeight),
          50,
          y(0) - y(height)
        );
        currentHeight += height;
        context.fill();
      });
    });
  }
  useEffect(() => {
    if (data.length > 0) {
      init(data, chartDim);
    }
  }, [data]);
  useEffect(() => {
    if (drawReady) {
      //  drawAxis(context, x, y, chartDim["chart"]);
      drawBars(context);
      drawAxisLabels(context);
      drawLegend();
    }
  }, [drawReady]);
  function drawAxisLabels(context) {
    context.beginPath();

    context.globalAlpha = 1;
    context.lineWidth = 1;
    context.fillStyle = "black";
    subtypes.map(subtype => {
      context.save();
      context.translate(x(subtype), y(0));
      context.rotate((322 * Math.PI) / 180);
      context.fillText(subtype, 0, 0);
      context.stroke();
      context.restore();
    });
  }
  function init(data, chartDim) {
    var canvas = d3.select("#barplotCanvas");

    var currContext = canvasInit(canvas, chartDim.width, chartDim.height);

    currContext.fillStyle = "white";
    currContext.fillRect(0, 0, chartDim.width, chartDim.height);
    saveContext(currContext);
    setDrawReady(true);
  }

  function drawLegend() {}
  return (
    <div>
      <div style={{ width: 600, height: 700, position: "relative" }}>
        <div
          id="barchart"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="barplotCanvas" />
        </div>
      </div>
    </div>
  );
};
export default Barplot;
