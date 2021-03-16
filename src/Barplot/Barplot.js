import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, changeFontSize } from "../DrawingUtils/utils.js";

const Barplot = ({ data, chartDim }) => {
  const [{ clonotypeParam, subtypeParam, fontSize }] = useDashboardState();

  const [drawReady, setDrawReady] = useState(false);
  const [context, saveContext] = useState(null);
  const barWidth = 50;
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
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"] - 30]);
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
    subtypes.forEach(subtype => {
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
          barWidth,
          y(0) - y(height)
        );
        currentHeight += height;
        context.fill();
      });
    });
  }
  function drawLegend(context) {
    changeFontSize(context, fontSize.legendFontSize);
    [...Array.from(Array(10).keys())]
      .sort((a, b) => b - a)
      .map((key, index) => {
        context.fillStyle = colors(key);
        context.fillRect(
          chartDim["chart"]["x2"] - 40,
          chartDim["chart"]["y1"] + index * 14 + index * 2,
          9,
          9
        );
        context.fillStyle = "#000000";
        const legendText = key + 1 >= 10 ? "≥10" : key + 1;
        context.fillText(
          legendText,
          chartDim["chart"]["x2"] - 25,
          chartDim["chart"]["y1"] + index * 14 + index * 2 + 8
        );
        context.fill();
      });

    context.fillRect(
      chartDim["chart"]["x2"] + 20,
      chartDim["chart"]["y1"],
      5,
      5
    );
  }
  useEffect(() => {
    if (data.length > 0) {
      init(data, chartDim);
    }
  }, [data]);
  useEffect(() => {
    if (drawReady) {
      //  drawAxis(context, x, y, chartDim["chart"]);
      drawLegend(context);
      drawBars(context);
      drawAxisLabels(context);
      drawYAxisLabels(context);
    }
  }, [drawReady]);

  function drawYAxisLabels(context) {
    context.fillStyle = "black";
    context.textAlign = "right";
    context.lineWidth = 1;
    context.textBaseline = "middle";
    const ticks = y.ticks(10);

    context.beginPath();
    ticks.forEach(function(d) {
      changeFontSize(context, fontSize["tickLabelFontSize"]);
      //  context.moveTo(chartDim["margin"]["left"] + 17, y(d));
      //  context.lineTo(chartDim["margin"]["left"] + 27, y(d));
      context.fillText(d, chartDim["margin"]["left"] + 15, y(d));
      context.stroke();
    });
  }
  function drawAxisLabels(context) {
    context.beginPath();
    context.globalAlpha = 1;
    context.lineWidth = 1;
    context.fillStyle = "black";
    context.textAlign = "right";

    changeFontSize(context, fontSize["axisLabelFontSize"]);
    subtypes.map(subtype => {
      context.save();
      context.translate(x(subtype) + barWidth / 2, y(0) + 7);
      context.rotate((322 * Math.PI) / 180);
      if (subtype.indexOf("/") !== -1) {
        context.fillText(subtype.split("/")[0] + "/", -5, 0);
        context.fillText(subtype.split("/")[1], -5, 10);
      } else {
        context.fillText(subtype, -5, 5);
      }
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

  return (
    <div>
      <div
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative"
        }}
      >
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
