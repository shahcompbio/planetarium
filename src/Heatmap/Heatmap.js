import React, { useState, useEffect } from "react";
import * as d3 from "d3";

import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit } from "../DrawingUtils/utils.js";
const clearAll = (context, chartDim) =>
  context.clearRect(0, 0, chartDim["chart"].x2, chartDim["chart"].y2);
const Heatmap = ({
  data,
  chartDim,
  selectedSubtype,
  selectedClonotype,
  setSelectedSubtype,
  setSelectedClonotype
}) => {
  const [
    { clonotypeParam, sampleTen, topTenNumbering, colors, subtypeParam }
  ] = useDashboardState();
  const [context, saveContext] = useState(null);

  const subTypes = data.reduce((final, current) => {
    var allSamples = final;
    const subtype = current[subtypeParam];
    allSamples[subtype] = final.hasOwnProperty(subtype)
      ? allSamples[subtype] + 1
      : 1;
    final = allSamples;
    return final;
  }, {});

  useEffect(() => {
    if (context) {
      drawHeatmap(context, chartDim, data, subTypes);
    }
  }, [context]);

  useEffect(() => {
    if (data.length > 0 && colors) {
      init(data, chartDim);
    }
  }, [data, colors]);

  useEffect(() => {
    if (context) {
      if (selectedSubtype !== null) {
        clearAll(context, chartDim);
        drawHeatmap(context, chartDim, data, selectedSubtype);
      } else {
        clearAll(context, chartDim);
        drawHeatmap(context, chartDim, data);
      }
    }
  }, [selectedSubtype, context]);

  function init(data, chartDim) {
    var canvas = d3.select("#heatmapCanvas");

    var context = canvasInit(canvas, chartDim.width, chartDim.height);

    context.fillStyle = "white";
    context.fillRect(0, 0, chartDim.width, chartDim.height);

    saveContext(context);
  }

  function drawHeatmap(context, allDim, data, selectedSubtype) {
    const dimensions = allDim;
    const allSubtypes = Object.keys(subTypes);

    var largestFreq = 0;
    const subtypeStats = data.reduce(
      (final, current) => {
        const seq = current[clonotypeParam];
        if (final[seq]) {
          const subtype = current[subtypeParam];

          if (final[seq].hasOwnProperty(subtype)) {
            final[seq][subtype] = final[seq][subtype] + 1;
          } else {
            final[seq][subtype] = 1;
          }
          largestFreq =
            final[seq][subtype] > largestFreq
              ? final[seq][subtype]
              : largestFreq;
        }
        return final;
      },
      {
        ...Object.keys(sampleTen)
          .sort(([, a], [, b]) => b - a)
          .reduce((final, entry) => {
            final[entry] = {};
            return final;
          }, {})
      }
    );
    const xAxis = d3
      .scaleBand()
      .domain(allSubtypes)
      .range([dimensions["chart"]["x1"], dimensions["chart"]["x2"] - 50]);

    const freqColouring = d3
      .scaleLinear()
      .range(["#ffec8b", "#d91e18"])
      .domain([0, 96]);

    const heatmapWidth =
      (dimensions["chart"]["x2"] - dimensions["chart"]["x1"] - 50) /
      allSubtypes.length;

    //  if (currSample === "NDVL") {
    //add labels once
    context.beginPath();
    context.fillStyle = "#000000";

    allSubtypes.forEach(function(d, index) {
      context.save();
      context.translate(
        xAxis(d) + heatmapWidth / 2,
        dimensions["chart"]["y1"] - 5
      );
      context.rotate((322 * Math.PI) / 180);
      context.fillText(d, 0, 0);
      context.restore();
    });
    //  }
    const sequenceLength = Object.keys(subtypeStats).length;
    const heatmapRowSpace = 3;
    const heatmapHeight =
      (dimensions["chart"]["y2"] -
        dimensions["chart"]["y1"] -
        sequenceLength * heatmapRowSpace) /
      sequenceLength;
    const startingY = dimensions["chart"]["y1"];

    const alphaIndexing = Object.entries(topTenNumbering)
      .map(entry => entry[1])
      .sort();

    Object.keys(subtypeStats)
      .sort((a, b) => {
        return (
          alphaIndexing.indexOf(topTenNumbering[a]) -
          alphaIndexing.indexOf(topTenNumbering[b])
        );
      })
      .map((sequence, index) => {
        const yPos =
          startingY + heatmapHeight * index + heatmapRowSpace * index;
        context.fillStyle = colors(sequence);

        context.fillText(
          topTenNumbering[sequence] + " - " + sequence,
          dimensions["chart"]["x2"] - 50,
          yPos + (3 * heatmapHeight) / 4
        );
        const seqSubtypes = subtypeStats[sequence];

        allSubtypes.map(subtype => {
          context.globalAlpha =
            selectedSubtype === null
              ? selectedSubtype === subtype
                ? 1
                : 0.2
              : 1;
          context.font = "20px";
          if (seqSubtypes.hasOwnProperty(subtype)) {
            context.fillStyle = freqColouring(subtypeStats[sequence][subtype]);
          } else {
            context.fillStyle = "#eeeeee";
          }
          context.fillRect(
            xAxis(subtype),
            yPos,
            heatmapWidth - 3,
            heatmapHeight
          );
          if (seqSubtypes.hasOwnProperty(subtype)) {
            context.fillStyle = "black";
            const freq = subtypeStats[sequence][subtype];
            const freqTextX =
              freq > 9
                ? xAxis(subtype) + (heatmapWidth - 15) / 2
                : xAxis(subtype) + (heatmapWidth - 6) / 2;
            context.fillText(freq, freqTextX, yPos + heatmapHeight / 2 + 5);
          }
        });
      });
  }
  return (
    <div>
      <div style={{ width: 600, height: 700, position: "relative" }}>
        <div
          id="heatmap"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="heatmapCanvas" />
          <svg id="legend" style={{ float: "right", width: 600 }} />
        </div>
      </div>
    </div>
  );
};
export default Heatmap;
