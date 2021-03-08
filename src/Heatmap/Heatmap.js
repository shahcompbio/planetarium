import React, { useEffect, useState, useRef, useCallback } from "react";
import * as d3 from "d3";
import * as d3Collection from "d3-collection";
import * as d3Array from "d3-array";
import * as d3Tsv from "d3-dsv";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import {
  canvasInit,
  drawAxisLabels,
  drawAxisTicks,
  drawAxis
} from "../DrawingUtils/utils.js";
const removeSubtypes = {
  C13: true,
  CD4: true,
  NK: true,
  "Cycling NK": true,
  "Cycling CD4": true,
  Treg: true,
  C15: true,
  C12: true
};
const Heatmap = ({ data, chartDim }) => {
  const [
    {
      xParam,
      yParam,
      cellIdParam,
      sampleType,
      clonotypeParam,
      sampleTen,
      topTenNumbering,
      colors,
      clonotypes
    }
  ] = useDashboardState();
  const subtypeParam = "seurat_clusters";

  useEffect(() => {
    drawAll(data, chartDim);
  }, [data]);

  function drawAll(data, chartDim) {
    var canvas = d3.select("#heatmapCanvas");

    var context = canvasInit(canvas, chartDim.width, chartDim.height);

    context.fillStyle = "white";
    context.fillRect(0, 0, chartDim.width, chartDim.height);

    const subTypes = data.reduce((final, current) => {
      var allSamples = final;
      const subtype = current[subtypeParam];
      allSamples[subtype] = final.hasOwnProperty(subtype)
        ? allSamples[subtype] + 1
        : 1;
      final = allSamples;
      return final;
    }, {});
    /*  const clonotypes = _.groupBy(data, clonotypeParam);
    const types = Object.keys(clonotypes);
    const colourList = [
      "#674172",
      "#098dde",
      "#fa832f",
      "#0e5702",
      "#c20c1e",
      "#911eb4",
      "#fc97bc",
      "#469990",
      "#b5762a",
      "#5aebed",
      "#8f8f3f",
      "#ed1a1a"
    ];
    var colors = d3
      .scaleOrdinal()
      .domain([...types])
      .range([...colourList]);*/

    drawHeatmap(context, chartDim, data, subTypes);
  }

  function drawHeatmap(context, allDim, data, subTypes) {
    const currSample = "NDVL";
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

    const yAxis = d3
      .scaleLinear()
      .domain([...Object.keys(sampleTen)])
      .range([dimensions["chart"]["y1"], dimensions["chart"]["y2"]]);

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
    const sortedKeys = Object.entries(sampleTen)
      .sort((a, b) => a - b)
      .map(entry => entry[0]);

    const alphaIndexing = Object.entries(topTenNumbering)
      .map(entry => entry[1])
      .sort();

    Object.keys(subtypeStats)
      .sort((a, b) => {
        //  console.log(topTenNumbering[a]);
        return (
          alphaIndexing.indexOf(topTenNumbering[a]) -
          alphaIndexing.indexOf(topTenNumbering[b])
        );
      })
      .map((sequence, index) => {
        const yPos =
          startingY + heatmapHeight * index + heatmapRowSpace * index;
        context.fillStyle = colors(sequence);
        //subtypeColours(sequence);
        //context.fillStyle = "#000000";
        context.fillText(
          topTenNumbering[sequence] + " - " + sequence,
          dimensions["chart"]["x2"] - 50,
          yPos + (3 * heatmapHeight) / 4
        );
        const seqSubtypes = subtypeStats[sequence];

        allSubtypes.map(subtype => {
          //  context.globalAlpha = 0.8;
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
          //context.globalAlpha = 1;
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
