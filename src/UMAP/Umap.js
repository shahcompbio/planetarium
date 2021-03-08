import React, { useEffect, useState, useRef, useCallback } from "react";
import * as d3 from "d3";
import * as d3Collection from "d3-collection";
import * as d3Array from "d3-array";
import _ from "lodash";

import { useDashboardState } from "../PlotState/dashboardState";

import {
  canvasInit,
  drawAxisLabels,
  drawAxisTicks,
  drawAxis
} from "../DrawingUtils/utils.js";

const Umap = ({ data, chartDim }) => {
  const [
    {
      xParam,
      yParam,
      cellIdParam,
      sampleType,
      clonotypeParam,
      sampleTen,
      topTen,
      colors,
      clonotypes
    }
  ] = useDashboardState();

  useEffect(() => {
    drawAll(data, "NDVL", chartDim);
  }, [data]);

  function drawOutline(context, x, y, data, colors) {
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    context.globalAlpha = 0.5;
    data.forEach(point => {
      context.beginPath();
      context.arc(
        x(point[xParam]),
        y(point[yParam]),
        1.5,
        0,
        Math.PI * 2,
        true
      );
      context.fillStyle = "#d2d7d3";

      context.fill();
    });
    context.globalAlpha = 1;
  }

  export function drawPoint(context, point, fill, isCountInsignificant, x, y) {
    const radius = isCountInsignificant ? 2 : point["radius"];

    context.beginPath();
    context.arc(
      x(point[xParam]),
      y(point[yParam]),
      radius,
      0,
      Math.PI * 2,
      true
    );

    context.fillStyle = fill;
    context.fill();
  }
  function drawLineGraph(
    context,
    x,
    y,
    dimensions,
    data,
    topTen,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    const maxValue = Math.max(...Object.entries(topTen).map(row => row[1]));

    const lineXFreq = d3
      .scaleLinear()
      .domain([0, maxValue / 2.5])
      .range([dimensions.x2, dimensions.x2 + 25]);

    const lineYFreq = d3
      .scaleLinear()
      .domain([0, maxValue / 2.5])
      .range([dimensions.y1, dimensions.y1 - 25]);

    const lineXaxis = d3
      .line()
      .x(function(d) {
        return x((d.x0 + d.x1) / 2);
      })
      .y(function(d) {
        const freq = Object.entries(d).length;
        return lineYFreq(freq);
      })
      .curve(d3.curveCatmullRom.alpha(0.5))
      .context(context);

    const lineYaxis = d3
      .line()
      .y(function(d) {
        return y((d.x0 + d.x1) / 2);
      })
      .x(function(d) {
        const freq = Object.entries(d).length;
        return lineXFreq(freq);
      })
      .curve(d3.curveCatmullRom.alpha(0.5))
      .context(context);

    var nestedSamples = Array.from(
      d3Array.group(data, d => d[clonotypeParam]),
      ([key, value]) => ({ key, value })
    );

    //if selected, move to end so it's drawn last
    if (selectedClonotype) {
      const keys = nestedSamples.map(row => row["key"]);

      nestedSamples = [
        ...nestedSamples.filter(row => row["key"] !== selectedClonotype),
        nestedSamples[keys.indexOf(selectedClonotype)]
      ];
    }

    nestedSamples.reduce((final, clonotype) => {
      const xBins = d3Array
        .bin()
        .value(d => d[xParam])
        .domain(x.domain())
        .thresholds(x.ticks(10))(clonotype["value"]);

      context.beginPath();
      lineXaxis(xBins);

      const isSelected =
        selectedClonotype && selectedClonotype === clonotype["key"];

      context.lineWidth = isSelected ? 2.5 : 1.5;

      context.strokeStyle = selectedClonotype
        ? isSelected
          ? colors(clonotype["key"])
          : "#e8e8e8"
        : colors(clonotype["key"]);

      context.stroke();

      const yBins = d3Array
        .bin()
        .value(d => d[yParam])
        .domain(y.domain())
        .thresholds(y.ticks(10))(clonotype["value"]);

      context.beginPath();
      lineYaxis(yBins);

      context.lineWidth = isSelected ? 2.5 : 1.5;

      context.strokeStyle = selectedClonotype
        ? isSelected
          ? colors(clonotype["key"])
          : "#e8e8e8"
        : colors(clonotype["key"]);

      context.stroke();
      final[clonotype["key"]] = { x: xBins, y: yBins };
      return final;
    }, []);
  }

  function drawPoints(
    context,
    x,
    y,
    dimensions,
    data,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    var nestedSamples = Array.from(
      d3Array.group(data, d => d[clonotypeParam]),
      ([key, value]) => ({ key, value })
    );

    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";

    const merge = nestedSamples.map((clonotype, i) => {
      const type = clonotype["key"];
      const values = clonotype["value"];

      const xBins = d3Array
        .bin()
        .value(d => d[xParam])
        .domain(x.domain())
        .thresholds(x.ticks(8))(clonotype["value"]);

      const yBins = d3Array
        .bin()
        .value(d => d[yParam])
        .domain(y.domain())
        .thresholds(y.ticks(8))(clonotype["value"]);

      const xBinned = xBins.reduce((final, xBin) => {
        const freq = Object.entries(xBin).length - 2;

        if (freq > 0) {
          const rows = Object.entries(xBin).filter((row, index) => {
            if (row[0] !== "x1" && row[1] !== "x0") {
              return true;
            } else {
              return false;
            }
          });

          final = {
            ...final,
            ...rows.reduce((finalRow, row) => {
              finalRow[row[1][cellIdParam]] = { ...row[1], xRadius: freq };
              return finalRow;
            }, {})
          };
        }
        return final;
      }, {});

      const yBinned = yBins.reduce((final, yBin) => {
        const freq = Object.entries(yBin).length - 2;

        if (freq > 0) {
          const rows = Object.entries(yBin).filter((row, index) => {
            if (row[0] !== "x1" && row[1] !== "x0") {
              return true;
            } else {
              return false;
            }
          });

          final = {
            ...final,
            ...rows.reduce((finalRow, row) => {
              finalRow[row[1][cellIdParam]] = { ...row[1], yRadius: freq };
              return finalRow;
            }, {})
          };
        }
        return final;
      }, {});

      const merged = Object.keys(yBinned).map((item, i) =>
        Object.assign({}, yBinned[item], xBinned[item])
      );
      return merged;
    });

    const sortedmerge = merge
      .flat(1)
      .filter(point => point.hasOwnProperty(cellIdParam))
      .map(point => ({
        ...point,
        radius: (point["xRadius"] + point["yRadius"]) / 3
      }))
      .sort((a, b) => b.radius - a.radius);

    const isCountInsignificant =
      Math.max(...sortedmerge.map(point => point["radius"])) < 1 ? true : false;

    sortedmerge.map(point => {
      const fill = selectedClonotype ? "grey" : colors(point[clonotypeParam]);
      context.globalAlpha = selectedClonotype ? 0.5 : 1;
      drawPoint(context, point, fill, isCountInsignificant, x, y);
    });

    context.globalAlpha = 1;
    //if selected, move to end so it's drawn last
    if (selectedClonotype) {
      sortedmerge
        .filter(row => row[clonotypeParam] === selectedClonotype)
        .map(point => {
          const fill = colors(point[clonotypeParam]);
          drawPoint(context, point, fill, isCountInsignificant, x, y);
        });
    }
  }
  export const clearAll = (context, chartDim) =>
    context.clearRect(
      0,
      0,
      chartDim["chart"].x2 + 30,
      chartDim["chart"].y2 + 30
    );

  function reDraw(
    context,
    x,
    y,
    dim,
    sampleData,
    sampleTen,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    drawAxis(context, x, y, dim);
    //  drawAxisLabels(context, x, y, chartDim);
    drawOutline(context, x, y, data, colors);
    drawLineGraph(
      context,
      x,
      y,
      dim,
      sampleData,
      sampleTen,
      colors,
      topTenNumbering,
      selectedClonotype
    );
    drawPoints(
      context,
      x,
      y,
      dim,
      sampleData,
      colors,
      topTenNumbering,
      selectedClonotype
    );
  }
  function drawAll(data, sampleType, chartDim) {
    var canvas = d3.select("#canvas");

    var context = canvasInit(canvas, chartDim.width, chartDim.height);

    context.fillStyle = "white";
    context.fillRect(0, 0, chartDim.width, chartDim.height);

    const yData = data.map(d => parseFloat(d[yParam]));
    const xData = data.map(d => parseFloat(d[xParam]));

    const yMin = Math.min(...yData);
    const yMax = Math.max(...yData);
    const xMin = Math.min(...xData);
    const xMax = Math.max(...xData);

    const sampleTen = topTen.reduce((final, curr) => {
      final[curr[0]] = curr[1];
      return final;
    }, {});

    const sampleData = data.filter(row =>
      sampleTen.hasOwnProperty(row[clonotypeParam])
    );

    const dim = chartDim["chart"];
    // X axis
    var x = d3
      .scaleLinear()
      .domain([xMin, xMax])
      .range([dim.x1, dim.x2]);

    // Y axis
    var y = d3
      .scaleLinear()
      .domain([yMin, yMax])
      .range([dim.y2, dim.y1]);
    const topTenNumbering = topTen.reduce((final, seq) => {
      const label =
        sampleType === "NDVL"
          ? "L" + (Object.keys(final).length + 1)
          : "R" + (Object.keys(final).length + 1);
      final[seq[0]] = label;
      return final;
    }, {});
    reDraw(context, x, y, dim, sampleData, sampleTen, colors, topTenNumbering);

    var legend = d3.select("#legend");
    legend = legend.append("g");
    legend
      .selectAll("circle")
      .data(topTen)
      .enter()
      .append("circle")
      .attr("r", 6)
      .attr("cx", function(d) {
        return chartDim["legend"].x1 + 5;
      })
      .attr("cy", function(d, i) {
        return i * 20 + chartDim["legend"].y1 + 30;
      })
      .attr("fill", function(d) {
        return colors(d[0]);
      })
      .on("mouseover", function(d, i) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(
          context,
          x,
          y,
          dim,
          sampleData,
          sampleTen,
          colors,
          topTenNumbering,
          d[0]
        );
      })
      .on("mouseout", function(event, d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(
          context,
          x,
          y,
          dim,
          sampleData,
          sampleTen,
          colors,
          topTenNumbering
        );
      });
    legend
      .selectAll("text")
      .data(topTen)
      .enter()
      .append("text")
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 20;
      })
      .attr("y", function(d, i) {
        return i * 20 + chartDim["legend"].y1 + 31;
      })
      .attr("dy", ".35em")
      .text(function(d) {
        return topTenNumbering[d[0]] + " - " + d[0] + " - " + d[1];
        //  return d["key"];
      })
      .attr("font-weight", "700")
      .attr("fill", function(d) {
        return colors(d[0]);
      })
      .attr("cursor", "pointer")
      .on("mouseover", function(d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(
          context,
          x,
          y,
          dim,
          sampleData,
          sampleTen,
          colors,
          topTenNumbering,
          d[0]
        );
      })
      .on("mouseout", function(event, d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(
          context,
          x,
          y,
          dim,
          sampleData,
          sampleTen,
          colors,
          topTenNumbering
        );
      });

    legend
      .append("g")
      .append("text")
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 17;
      })
      .attr("y", function(d, i) {
        return chartDim["legend"].y1 + 15;
      })
      .attr("font-size", 20)
      .text(sampleType);
  }

  return (
    <div>
      <div style={{ width: 600, height: 700, position: "relative" }}>
        <div
          id="scatterplot"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="canvas" />
          <svg id="legend" style={{ float: "right", width: 600 }} />
        </div>
      </div>
    </div>
  );
};
export default Umap;
