import React, { useEffect } from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useDashboardState } from "../PlotState/dashboardState";
import { drawPoint, clearAll } from "./Umap.js";
import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";

const SubtypeUmap = ({
  data,
  chartDim,
  selectedSubtype,
  selectedClonotype,
  setSelectedSubtype,
  setSelectedClonotype
}) => {
  const [{ xParam, yParam, subtypeParam }] = useDashboardState();
  useEffect(() => {
    if (data.length > 0) {
      drawAll(data, chartDim);
    }
  }, [data]);

  function drawAll(data, chartDim) {
    var canvas = d3.select("#subTypeUmapCanvas");

    var context = canvasInit(canvas, chartDim.width, chartDim.height);

    context.fillStyle = "white";
    context.fillRect(0, 0, chartDim.width, chartDim.height);

    const yData = data.map(d => parseFloat(d[yParam]));
    const xData = data.map(d => parseFloat(d[xParam]));

    const yMin = Math.min(...yData);
    const yMax = Math.max(...yData);
    const xMin = Math.min(...xData);
    const xMax = Math.max(...xData);

    // X axis
    var x = d3
      .scaleLinear()
      .domain([xMin, xMax])
      .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"]]);

    // Y axis
    var y = d3
      .scaleLinear()
      .domain([yMax, yMin])
      .range([chartDim["chart"]["y1"], chartDim["chart"]["y2"]]);
    drawAxis(context, x, y, chartDim["chart"]);
    drawPoints(data, chartDim, context, x, y);
  }
  function drawPoints(data, chartDim, context, x, y, selectedSubtype) {
    const subTypes = _.groupBy(data, subtypeParam);

    const types = Object.keys(subTypes).map(type => type.replace(/\s/g, ""));
    var colors = d3
      .scaleOrdinal()
      .domain([...types])
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

    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";

    data.forEach(point => {
      const type = point[subtypeParam].replace(/\s/g, "");
      const fill =
        selectedSubtype && selectedSubtype !== type
          ? "#e8e8e8"
          : colors(point[subtypeParam]);
      drawPoint(context, point, fill, true, x, y, xParam, yParam);
    });

    appendSubtypeLabels(subTypes, x, y, yParam, xParam, context, colors);
    appendLegend(colors, types, context, x, y);
  }
  function appendSubtypeLabels(
    subTypes,
    x,
    y,
    yParam,
    xParam,
    context,
    colors
  ) {
    var boundingBox = {};
    const allBoxes = Object.keys(subTypes).map(subtype => {
      const type = subtype.replace(/\s/g, "");
      const ySubtypeData = subTypes[subtype].map(d => x(parseFloat(d[yParam])));
      const xSubtypeData = subTypes[subtype].map(d => y(parseFloat(d[xParam])));

      const xDataNoOutliers = filterOutliers(xSubtypeData, [50, 75]);
      const yDataNoOutliers = filterOutliers(ySubtypeData, [50, 75]);

      const intersection = _.intersectionBy(
        xDataNoOutliers,
        yDataNoOutliers,
        "index"
      ).reduce((final, curr) => {
        final["i-" + curr["index"]] = curr["value"];
        return final;
      }, {});

      const remainingPoints = subTypes[subtype].filter(
        (d, i) => intersection["i-" + i]
      );

      return {
        type: type,
        subtype: subtype,
        box: getBoundingBox(remainingPoints)
      };
    });
    const avgBoxSize =
      allBoxes
        .map(box => x(box["box"]["right"]) - x(box["box"]["left"]))
        .reduce((a, b) => a + b) / allBoxes.length;

    //if the bounding box is larger than 1.5 the  avg do not add label
    allBoxes.map(box => {
      //    if (
      //      !(x(box["box"]["right"]) - x(box["box"]["left"]) * 1.5 > avgBoxSize)
      //    ) {
      drawBoundingBox(context, box, x, y, colors);
      //    }
    });
  }
  function reDraw(data, chartDim, context, x, y, selectedSubtype) {
    drawAxis(context, x, y, chartDim["chart"]);
    drawPoints(data, chartDim, context, x, y, selectedSubtype);
  }
  function appendLegend(colors, subTypes, context, x, y) {
    var legend = d3.select("#subTypeUmapLegend");
    legend = legend.append("g");
    legend
      .selectAll("circle")
      .data(subTypes)
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
      .attr("cursor", "pointer")
      .on("mouseover", function(d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y, d[0]);
      })
      .on("mouseout", function(event, d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y);
      });

    legend
      .selectAll("text")
      .data(subTypes)
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
        return d;
      })
      .attr("font-weight", "700")
      .attr("fill", function(d) {
        return colors(d[0]);
      })
      .attr("cursor", "pointer")
      .on("mouseover", function(d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y, d[0]);
      })
      .on("mouseout", function(event, d) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y);
      });
  }
  function filterOutliers(coords, quantiles) {
    const sortedCollection = coords.slice().sort((a, b) => a - b); //copy array fast and sort
    const size = sortedCollection.length;

    let q1 = getQuantile(sortedCollection, quantiles[0]);
    let q3 = getQuantile(sortedCollection, quantiles[1]);

    const iqr = q3 - q1;
    const maxValue = q3 + iqr * 1.5;
    const minValue = q1 - iqr * 1.5;

    return coords
      .map((value, index) => ({ value: value, index: index }))
      .filter(d => d.value >= minValue && d.value <= maxValue);
  }
  function getQuantile(array, quantile) {
    let index = (quantile / 100.0) * (array.length - 1);

    if (index % 1 === 0) {
      return array[index];
    } else {
      let lowerIndex = Math.floor(index);
      let remainder = index - lowerIndex;
      return (
        array[lowerIndex] +
        remainder * (array[lowerIndex + 1] - array[lowerIndex])
      );
    }
  }
  function drawBoundingBox(context, box, x, y, colors) {
    const type = box["type"];
    const title = box["subtype"];
    const boxCords = box["box"];
    context.restore();
    context.fillStyle = colors(type);
    context.strokeStyle = "black";
    context.lineWidth = 1;
    context.beginPath();

    const width =
      boxCords.right > boxCords.left
        ? x(boxCords.right) - x(boxCords.left)
        : x(boxCords.left) - x(boxCords.right);
    const height = x(boxCords.bottom) - x(boxCords.top);
    const subtypeFontSize = title.length === 1 ? 18 : 13;
    context.font = "normal bold " + subtypeFontSize + "px Droid";

    const textWidth = context.measureText(title).width + 2;

    const textHeight = context.measureText(title);

    context.fillStyle = "white";
    context.globalAlpha = 1;

    if (title.indexOf("/") !== -1) {
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10,
        textWidth / 2,
        subtypeFontSize
      );
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10 + subtypeFontSize,
        textWidth / 2,
        subtypeFontSize
      );
    } else {
      context.fillRect(
        x(boxCords.left) + width / 2,
        y(boxCords.top) - height / 2 - 14,
        textWidth,
        subtypeFontSize
      );
      /*context.fillStyle = "black";
      context.globalAlpha = 0.5;
      context.fillRect(
        x(boxCords.left),
        y(boxCords.top) - height,
        width,
        height
      );*/
    }
    context.fill();
    context.save();

    context.globalAlpha = 1;
    context.fillStyle = "black";

    //  context.font = "normal bold " + subtypeFontSize + "px Droid";
    //iff titlee is too long, 2 lines
    if (title.indexOf("/") !== -1) {
      const splitTitle = title.split("/");

      context.fillText(
        splitTitle[0] + "/",
        x(boxCords.left),
        y(boxCords.top) + height / 2
      );
      context.fillText(
        splitTitle[1].replace(/\s/g, ""),
        x(boxCords.left),
        y(boxCords.top) + height / 2 + subtypeFontSize
      );
    } else {
      context.fillText(
        title,
        x(boxCords.left) + width / 2,
        y(boxCords.top) - height / 2
      );
    }

    context.fill();
    context.stroke();
    context.save();
  }
  function getBoundingBox(data) {
    return data.reduce((final, point) => {
      const xPoint = parseFloat(point[xParam]);
      const yPoint = parseFloat(point[yParam]);
      if (final["top"]) {
        if (yPoint < final["top"]) {
          final["top"] = yPoint;
        }
        if (xPoint < final["left"]) {
          final["left"] = xPoint;
        }
        if (yPoint > final["bottom"]) {
          final["bottom"] = yPoint;
        }
        if (xPoint > final["right"]) {
          final["right"] = xPoint;
        }
      } else {
        final["top"] = yPoint;
        final["left"] = xPoint;
        final["bottom"] = yPoint;
        final["right"] = xPoint;
      }
      return final;
    }, {});
  }
  return (
    <div>
      <div style={{ width: 600, height: 700, position: "relative" }}>
        <div
          id="subTypeUmap"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="subTypeUmapCanvas" />
          <svg id="subTypeUmapLegend" style={{ float: "right", width: 600 }} />
        </div>
      </div>
    </div>
  );
};

export default SubtypeUmap;
