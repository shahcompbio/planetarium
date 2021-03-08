import React, { useEffect } from "react";
import * as d3 from "d3";

import { useDashboardState } from "../PlotState/dashboardState";
//import { drawPoint, clearAll } from "./Umap.js";
import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";
const subtypeFontSize = 11;
const SubtypeUmap = ({ data, chartDim }) => {
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
    drawAll(data, chartDim);
  }, [data]);

  function drawAll(data, chartDim) {
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

    // X axis
    var x = d3
      .scaleLinear()
      .domain([xMin, xMax])
      .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"]]);

    // Y axis
    var y = d3
      .scaleLinear()
      .domain([yMin, yMax])
      .range([chartDim["chart"]["y1"], chartDim["chart"]["y2"]]);
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
    context.font = "normal bold 13px Droid";
    //    context.fillText(title, x(boxCords.left), y(boxCords.top) + height / 2);

    const textWidth = context.measureText(title).width + 2;

    context.fillStyle = "white";
    context.globalAlpha = 1;

    if (title.indexOf("/") !== -1) {
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10,
        textWidth / 2,
        14
      );
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10 + subtypeFontSize,
        textWidth / 2,
        14
      );
    } else {
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10,
        textWidth,
        14
      );
    }
    context.fill();
    context.save();

    context.globalAlpha = 1;
    context.fillStyle = "black";

    context.font = "normal bold " + subtypeFontSize + "px Droid";
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
      context.fillText(title, x(boxCords.left), y(boxCords.top) + height / 2);
    }

    context.fill();
    context.stroke();
    context.save();

    /*  if (markers[title]) {
      context.fillStyle = "white";
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10 + subtypeFontSize,
        40,
        geneFontSize * markers[title].length + 2
      );

      context.fillStyle = "grey";
      context.font = "normal bold " + geneFontSize + "px Droid";
      context.fill();
      markers[title].map((gene, index) => {
        const yCord =
          y(boxCords.top) + height / 2 + geneFontSize * index + subtypeFontSize;
        context.fillText(gene, x(boxCords.left), yCord);
        context.fill();
        context.stroke();
      });
    }*/
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
          <svg id="legend" style={{ float: "right", width: 600 }} />
        </div>
      </div>
    </div>
  );
};

export default SubtypeUmap;
