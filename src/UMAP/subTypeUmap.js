import React, { useState, useEffect, useRef, useCallback } from "react";
import * as d3 from "d3";
import _ from "lodash";
import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";

import { useDashboardState } from "../PlotState/dashboardState";
import { drawPoint, clearAll } from "./Umap.js";
import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";

const SubtypeUmap = ({
  chartName,
  data,
  metadata,
  chartDim,
  selectedSubtype,
  hoveredSubtype,
  setSelectedSubtype
}) => {
  const [{ xParam, yParam, subtypeParam, fontSize }] = useDashboardState();

  const [drawReady, setDrawReady] = useState(false);
  const [context, saveContext] = useState(null);

  const yData = data.map(d => parseFloat(d[yParam]));
  const xData = data.map(d => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  // X axis
  const x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"]]);

  // Y axis
  const y = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([chartDim["chart"]["y1"], chartDim["chart"]["y2"]]);

  const subTypes = _.groupBy(data, subtypeParam);

  const types = Object.keys(subTypes);
  var colors = d3
    .scaleOrdinal()
    .domain([...types])
    .range([
      "#5E4FA2",
      "#3288BD",
      "#66C2A5",
      "#FEE08B",
      "#FDAE61",
      "#F46D43",
      "#D53E4F",
      "#9E0142"
    ]);

  useEffect(() => {
    if (data.length > 0) {
      init(data, chartDim);
    }
  }, [data]);

  useEffect(() => {
    if (context) {
      const selection =
        hoveredSubtype !== null
          ? hoveredSubtype
          : selectedSubtype !== null
          ? selectedSubtype
          : null;
      console.log("sle", selection);
      if (selection !== null) {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y, selection);
      } else {
        clearAll(context, chartDim);
        context.beginPath();
        reDraw(data, chartDim, context, x, y);
      }
    }
  }, [selectedSubtype, hoveredSubtype, context]);

  function init(data, chartDim) {
    var canvas = d3.select("#subTypeUmapCanvas");

    var currContext = canvasInit(canvas, chartDim.width, chartDim.height);

    currContext.fillStyle = "white";
    currContext.fillRect(0, 0, chartDim.width, chartDim.height);
    saveContext(currContext);
  }

  function drawPoints(data, chartDim, context, x, y, selection) {
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";

    data.forEach(point => {
      const fill =
        selection && selection !== point[subtypeParam]
          ? "#e8e8e8"
          : colors(point[subtypeParam]);
      drawPoint(context, point, fill, true, x, y, xParam, yParam);
    });
    if (selection) {
      data
        .filter(point => point[subtypeParam] === selection)
        .forEach(point => {
          drawPoint(
            context,
            point,
            colors(point[subtypeParam]),
            true,
            x,
            y,
            xParam,
            yParam
          );
        });
    }
  }
  function appendSubtypeLabels(
    subTypes,
    x,
    y,
    yParam,
    xParam,
    context,
    colors,
    selection
  ) {
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
    /*  const avgBoxSize =
      allBoxes
        .map(box => x(box["box"]["right"]) - x(box["box"]["left"]))
        .reduce((a, b) => a + b) / allBoxes.length;
*/
    //if the bounding box is larger than 1.5 the  avg do not add label
    allBoxes.forEach(box => {
      //    if (
      //      !(x(box["box"]["right"]) - x(box["box"]["left"]) * 1.5 > avgBoxSize)
      //    ) {
      drawBoundingBox(context, box, x, y, colors, selection);
      //    }
    });
  }
  function reDraw(data, chartDim, context, x, y, selection) {
    drawAxis(context, x, y, chartDim["chart"]);
    drawPoints(data, chartDim, context, x, y, selection);
    appendSubtypeLabels(
      subTypes,
      x,
      y,
      yParam,
      xParam,
      context,
      colors,
      selection
    );
    appendLegend(colors, types, context, x, y);
  }
  function appendLegend(colors, subTypes, context, x, y) {
    const mouseInteractions = element =>
      element
        .on("mouseenter", function(d) {
          d3.event.stopPropagation();
          setSelectedSubtype({
            hover: d,
            selected: selectedSubtype
          });
        })
        .on("mousedown", function(d, i) {
          d3.event.stopPropagation();
          setSelectedSubtype({
            hover: null,
            selected: d
          });
        })
        .on("mouseout", function(event, d) {
          d3.event.stopPropagation();
          setSelectedSubtype({
            hover: null,
            selected: selectedSubtype
          });
        });

    var legend = d3.select("#subTypeUmapLegend");
    legend.selectAll("*").remove();
    const legendRect = legend.selectAll("rect").data(subTypes);

    const legendEnter = legendRect
      .append("rect")
      .attr("width", fontSize.legendSquare)
      .attr("height", fontSize.legendSquare)
      .attr("x", function(d) {
        return chartDim["legend"]["x1"] + 20;
      })
      .attr("y", function(d, i) {
        return chartDim["legend"].y1 + i * 20;
      })
      .attr("fill", function(d) {
        return colors(d);
      })
      .attr("cursor", "pointer");

    legendRect
      .enter()
      .append("rect")
      .attr("width", fontSize.legendSquare)
      .attr("height", fontSize.legendSquare)
      .attr("x", function(d) {
        return chartDim["legend"]["x1"] + 20;
      })
      .attr("y", function(d, i) {
        return chartDim["legend"].y1 + i * 20;
      })
      .attr("fill", function(d) {
        return colors(d);
      })
      .attr("cursor", "pointer")
      .call(mouseInteractions);

    const legendText = legend.selectAll("text").data(subTypes);
    const legendTextEnter = legendText
      .append("text")
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 13 + 20;
      })
      .attr("y", function(d, i) {
        return chartDim["legend"].y1 + i * 20 + 4;
      })
      .attr("dy", ".35em")
      .text(function(d) {
        return d;
      })
      .attr("font", "Helvetica")
      .attr("font-size", fontSize.legendFontSize)
      .attr("font-weight", "500")
      .attr("fill", function(d) {
        return "#000000";
      })
      .attr("cursor", "pointer");
    //  .call(mouseInteractions);

    legendText
      .enter()
      .append("text")
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 13 + 20;
      })
      .attr("y", function(d, i) {
        return chartDim["legend"].y1 + i * 20 + 4;
      })
      .attr("dy", ".35em")
      .text(function(d) {
        return d;
      })
      .attr("font", "Helvetica")
      .attr("font-size", fontSize.legendFontSize)
      .attr("font-weight", "500")
      .attr("fill", function(d) {
        return "#000000";
      })
      .attr("cursor", "pointer")
      .call(mouseInteractions);
  }

  function filterOutliers(coords, quantiles) {
    const sortedCollection = coords.slice().sort((a, b) => a - b); //copy array fast and sort

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
  function drawBoundingBox(context, box, x, y, colors, selection) {
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
    const subtypeFontSize = title.length === 1 ? 18 : 12;
    context.font = "500 " + subtypeFontSize + "px Helvetica";

    const textWidth = context.measureText(title).width + 2;

    context.fillStyle = "white";
    context.globalAlpha = selection ? (selection === title ? 1 : 0.2) : 1;

    if (title.indexOf("/") !== -1) {
      const firstTextWidth = context.measureText(title.split("/")[0]).width + 4;
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10,
        firstTextWidth,
        subtypeFontSize
      );
      const secondTextWidth = context.measureText(title.split("/")[1]).width;
      context.fillRect(
        x(boxCords.left) - 2,
        y(boxCords.top) + height / 2 - 10 + subtypeFontSize,
        secondTextWidth,
        subtypeFontSize
      );
    } else {
      context.fillRect(
        x(boxCords.left) + width / 2,
        y(boxCords.top) - height / 2 - 12,
        textWidth,
        subtypeFontSize + 2
      );
    }
    context.fill();
    context.save();

    context.globalAlpha = selection ? (selection === title ? 1 : 0.2) : 1;
    context.fillStyle = "black";

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
    context.globalAlpha = 1;
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
    <div class="card" style={{ margin: 10 }}>
      <div
        class="container"
        style={{
          width: chartDim["width"] + 200,
          height: chartDim["height"],
          position: "relative"
        }}
      >
        <div class="row">
          <div
            class="col-9"
            id="subTypeUmap"
            style={{
              pointerEvents: "all",
              display: "flex",
              paddingRight: 0
            }}
          >
            <canvas id="subTypeUmapCanvas" />
          </div>
          <div class="col-3" style={{ paddingLeft: 0 }}>
            <div
              class="card-title"
              style={{
                marginTop: chartDim["chart"]["x1"],
                width: "100%",
                height: 80,
                paddingTop: 40,
                paddingLeft: -50,
                textAlign: "left"
              }}
            >
              {infoText[chartName]["title"] + "    "}

              <Info name={chartName} direction="s" />
            </div>
            <div style={{ marginLeft: -50, height: 250 }}>
              <svg id="subTypeUmapLegend" style={{ float: "right" }} />
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default SubtypeUmap;
