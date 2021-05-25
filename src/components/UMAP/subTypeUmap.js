import React, { useState } from "react";
import * as d3 from "d3";
import { quantileSorted } from "d3-array";
import _ from "lodash";
import Info from "../../Info/Info.js";
import infoText from "../../Info/InfoText.js";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useDashboardState } from "../../PlotState/dashboardState";

import { useCanvas } from "../utils/useCanvas";
import { useD3 } from "../utils/useD3";

const PADDING = 10;
const TITLE_HEIGHT = 30;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const AXIS_COLOR = "#000000";

const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470"
];

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;
const PERCENTILE_RANGE = [0.25, 0.75];
const LEGEND_SQUARE_LENGTH = 10;
const LEGEND_SQUARE_SPACING = 8;

const DataWrapper = ({
  chartName,
  data,
  chartDim,
  selectedSubtype,
  hoveredSubtype,
  setSelectedSubtype
}) => {
  const [{ xParam, yParam, subtypeParam }] = useDashboardState();

  const setHighlighted = (event, value) => {
    if (event === "mouseenter") {
      setSelectedSubtype({ hover: value });
    } else if (event === "mousedown") {
      setSelectedSubtype({
        hover: null,
        selected: value
      });
    } else if (event === "mouseout") {
      setSelectedSubtype({
        hover: null
      });
    }
  };
  return (
    <UMAP
      chartDim={chartDim}
      chartName={chartName}
      data={data}
      highlighted={hoveredSubtype || selectedSubtype}
      xParam={xParam}
      yParam={yParam}
      subsetParam={subtypeParam}
      setHighlighted={setHighlighted}
    />
  );
};

const UMAP = ({
  chartDim,
  chartName,
  data,
  highlighted,
  xParam,
  yParam,
  subsetParam,
  setHighlighted
}) => {
  var mousePos = { x: 0, y: 0 };
  var lassPath = "";
  var polyList = [];
  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - LEGEND_WIDTH - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const yData = data.map(d => parseFloat(d[yParam]));
  const xData = data.map(d => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([PADDING, PADDING + chartWidth]);

  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([PADDING, PADDING + chartHeight]);

  const subsetGroups = _.groupBy(data, subsetParam);
  const subsetValues = Object.keys(subsetGroups).sort();

  const subsetColors = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(
      COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
    );

  const canvasRef = useCanvas(
    canvas => {
      const context = canvas.getContext("2d");
      saveContext(context);
      drawUMAPAxis(context, chartHeight, xParam, yParam);
      drawPoints(
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        highlighted,
        subsetColors
      );

      drawSubsetLabels(
        context,
        subsetGroups,
        xScale,
        yScale,
        xParam,
        yParam,
        highlighted
      );
      appendEventListenersToCanvas(context, data);
    },
    canvasWidth,
    canvasHeight,
    [highlighted]
  );
  //taken from
  //https://gist.github.com/maxogden/574870
  function isPointInPoly(poly, pt) {
    for (var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
      ((poly[i][1] <= pt.y && pt.y < poly[j][1]) ||
        (poly[j][1] <= pt.y && pt.y < poly[i][1])) &&
        pt.x <
          ((poly[j][0] - poly[i][0]) * (pt.y - poly[i][1])) /
            (poly[j][1] - poly[i][1]) +
            poly[i][0] &&
        (c = !c);
    return c;
  }
  function setMousePosition(e, boundingRect) {
    mousePos.x = e.clientX - boundingRect.left;
    mousePos.y = e.clientY - boundingRect.top;
  }
  function setD3MousePosition(event, boundingRect) {
    mousePos.x = event.pageX - boundingRect.left;
    mousePos.y = event.pageY - boundingRect.top;
  }
  function drawLasso(context, boundingRect, x, y) {
    // mouse left button must be pressed
    if (d3.event.buttons !== 1) return;

    context.lineCap = "round";
    context.strokeStyle = "#c0392b";
    context.lineWidth = 3;

    context.moveTo(mousePos.x, mousePos.y);
    polyList = [...polyList, [mousePos.x, mousePos.y]];
    lassPath +=
      lassPath === ""
        ? "M " + mousePos.x + " " + mousePos.y + " "
        : "L " + mousePos.x + " " + mousePos.y + " ";

    //setMousePosition(e, boundingRect);
    setD3MousePosition(d3.event, boundingRect);
    context.lineTo(mousePos.x, mousePos.y);
    lassPath += "L " + mousePos.x + " " + mousePos.y + " ";
    polyList = [...polyList, [mousePos.x, mousePos.y]];

    context.stroke();
  }
  const appendEventListenersToCanvas = (context, data) => {
    const scatterSelection = d3.select("#umapSelection");
    console.log(data);
    d3.select("#umapCanvas")
      .on("mousemove", function mousemove(e) {
        drawLasso(
          context,
          this.getBoundingClientRect(),
          d3.event.pageX,
          d3.event.pageY
        );
      })
      .on("mousedown", function mousedown() {
        context.restore();
        context.beginPath();
        polyList = [];
        setD3MousePosition(d3.event, this.getBoundingClientRect());
        //setMousePosition(e, this.getBoundingClientRect());
      })
      .on("mouseenter", function mouseenter(e) {
        setD3MousePosition(d3.event, this.getBoundingClientRect());
        //setMousePosition(e, this.getBoundingClientRect());
      })
      .on("mouseup", function mouseup(e) {
        context.save();

        context.strokeWidth = 1;
        context.globalAlpha = 0.5;
        context.strokeStyle = "purple";
        context.fillStyle = "black";

        context.fill();
        var lassoSvgPath = new Path2D(lassPath);
        lassoSvgPath.closePath();
        context.fill(lassoSvgPath, "evenodd");
        context.fill();
        context.closePath();

        var selectedNodes = [];
        console.log(data);
        data.map(point => {
          const cords = {
            x: xScale(parseFloat(point[xParam])),
            y: yScale(parseFloat(point[yParam]))
          };
          if (isPointInPoly(polyList, cords)) {
            selectedNodes = [...selectedNodes, point];
          }
        });
        const lassoIDs = selectedNodes.reduce((final, curr) => {
          final[curr["cell_id"]] = true;
          return final;
        }, {});
        console.log(selectedNodes);
        /*scatterSelection.selectAll("rect").each(function(d) {
            const xCord = d3.select(this).attr("x");
            const yCord = d3.select(this).attr("y");
            if (isPointInPoly(polyList, { x: xCord, y: yCord })) {
              selectedNodes = [...selectedNodes, d["heatmapOrder"]];
            }
          });*/
        context.restore();

        lassPath = "";
        if (selectedNodes.length > 0) {
          drawPoints(
            context,
            data,
            xScale,
            yScale,
            xParam,
            yParam,
            subsetParam,
            highlighted,
            subsetColors,
            lassoIDs
          );
        }
      });
  };
  const svgRef = useD3(svg => {
    drawLegend(
      svg,
      subsetValues,
      subsetColors,
      canvasHeight,
      highlighted,
      setHighlighted
    );
  }, []);

  return (
    <Paper
      style={{
        margin: 10,
        padding: "10px 0px",
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
        <Grid container direction="row" style={{ padding: 0 }}>
          <Grid item>
            <canvas ref={canvasRef} id="umapCanvas" />
          </Grid>
          <Grid item>
            <svg ref={svgRef} id="umapSelection" />
          </Grid>
        </Grid>
      </Grid>
    </Paper>
  );
};

const drawPoints = (
  context,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  highlightedSubtype,
  colorScale,
  lassoIDs
) => {
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;

  data.forEach(point => {
    if (lassoIDs) {
      context.fillStyle = lassoIDs[point["cell_id"]]
        ? colorScale(point[subsetParam])
        : NULL_POINT_COLOR;
    } else {
      context.fillStyle = isHighlighted(point[subsetParam], highlightedSubtype)
        ? colorScale(point[subsetParam])
        : NULL_POINT_COLOR;
    }

    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      POINT_RADIUS,
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });
  if (highlightedSubtype || lassoIDs) {
    data
      .filter(
        row =>
          row[subsetParam] === highlightedSubtype ||
          (lassoIDs && lassoIDs[row["cell_id"]])
      )
      .forEach(point => {
        context.fillStyle = colorScale(point[subsetParam]);
        context.beginPath();
        context.arc(
          xScale(point[xParam]),
          yScale(point[yParam]),
          POINT_RADIUS,
          0,
          Math.PI * 2,
          true
        );
        context.fill();
      });
  }
};

const drawUMAPAxis = (context, chartHeight, xParam, yParam) => {
  context.beginPath();

  const START_X = AXIS_SPACE / 2;
  const START_Y = chartHeight + AXIS_SPACE / 2;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.moveTo(START_X, START_Y);
  context.lineTo(START_X, START_Y - 50);
  context.stroke();

  context.beginPath();
  context.moveTo(START_X, START_Y);
  context.lineTo(START_X + 50, START_Y);
  context.stroke();

  context.textAlign = "left";
  context.textBaseline = "middle";
  context.fillText(xParam, START_X + 52, START_Y);
  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText(yParam, -(START_Y - 52), START_X);
  context.restore();
};

const drawSubsetLabels = (
  context,
  subsetGroups,
  xScale,
  yScale,
  xParam,
  yParam,
  highlighted
) => {
  const subsetValues = Object.keys(subsetGroups);

  subsetValues.forEach(subset => {
    const subsetData = subsetGroups[subset];

    const filteredData = filterOutliers(subsetData, xParam, yParam);
    const { xMin, xMax, yMin, yMax } = getBoxBounds(
      filteredData,
      xParam,
      yParam
    );

    const [x1, x2, y1, y2] = [
      xScale(xMin),
      xScale(xMax),
      yScale(yMin),
      yScale(yMax)
    ];
    const width = Math.abs(x2 - x1);
    const height = Math.abs(y2 - y1);

    context.font = "500 12px Helvetica";
    const textWidth = context.measureText(subset).width;
    context.globalAlpha = isHighlighted(subset, highlighted) ? 0.8 : 0.2;
    context.fillStyle = "white";
    context.fillRect(
      x1 + (width - textWidth) / 2 - 1,
      y1 - height / 2 - 7,
      textWidth + 2,
      14
    );

    context.globalAlpha = isHighlighted(subset, highlighted) ? 1 : 0.2;
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.fillStyle = "black";
    context.fillText(subset, x1 + width / 2, y1 - height / 2);
  });
};

const filterOutliers = (data, xParam, yParam) => {
  const xValues = data
    .map(datum => parseFloat(datum[xParam]))
    .sort((a, b) => a - b);
  const yValues = data
    .map(datum => parseFloat(datum[yParam]))
    .sort((a, b) => a - b);

  const [xMin, xMax] = PERCENTILE_RANGE.map(range =>
    quantileSorted(xValues, range)
  );
  const [yMin, yMax] = PERCENTILE_RANGE.map(range =>
    quantileSorted(yValues, range)
  );

  const xiqr = Math.abs(xMax - xMin) * 1.5;
  const yiqr = Math.abs(yMax - yMin) * 1.5;

  return data.filter(datum => {
    const x = datum[xParam];
    const y = datum[yParam];

    return (
      xMin - xiqr <= x &&
      x <= xMax + xiqr &&
      (yMin - yiqr <= y && y <= yMax + yiqr)
    );
  });
};

const getBoxBounds = (data, xParam, yParam) => {
  const xValues = data.map(datum => datum[xParam]);
  const yValues = data.map(datum => datum[yParam]);

  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const yMin = Math.min(...yValues);
  const yMax = Math.max(...yValues);

  return { xMin, xMax, yMin, yMax };
};

const drawLegend = (
  svg,
  subsetValues,
  colors,
  chartHeight,
  highlighted,
  setHighlighted
) => {
  const mouseEvents = element =>
    element.call(element =>
      element
        .on("mouseenter", function(d) {
          d3.event.stopPropagation();
          setHighlighted("mouseenter", d);
        })
        .on("mousedown", function(d, i) {
          d3.event.stopPropagation();
          setHighlighted("mousedown", d);
        })
        .on("mouseout", function(d, i) {
          d3.event.stopPropagation();
          setHighlighted("mouseout", d);
        })
    );
  svg.attr("width", LEGEND_WIDTH).attr("height", chartHeight);

  const subsets = svg.selectAll("g").data(subsetValues);

  subsets
    .enter()
    .append("rect")
    .merge(subsets)
    .attr("width", LEGEND_SQUARE_LENGTH)
    .attr("height", LEGEND_SQUARE_LENGTH)
    .attr("x", 5)
    .attr("y", (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5)
    .attr("fill", d => colors(d))
    .call(mouseEvents);

  subsets
    .enter()
    .append("text")
    .merge(subsets)
    .attr("alignment-baseline", "hanging")
    .attr("text-align", "left")
    .attr("font", "Helvetica")
    .attr("font-weight", "500")
    .attr("font-size", "12px")
    .attr("fill", "#000000")
    .attr("x", LEGEND_SQUARE_LENGTH + 10)
    .attr("y", (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5)
    .text(d => d)
    .call(mouseEvents);
};

const isHighlighted = (datumValue, highlighted) =>
  highlighted === null || datumValue === highlighted;

export default DataWrapper;
