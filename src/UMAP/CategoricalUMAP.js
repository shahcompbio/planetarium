import React, { useState } from "react";
import * as d3 from "d3";
import _ from "lodash";

import { quantileSorted } from "d3-array";

import { Grid } from "@mui/material";
import { useCanvas } from "../utils/useCanvas";
import VerticalLegend from "../Legend/Vertical";
import isHighlighted from "../utils/isHighlighted";
import drawAxis from "./utils/drawAxis";

const PADDING = 10;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const LABEL_FONT = "500 12px Helvetica";

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
  "#ffd470",
];

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;

const PERCENTILE_RANGE = [0.25, 0.75];

const LASSO_COLOR = "#c0392b";

const drawPoints = (
  context,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  idParam,
  highlighted,
  colorScale
) => {
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;
  data.forEach((point) => {
    context.fillStyle = highlighted.includes(point[idParam])
      ? colorScale(point[subsetParam])
      : NULL_POINT_COLOR;

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
  if (highlighted) {
    data
      .filter((row) => row[subsetParam] === highlighted)
      .forEach((point) => {
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

const drawLasso = (context, polys) => {
  context.globalAlpha = 0.2;
  context.fillStyle = "black";
  context.fill();

  const lassoPath = polys
    .map((poly, idx) => `${idx === 0 ? "M" : " L"} ${poly[0]} ${poly[1]}`)
    .reduce((str, poly) => str + poly, "");

  const lassoSvgPath = new Path2D(lassoPath);
  lassoSvgPath.closePath();

  context.fill(lassoSvgPath, "evenodd");
  context.closePath();
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

  subsetValues.forEach((subset) => {
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
      yScale(yMax),
    ];
    const width = Math.abs(x2 - x1);
    const height = Math.abs(y2 - y1);

    context.font = LABEL_FONT;
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
    .map((datum) => parseFloat(datum[xParam]))
    .sort((a, b) => a - b);
  const yValues = data
    .map((datum) => parseFloat(datum[yParam]))
    .sort((a, b) => a - b);

  const [xMin, xMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(xValues, range)
  );
  const [yMin, yMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(yValues, range)
  );

  const xiqr = Math.abs(xMax - xMin) * 1.5;
  const yiqr = Math.abs(yMax - yMin) * 1.5;

  return data.filter((datum) => {
    const x = datum[xParam];
    const y = datum[yParam];

    return (
      xMin - xiqr <= x &&
      x <= xMax + xiqr &&
      yMin - yiqr <= y &&
      y <= yMax + yiqr
    );
  });
};
const getBoxBounds = (data, xParam, yParam) => {
  const xValues = data.map((datum) => datum[xParam]);
  const yValues = data.map((datum) => datum[yParam]);

  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const yMin = Math.min(...yValues);
  const yMax = Math.max(...yValues);

  return { xMin, xMax, yMin, yMax };
};

function isPointInPoly(poly, x, y) {
  for (var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
    ((poly[i][1] <= y && y < poly[j][1]) ||
      (poly[j][1] <= y && y < poly[i][1])) &&
      x <
        ((poly[j][0] - poly[i][0]) * (y - poly[i][1])) /
          (poly[j][1] - poly[i][1]) +
          poly[i][0] &&
      (c = !c);
  return c;
}

const UMAP = ({
  width,
  height,
  data,
  subset = null,
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  disable = false,
  colorScale = null,
  onLasso = (data) => {},
  onLegendHover = (value) => {},
  onLegendClick = (value) => {},
  MoreInfoComponent = () => null,
}) => {
  const [highlightedSubset, setHighlightedSubset] = useState(null);
  const [selectedSubset, setSelectedSubset] = useState(null);

  const [lassoPolys, setLassoPolys] = useState([]);

  const highlightedOverall = highlightedSubset || subset || selectedSubset;

  const canvasWidth = width - LEGEND_WIDTH;
  const canvasHeight = height;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

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

  const subsetColors =
    colorScale ||
    d3
      .scaleOrdinal()
      .domain(subsetValues)
      .range(
        COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
      );

  const getHighlighted = (data, polys, subsetValue) => {
    const subsetFiltered = data.filter((datum) =>
      isHighlighted(datum[subsetParam], subsetValue)
    );

    const polyFiltered =
      polys.length > 0
        ? subsetFiltered.filter((datum) =>
            isPointInPoly(polys, xScale(datum[xParam]), yScale(datum[yParam]))
          )
        : subsetFiltered;

    return polyFiltered;
  };

  const addLassoHandler = (
    canvas,
    setLassoPolys,
    data,
    highlightedOverall,
    disable
  ) => {
    const context = canvas.getContext("2d");

    const disableLasso = disable || highlightedOverall !== null;
    let polys = [];
    d3.select(canvas)
      .on("mousemove", (d, i, e) => {
        if (disableLasso) {
          return;
        }
        if (d3.event.buttons === 1) {
          const poly = d3.mouse(e[0]);
          polys.push(poly);
          context.lineTo(poly[0], poly[1]);
          context.stroke();
        }
      })
      .on("mousedown", function mousedown() {
        if (disableLasso) {
          return;
        }
        polys = [];
        context.lineCap = "round";
        context.strokeStyle = LASSO_COLOR;
        context.lineWidth = 3;

        context.restore();
        context.beginPath();
      })
      .on("mouseup", () => {
        if (disableLasso) {
          return;
        }
        setLassoPolys(polys);
        const lassoedData = getHighlighted(data, polys, highlightedOverall);
        onLasso(lassoedData);
      });
  };

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawAxis(
        context,
        AXIS_SPACE / 2,
        canvasHeight - AXIS_SPACE / 2,
        xParam,
        yParam
      );

      const highlightedIDs = getHighlighted(
        data,
        lassoPolys,
        highlightedOverall
      ).map((datum) => datum[idParam]);

      drawPoints(
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        idParam,
        highlightedIDs,
        subsetColors
      );

      drawSubsetLabels(
        context,
        subsetGroups,
        xScale,
        yScale,
        xParam,
        yParam,
        highlightedOverall
      );

      drawLasso(context, lassoPolys);

      addLassoHandler(canvas, setLassoPolys, data, highlightedOverall, disable);
    },
    canvasWidth,
    canvasHeight,
    [data, highlightedOverall, lassoPolys, disable]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid item>
        <VerticalLegend
          title={subsetParam}
          width={LEGEND_WIDTH}
          ticks={subsetValues}
          colorScale={subsetColors}
          disable={disable || lassoPolys.length > 0}
          onHover={(value) => {
            setHighlightedSubset(value);
            onLegendHover(value);
          }}
          onClick={(value) => {
            setSelectedSubset(value);
            onLegendClick(value);
          }}
        />
        <MoreInfoComponent />
      </Grid>
    </Grid>
  );
};

export default UMAP;
