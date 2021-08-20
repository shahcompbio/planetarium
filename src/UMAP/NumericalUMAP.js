import React, { useState } from "react";
import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@material-ui/core";
import { useCanvas } from "../utils/useCanvas";
import Legend from "../Legend/Vertical";
import drawAxis from "./utils/drawAxis";

const PADDING = 10;

const LEGEND_WIDTH = 70;
const AXIS_SPACE = 20;

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;

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
  if (highlighted) {
  }
  data.forEach((point) => {
    context.fillStyle =
      point.hasOwnProperty(subsetParam) && highlighted.includes(point[idParam])
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
  xParam,
  yParam,
  subsetParam,
  idParam = "id",
  highlighted = null,
  disable = false,
  colorScale = null,
  onLasso = (data) => {},
}) => {
  const [lassoPolys, setLassoPolys] = useState([]);

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

  const subsetData = data
    .filter((d) => d.hasOwnProperty(subsetParam))
    .map((d) => parseFloat(d[subsetParam]));
  const subsetMax = Math.max(...subsetData);
  const subsetColors =
    colorScale ||
    d3.scaleSequential(d3.interpolateViridis).domain([0, subsetMax]).nice();

  const getHighlighted = (data, polys) => {
    const polyFiltered =
      polys.length > 0
        ? data.filter((datum) =>
            isPointInPoly(polys, xScale(datum[xParam]), yScale(datum[yParam]))
          )
        : data;

    return polyFiltered;
  };

  const addLassoHandler = (
    canvas,
    setLassoPolys,
    data,
    disable,
    highlighted
  ) => {
    const context = canvas.getContext("2d");

    const disableLasso = disable || highlighted !== null;
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
        const lassoedData = getHighlighted(data, polys);
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

      const highlightedIDs =
        highlighted !== null
          ? highlighted.map((datum) => datum[idParam])
          : getHighlighted(data, lassoPolys).map((datum) => datum[idParam]);

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

      drawLasso(context, lassoPolys);

      addLassoHandler(canvas, setLassoPolys, data, disable, highlighted);
    },
    canvasWidth,
    canvasHeight,
    [data, lassoPolys, disable, subsetParam, highlighted]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }} alignItems="center">
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      <Grid item style={{ paddingLeft: "40px" }}>
        <Legend
          title={subsetParam}
          width={LEGEND_WIDTH}
          height={height / 2}
          colorScale={subsetColors}
          ticks={10}
        />
      </Grid>
    </Grid>
  );
};

export default UMAP;
