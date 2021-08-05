import React, { useState } from "react";
import * as d3 from "d3";
import _ from "lodash";

import { Grid } from "@material-ui/core";
import { useCanvas } from "@shahlab/planetarium";

const PADDING = 10;

const LEGEND_WIDTH = 0;
const AXIS_SPACE = 20;
const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";

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
  // if (highlighted) {
  //   data
  //     .filter((row) => row[subsetParam] === highlighted)
  //     .forEach((point) => {
  //       context.fillStyle = colorScale(point[subsetParam]);
  //       context.beginPath();
  //       context.arc(
  //         xScale(point[xParam]),
  //         yScale(point[yParam]),
  //         POINT_RADIUS,
  //         0,
  //         Math.PI * 2,
  //         true
  //       );
  //       context.fill();
  //     });
  // }
};

const drawUMAPAxis = (context, chartHeight, xParam, yParam) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  const START_X = AXIS_SPACE / 2;
  const START_Y = chartHeight + AXIS_SPACE / 2;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.lineWidth = 1;
  context.lineCap = "butt";
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
    d3.scaleSequential(d3.interpolateViridis).domain([0, subsetMax]);

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
      drawUMAPAxis(context, canvasHeight - AXIS_SPACE, xParam, yParam);

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
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <canvas ref={canvasRef} />
      </Grid>
      {/* <Grid item>
        <VerticalLegend
          title={subsetParam}
          width={LEGEND_WIDTH}
          labels={subsetLabels}
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
      </Grid> */}
    </Grid>
  );
};

export default UMAP;
