import React from "react";
import * as d3 from "d3";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useCanvas } from "../utils/useCanvas";
import Info from "../../Info/Info.js";
import infoText from "../../Info/InfoText.js";

/*

Stacked horizontal bar chart of proportions (so they all add to 100%)

*/

const BAR_COLORS = [
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
  "#9E0142",
];

const LEGEND_HEIGHT = 50;

const PROP_AXIS_FONT = "normal 10px Helvetica";
const CAT_LABEL_SPACE = 150;
const CAT_LABEL_FONT = "normal 12px Helvetica";
const PADDING = 10;
const TITLE_HEIGHT = 30;

const LEGEND_SQUARE_LENGTH = 12;
const LEGEND_SQUARE_PADDING = 10;

const StackedBarProportion = ({ data, chartDim, barLabels, chartName }) => {
  const categoryValues = Object.keys(data).sort();
  const barValues =
    typeof barLabels[0] === "string"
      ? barLabels
      : barLabels.map((obj) => obj["value"]);
  const legendLabels =
    typeof barLabels[0] === "string"
      ? barLabels.map((str) => ({ value: str, label: str }))
      : barLabels;

  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - PADDING - PADDING - TITLE_HEIGHT;

  const chartWidth = canvasWidth - CAT_LABEL_SPACE;
  const chartHeight = canvasHeight - LEGEND_HEIGHT;

  const catScale = d3
    .scaleBand()
    .domain(categoryValues)
    .range([LEGEND_HEIGHT, LEGEND_HEIGHT + chartHeight])
    .paddingInner(0.03)
    .paddingOuter(0.25);

  const barPosScale = d3
    .scaleLinear()
    .domain([0, 1])
    .range([PADDING, PADDING + chartWidth]);

  const barSizeScale = d3
    .scaleLinear()
    .domain([0, 1])
    .range([0, chartWidth]);

  const colors = d3
    .scaleOrdinal()
    .domain(barValues)
    .range(BAR_COLORS.slice(0, barValues.length));

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      drawBars(
        context,
        data,
        categoryValues,
        barValues,
        catScale,
        barPosScale,
        barSizeScale,
        colors
      );
      drawLabels(
        context,
        categoryValues,
        catScale,
        barPosScale,
        chartWidth + PADDING + 2,
        chartHeight + LEGEND_HEIGHT
      );

      drawLegend(context, legendLabels, colors, canvasWidth);
    },

    canvasWidth,
    canvasHeight,
    []
  );

  return (
    <Paper
      style={{
        margin: 10,
        padding: PADDING,
        height: chartDim["height"],
        width: chartDim["width"],
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="stretch"
      >
        <Grid
          item
          style={{
            textAlign: "right",
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
        <Grid item>
          <canvas ref={ref} />
        </Grid>
      </Grid>
    </Paper>
  );
};

const drawBars = (
  context,
  data,
  categoryValues,
  barValues,
  catScale,
  barPosScale,
  barSizeScale,
  colors
) => {
  categoryValues.forEach((cValue) => {
    const categoryData = data[cValue];
    const total = Object.values(categoryData).reduce((sum, x) => sum + x, 0);

    var currPos = barPosScale(0);
    const yPos = catScale(cValue);

    barValues.forEach((bValue) => {
      if (categoryData.hasOwnProperty(bValue)) {
        context.fillStyle = colors(bValue);
        const barSize = barSizeScale(categoryData[bValue] / total);

        context.fillRect(currPos, yPos, barSize, catScale.bandwidth());

        currPos += barSize;
      }
    });
  });
};

const drawLabels = (
  context,
  categoryValues,
  catScale,
  barScale,
  xAxisPos,
  yAxisPos
) => {
  const values = barScale.ticks(10);
  context.font = PROP_AXIS_FONT;
  context.fillStyle = "black";
  context.textAlign = "center";
  context.lineWidth = 1;
  context.textBaseline = "bottom";

  values.forEach((value) => {
    context.fillText(value * 100, barScale(value), yAxisPos);
  });

  context.font = CAT_LABEL_FONT;
  context.textAlign = "left";
  context.textBaseline = "middle";

  categoryValues.forEach((cValue) => {
    context.fillText(
      cValue,
      xAxisPos,
      catScale(cValue) + catScale.bandwidth() / 2
    );
  });
};

const drawLegend = (context, barValues, colors, canvasWidth) => {
  context.textAlign = "center";
  context.textBaseline = "top";
  context.font = PROP_AXIS_FONT;

  const START_X = canvasWidth - 30;
  const START_Y = 20;

  barValues.reverse().forEach((labelObj, index) => {
    const { value, label } = labelObj;
    context.fillStyle = colors(value);

    context.fillRect(
      START_X - (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_PADDING) * (index + 1),
      START_Y,
      LEGEND_SQUARE_LENGTH,
      LEGEND_SQUARE_LENGTH
    );

    context.fillStyle = "black";
    context.fillText(
      label,
      START_X -
        (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_PADDING) * (index + 1) +
        LEGEND_SQUARE_LENGTH / 2,
      START_Y + LEGEND_SQUARE_LENGTH + 3
    );
  });
};

export default StackedBarProportion;
