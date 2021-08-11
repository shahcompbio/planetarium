import React from "react";
import PropTypes from "prop-types";
import * as d3 from "d3";

import { useCanvas } from "../utils/useCanvas";
import drawAxis from "../utils/canvas/drawAxis";
import { isValueHighlighted } from "../utils/isHighlighted";

/*

Proportional stacked horizontal bar chart, divided by subgroups

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

const LEGEND_HEIGHT = 70;

const PROP_AXIS_FONT = "normal 10px Helvetica";
const CAT_LABEL_SPACE = 150;
const CAT_LABEL_FONT = "normal 12px Helvetica";
const PADDING = 10;
const LABEL_PADDING = 20;

const LEGEND_SQUARE_LENGTH = 12;
const LEGEND_SQUARE_PADDING = 10;

const StackedHorizontalBar = ({
  data,
  width = 400,
  height = 400,
  barLabels,
  highlightedRow = null,
}) => {
  const categoryValues = Object.keys(data).sort();
  const barValues =
    typeof barLabels[0] === "string"
      ? barLabels
      : barLabels.map((obj) => obj["value"]);
  const legendLabels =
    typeof barLabels[0] === "string"
      ? barLabels.map((str) => ({ value: str, label: str }))
      : barLabels;

  const chartWidth = width - CAT_LABEL_SPACE;
  const chartHeight = height - LEGEND_HEIGHT - LABEL_PADDING;

  const catScale = d3
    .scaleBand()
    .domain(categoryValues)
    .range([LEGEND_HEIGHT, LEGEND_HEIGHT + chartHeight])
    .paddingInner(0.03);

  const barPosScale = d3
    .scaleLinear()
    .domain([0, 1])
    .range([PADDING, PADDING + chartWidth]);

  const barSizeScale = d3.scaleLinear().domain([0, 1]).range([0, chartWidth]);

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
        colors,
        highlightedRow
      );
      drawLabels(
        context,
        categoryValues,
        catScale,
        barPosScale,
        chartWidth + PADDING + 2,
        highlightedRow
      );

      drawLegend(context, legendLabels, colors, width);
    },

    width,
    height,
    [highlightedRow]
  );

  return <canvas ref={ref} />;
};

const drawBars = (
  context,
  data,
  categoryValues,
  barValues,
  catScale,
  barPosScale,
  barSizeScale,
  colors,
  highlightedRow
) => {
  categoryValues.forEach((cValue) => {
    const categoryData = data[cValue];
    const total = Object.values(categoryData).reduce((sum, x) => sum + x, 0);

    var currPos = barPosScale(0);
    const yPos = catScale(cValue);

    barValues.forEach((bValue) => {
      if (categoryData.hasOwnProperty(bValue)) {
        context.fillStyle = colors(bValue);
        context.globalAlpha = isValueHighlighted(cValue, highlightedRow)
          ? 1
          : 0.2;
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
  highlightedRow
) => {
  drawAxis({
    context,
    xScale: barScale,
    yScale: catScale,
    ticks: 10,
    orientation: "horizontal",
    gridlines: false,
    format: (tick) => tick * 100,
  });

  context.font = CAT_LABEL_FONT;
  context.textAlign = "left";
  context.textBaseline = "middle";

  categoryValues.forEach((cValue) => {
    context.globalAlpha = isValueHighlighted(cValue, highlightedRow) ? 1 : 0.2;
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
  context.globalAlpha = 1;
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

StackedHorizontalBar.propTypes = {
  /**
   * object of keys of rows with counts
   */
  data: PropTypes.object.isRequired,
  /**
   * width of plot
   */
  width: PropTypes.number,
  /**
   * height of plot
   */
  height: PropTypes.number,
  /**
   * List of row names
   */
  barLabels: PropTypes.arrayOf(
    PropTypes.oneOfType([PropTypes.string, PropTypes.object])
  ).isRequired,
  /**
   * Name of row to highlight
   */
  highlightedRow: PropTypes.string,
};
export default StackedHorizontalBar;
