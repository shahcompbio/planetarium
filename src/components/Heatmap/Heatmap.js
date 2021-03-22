import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useCanvas } from "../utils/useCanvas";

import Info from "../../Info/Info";
import infoText from "../../Info/InfoText";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

const HEATMAP_NULL_COLOR = "#eeeeee";
const HEATMAP_COLOR = ["#ffec8b", "#d91e18"];
const CELL_FONT = "normal 12px Helvetica";

const COLUMN_LABEL_SPACE = 150;
const ROW_LABEL_SPACE = 150;
const DEFAULT_LABEL_COLOR = "#000000";
const LABEL_FONT = "12px Helvetica";

const PADDING = 10;
const TITLE_HEIGHT = 30;

/*

This heatmap calculates it by total row value (or given total row value)

*/
const Heatmap = ({
  data,
  chartDim,
  column,
  row,
  highlightedColumn,
  highlightedRow,
  columnLabels,
  rowLabels,
  rowTotal,
}) => {
  const columnValues =
    columnLabels.map((col) => col["value"]) ||
    _.uniq(data.map((record) => record[column])).sort();
  const rowValues =
    rowLabels || _.uniq(data.map((record) => record[row])).sort();

  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - PADDING - PADDING - TITLE_HEIGHT;

  const chartWidth = canvasWidth - ROW_LABEL_SPACE;
  const chartHeight = canvasHeight - COLUMN_LABEL_SPACE;

  const columnScale = d3
    .scaleBand()
    .domain(columnValues)
    .range([0, chartWidth])
    .paddingInner(0.03);

  const rowScale = d3
    .scaleBand()
    .domain(rowValues)
    .range([COLUMN_LABEL_SPACE, chartHeight + COLUMN_LABEL_SPACE])
    .paddingInner(0.03);

  const groupedRow = _.groupBy(data, row);

  const freqMap = Object.keys(groupedRow).reduce(
    (currMap, rowName) => ({
      ...currMap,
      [rowName]: {
        ..._.countBy(groupedRow[rowName], column),
        total: rowTotal ? rowTotal[rowName] : groupedRow[rowName].length,
      },
    }),
    {}
  );
  const mostFreqCount = _.reduce(
    freqMap,
    (currMax, rowData, key) => {
      return Math.max(
        currMax,
        _.reduce(
          rowData,
          (columnMax, value, key) => {
            return key === "total" ? columnMax : Math.max(value, columnMax);
          },
          0
        )
      );
    },
    0
  );

  const heatmapColor = d3
    .scaleLinear()
    .range(HEATMAP_COLOR)
    .domain([0, mostFreqCount]);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawHeatmap(
        context,
        freqMap,
        rowValues,
        columnValues,
        rowScale,
        columnScale,
        highlightedRow,
        highlightedColumn,
        heatmapColor
      );

      const formattedColumnLabels = formatLabelData(
        columnLabels || columnValues
      );
      const formattedRowLabels = formatLabelData(rowLabels || rowValues);
      drawLabels(
        context,
        formattedRowLabels,
        formattedColumnLabels,
        rowScale,
        columnScale,
        highlightedRow,
        highlightedColumn,
        chartWidth,
        chartDim
      );
    },
    canvasWidth,
    canvasHeight,
    [highlightedColumn, highlightedRow]
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
          {infoText["HEATMAP"]["title"] + "    "}
          <Info name={"HEATMAP"} direction="s" />
        </Grid>
        <Grid item>
          <canvas ref={ref} />
        </Grid>
      </Grid>
    </Paper>
  );
};

const formatLabelData = (values) => {
  if (typeof values[0] === "string") {
    return values.map((value) => ({
      value,
      label: value,
      color: DEFAULT_LABEL_COLOR,
    }));
  }
  return values;
};

const drawLabels = (
  context,
  rowValues,
  columnValues,
  rowScale,
  columnScale,
  highlightedRow,
  highlightedColumn,
  chartWidth,
  chartDim
) => {
  context.font = LABEL_FONT;
  columnValues.forEach((columnData) => {
    const { value, label, color } = columnData;

    context.save();
    context.translate(
      columnScale(value) + columnScale.bandwidth() / 2,
      COLUMN_LABEL_SPACE
    );

    context.rotate((322 * Math.PI) / 180);

    context.textAlign = "left";
    context.textBaseline = "bottom";
    context.fillStyle = color;
    context.globalAlpha = isHighlighted(
      highlightedColumn,
      highlightedRow,
      value,
      undefined
    )
      ? 1
      : 0.2;

    // artifact from NDV :(
    if (label.indexOf("/") !== -1) {
      context.fillText(label.split("/")[0] + "/", 0, 0);
      context.fillText(label.split("/")[1], 10, 10);
    } else {
      context.fillText(label, 0, 0);
    }

    context.restore();
    context.fill();
  });

  rowValues.forEach((rowData) => {
    const { value, label, color } = rowData;
    context.font = "bold 12px Helvetica";
    context.fillStyle = color;
    context.globalAlpha = isHighlighted(
      highlightedColumn,
      highlightedRow,
      undefined,
      value
    )
      ? 1
      : 0.2;

    context.textAlign = "left";
    context.textBaseline = "middle";

    context.fillText(
      label,
      chartWidth + 5,
      rowScale(value) + rowScale.bandwidth() / 2
    );
  });
};

const drawHeatmap = (
  context,
  freqMap,
  rowValues,
  columnValues,
  rowScale,
  columnScale,
  highlightedRow,
  highlightedColumn,
  heatmapColor
) => {
  const cellWidth = columnScale.bandwidth();
  const cellHeight = rowScale.bandwidth();

  rowValues.forEach((rowName) => {
    const rowData = freqMap[rowName];
    const yPos = rowScale(rowName);
    const total = rowData["total"];

    columnValues.forEach((columnName) => {
      const colFreq = rowData[columnName];
      const xPos = columnScale(columnName);

      context.globalAlpha = isHighlighted(
        highlightedColumn,
        highlightedRow,
        rowName,
        columnName
      )
        ? 1
        : 0.2;

      if (colFreq) {
        context.fillStyle = heatmapColor(colFreq);
      } else {
        context.fillStyle = HEATMAP_NULL_COLOR;
      }
      context.fillRect(xPos, yPos, cellWidth, cellHeight);

      if (colFreq) {
        context.fillStyle = "black";
        context.font = CELL_FONT;
        context.textAlign = "center";
        context.textBaseline = "middle";
        context.fillText(
          `${colFreq} (${Math.round((colFreq * 100) / total)}%)`,
          xPos + cellWidth / 2,
          yPos + cellHeight / 2
        );
      }
    });
  });
};

const isHighlighted = (highlightedColumn, highlightedRow, column, row) => {
  if (highlightedColumn !== null || highlightedRow !== null) {
    if (highlightedColumn === column || highlightedRow === row) {
      return true;
    } else {
      return false;
    }
  } else {
    return true;
  }
};
export default Heatmap;
