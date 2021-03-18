import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useCanvas } from "../utils/useCanvas";

import Info from "../../Info/Info";
import infoText from "../../Info/InfoText";

const HEATMAP_NULL_COLOR = "#eeeeee";
const HEATMAP_COLOR = ["#ffec8b", "#d91e18"];
const CELL_FONT = "normal 12px Helvetica";

const COLUMN_LABEL_SPACE = 100;
const ROW_LABEL_SPACE = 300;
const DEFAULT_LABEL_COLOR = "#000000";
const LABEL_FONT = "bold 12px Helvetica";

/*

This heatmap calculates it by total column value (or given total column value)

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
  columnTotal
}) => {
  const columnValues =
    columnLabels || _.uniq(data.map(record => record[column])).sort();
  const rowValues =
    rowLabels.map(row => row["value"]) ||
    _.uniq(data.map(record => record[row])).sort();
  const chartWidth = chartDim["width"] - ROW_LABEL_SPACE;

  const columnScale = d3
    .scaleBand()
    .domain(columnValues)
    .range([0, chartWidth])
    .paddingInner(0.03);

  const rowScale = d3
    .scaleBand()
    .domain(rowValues)
    .range([COLUMN_LABEL_SPACE, chartDim["height"]])
    .paddingInner(0.03);

  const groupedColumn = _.groupBy(data, column);
  const freqMap = Object.keys(groupedColumn).reduce(
    (currMap, columnName) => ({
      ...currMap,
      [columnName]: {
        ..._.countBy(groupedColumn[columnName], row),
        total: columnTotal
          ? columnTotal[columnName]
          : groupedColumn[columnName].length
      }
    }),
    {}
  );

  const mostFreqCount = _.reduce(
    freqMap,
    (currMax, columnData, key) => {
      return Math.max(
        currMax,
        _.reduce(
          columnData,
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
    canvas => {
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
        chartWidth
      );
    },
    chartDim["width"],
    chartDim["height"],
    [highlightedColumn, highlightedRow]
  );

  return (
    <div
      class="card"
      style={{ margin: 10, paddingBottom: 10, paddingLeft: 10 }}
    >
      <div
        class="container"
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative"
        }}
      >
        <div class="row">
          <div class="col-10">
            <div
              id="heatmap"
              style={{
                position: "absolute",
                pointerEvents: "all",
                display: "flex"
              }}
            >
              <canvas ref={ref} />
            </div>
          </div>
          <div class="col-2">
            <div
              class="card-title"
              style={{
                width: "100%",
                height: 80,
                paddingTop: 40,
                paddingLeft: -50,
                textAlign: "left"
              }}
            >
              <h6 class="card-title">
                {infoText["HEATMAP"]["title"] + "    "}

                <Info name={"HEATMAP"} direction="s" />
              </h6>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

const formatLabelData = values => {
  if (typeof values[0] === "string") {
    return values.map(value => ({
      value,
      label: value,
      color: DEFAULT_LABEL_COLOR
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
  chartWidth
) => {
  context.font = LABEL_FONT;

  columnValues.forEach(columnData => {
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

  rowValues.forEach(rowData => {
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

  columnValues.forEach(columnName => {
    const columnData = freqMap[columnName];
    const xPos = columnScale(columnName);
    const total = columnData["total"];

    rowValues.forEach(rowName => {
      const rowFreq = columnData[rowName];
      const yPos = rowScale(rowName);

      context.globalAlpha = isHighlighted(
        highlightedColumn,
        highlightedRow,
        columnName,
        rowName
      )
        ? 1
        : 0.2;

      if (rowFreq) {
        context.fillStyle = heatmapColor(rowFreq);
      } else {
        context.fillStyle = HEATMAP_NULL_COLOR;
      }
      context.fillRect(xPos, yPos, cellWidth, cellHeight);

      if (rowFreq) {
        context.fillStyle = "black";
        context.font = CELL_FONT;
        context.textAlign = "center";
        context.textBaseline = "middle";
        context.fillText(
          `${rowFreq} (${Math.round((rowFreq * 100) / total)}%)`,
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
