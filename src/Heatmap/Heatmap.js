import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useDashboardState } from "../PlotState/dashboardState";
import { useCanvas } from "../components/utils/useCanvas";

const HEATMAP_NULL_COLOR = "#eeeeee";
const HEATMAP_COLOR_SCALE = d3
  .scaleLinear()
  .range(["#ffec8b", "#d91e18"])
  .domain([0, 1]);
const CELL_FONT = "normal 12px Helvetica";

const COLUMN_LABEL_SPACE = 100;
const ROW_LABEL_SPACE = 300;
const DEFAULT_LABEL_COLOR = "#000000";
const LABEL_FONT = "bold 12px Helvetica";

const DataWrapper = ({
  data,
  chartDim,
  selectedSubtype,
  selectedClonotype,
}) => {
  const [
    { clonotypeParam, sampleTen, colors, subtypeParam },
  ] = useDashboardState();

  const rowLabels = Object.keys(sampleTen)
    .sort(([, a], [, b]) => b - a)
    .map((clonotype, index) => ({
      value: clonotype,
      label: `SEQ${index + 1} - ${clonotype}`,
      color: colors(clonotype),
    }));

  return (
    <Heatmap
      data={data}
      chartDim={chartDim}
      column={subtypeParam}
      row={clonotypeParam}
      highlightedColumn={selectedSubtype}
      highlightedRow={selectedClonotype}
      rowLabels={rowLabels}
    />
  );
};

const Heatmap = ({
  data,
  chartDim,
  column,
  row,
  highlightedColumn,
  highlightedRow,
  columnLabels,
  rowLabels,
}) => {
  const columnValues =
    columnLabels || _.uniq(data.map((record) => record[column])).sort();
  const rowValues =
    rowLabels.map((row) => row["value"]) ||
    _.uniq(data.map((record) => record[row])).sort();
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

  const columnMap = rowValues.reduce(
    (currMap, rowName) => ({ ...currMap, [rowName]: 0 }),
    { total: 0 }
  );
  const initMap = columnValues.reduce(
    (currMap, columnName) => ({ ...currMap, [columnName]: columnMap }),
    {}
  );

  const freqMap = data.reduce((currMap, record) => {
    const columnValue = record[column];
    const rowValue = record[row];

    if (currMap.hasOwnProperty(columnValue)) {
      const currColumn = currMap[columnValue];
      if (currColumn.hasOwnProperty(rowValue)) {
        return {
          ...currMap,
          [columnValue]: {
            ...currColumn,
            [rowValue]: currColumn[rowValue] + 1,
            total: currColumn["total"] + 1,
          },
        };
      } else {
        return {
          ...currMap,
          [columnValue]: { ...currColumn, total: currColumn["total"] + 1 },
        };
      }
    }

    return currMap;
  }, initMap);

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
        highlightedColumn
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
      id="heatmap"
      style={{
        margin: 30,
        pointerEvents: "all",
        display: "flex",
      }}
    >
      <canvas ref={ref} />
    </div>
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
  chartWidth
) => {
  context.font = LABEL_FONT;

  columnValues.map((columnData) => {
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

  rowValues.map((rowData) => {
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
  highlightedColumn
) => {
  const cellWidth = columnScale.bandwidth();
  const cellHeight = rowScale.bandwidth();

  columnValues.map((columnName) => {
    const columnData = freqMap[columnName];
    const xPos = columnScale(columnName);
    const total = columnData["total"];

    rowValues.map((rowName) => {
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
        context.fillStyle = HEATMAP_COLOR_SCALE(rowFreq / total);
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
export default DataWrapper;
