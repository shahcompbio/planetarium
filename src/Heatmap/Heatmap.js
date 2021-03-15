import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import { useDashboardState } from "../PlotState/dashboardState";
import { useCanvas } from "../components/utils/useCanvas";

const HEATMAP_COLOR = ["#ffec8b", "#d91e18"];
const HEATMAP_NULL_COLOR = "#eeeeee";

const X_PADDING = 50;

const DataWrapper = ({
  data,
  chartDim,
  selectedSubtype,
  selectedClonotype,
}) => {
  const [
    {
      clonotypeParam,
      sampleTen,
      topTenNumbering,
      colors,
      subtypeParam,
      fontSize,
    },
  ] = useDashboardState();

  const topTenData = data.filter((record) =>
    Object.keys(sampleTen).includes(record[clonotypeParam])
  );

  console.log(sampleTen);

  // proper should be
  // data, chartDim, column, row, highlightedhighlightedColumn, highlightedhighlightedRow

  return (
    <Heatmap
      data={data}
      chartDim={chartDim}
      column={subtypeParam}
      row={clonotypeParam}
      highlightedColumn={selectedSubtype}
      highlightedRow={selectedClonotype}
      rowLabels={Object.keys(sampleTen).sort(([, a], [, b]) => b - a)}
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
  const [{ topTenNumbering, colors, fontSize }] = useDashboardState();

  const columnValues =
    columnLabels || _.uniq(data.map((record) => record[column])).sort();
  const rowValues =
    rowLabels || _.uniq(data.map((record) => record[row])).sort();
  const xAxis = d3
    .scaleBand()
    .domain(columnValues)
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"] - X_PADDING]);
  const heatmapWidth =
    (chartDim["chart"]["x2"] - chartDim["chart"]["x1"] - 50) /
    columnValues.length;

  const subTypes = data.reduce((final, current) => {
    var allSamples = final;
    const subtype = current[column];
    allSamples[subtype] = final.hasOwnProperty(subtype)
      ? allSamples[subtype] + 1
      : 1;
    final = allSamples;
    return final;
  }, {});

  // var largestFreq = 0;
  const subtypeStats = data.reduce(
    (final, current) => {
      const seq = current[row];

      if (rowValues.includes(seq)) {
        const subtype = current[column];

        if (final[seq].hasOwnProperty(subtype)) {
          final[seq][subtype] = final[seq][subtype] + 1;
        } else {
          final[seq][subtype] = 1;
        }
      }
      // largestFreq =
      //   final[seq][subtype] > largestFreq ? final[seq][subtype] : largestFreq;

      return final;
    },
    {
      // this initializes an empty map
      ...rowValues.reduce((final, entry) => {
        final[entry] = {};
        return final;
      }, {}),
    }
  );

  const alphaIndexing = Object.entries(topTenNumbering)
    .map((entry) => entry[1])
    .sort((a, b) => a - b);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawLabels(context, highlightedRow);
      drawHeatmap(context, chartDim, data, highlightedColumn, highlightedRow);
    },
    chartDim["width"],
    chartDim["height"],
    [highlightedColumn, highlightedRow]
  );

  function drawLabels(context, highlightedRow) {
    context.beginPath();
    context.globalAlpha = 1;
    context.fillStyle = "#000000";
    context.font = "normal " + fontSize.axisLabelFontSize + "px Helvetica";

    columnValues.forEach(function(d, index) {
      context.save();
      context.translate(
        xAxis(d) + heatmapWidth / 2,
        chartDim["chart"]["y1"] - 5
      );
      context.rotate((322 * Math.PI) / 180);
      if (d.indexOf("/") !== -1) {
        context.fillText(d.split("/")[0] + "/", 5, 0);
        context.fillText(d.split("/")[1], 5, 10);
      } else {
        context.fillText(d, 5, 5);
      }
      context.restore();
      context.fill();
    });
    const startingY = chartDim["chart"]["y1"];

    const sequenceLength = Object.keys(subtypeStats).length;
    const heatmaphighlightedRowSpace = 3;
    const heatmapHeight =
      (chartDim["chart"]["y2"] -
        chartDim["chart"]["y1"] -
        sequenceLength * heatmaphighlightedRowSpace) /
      sequenceLength;

    Object.keys(subtypeStats)
      .sort((a, b) => {
        return (
          alphaIndexing.indexOf(topTenNumbering[a]) -
          alphaIndexing.indexOf(topTenNumbering[b])
        );
      })
      .map((sequence, index) => {
        const yPos =
          startingY +
          heatmapHeight * index +
          heatmaphighlightedRowSpace * index;
        context.fillStyle = colors(sequence);
        context.globalAlpha =
          highlightedRow !== null && highlightedRow !== undefined
            ? highlightedRow === sequence
              ? 1
              : 0.2
            : 1;

        context.font = "bold " + fontSize.axisLabelFontSize + "px Helvetica";

        context.fillText(
          topTenNumbering[sequence] + " - " + sequence,
          chartDim["chart"]["x2"] - 50,
          yPos + (3 * heatmapHeight) / 4
        );
      });
  }

  function drawHeatmap(
    context,
    allDim,
    data,
    highlightedColumn,
    highlightedRow
  ) {
    const dimensions = allDim;
    const columnValues = Object.keys(subTypes);
    const freqColouring = d3
      .scaleLinear()
      .range(HEATMAP_COLOR)
      .domain([0, 96]);

    const sequenceLength = Object.keys(subtypeStats).length;
    const heatmaphighlightedRowSpace = 3;
    const heatmapHeight =
      (dimensions["chart"]["y2"] -
        dimensions["chart"]["y1"] -
        sequenceLength * heatmaphighlightedRowSpace) /
      sequenceLength;
    const startingY = dimensions["chart"]["y1"];

    Object.keys(subtypeStats)
      .sort((a, b) => {
        return (
          alphaIndexing.indexOf(topTenNumbering[a]) -
          alphaIndexing.indexOf(topTenNumbering[b])
        );
      })
      .map((sequence, index) => {
        const yPos =
          startingY +
          heatmapHeight * index +
          heatmaphighlightedRowSpace * index;

        context.font = "normal " + fontSize.axisLabelFontSize + "px Helvetica";

        const seqSubtypes = subtypeStats[sequence];
        columnValues.map((subtype) => {
          context.globalAlpha =
            (highlightedColumn !== null || highlightedRow !== null) &&
            (highlightedColumn !== undefined || highlightedRow !== undefined)
              ? highlightedColumn === subtype || highlightedRow === sequence
                ? 1
                : 0.2
              : 1;

          if (seqSubtypes.hasOwnProperty(subtype)) {
            context.fillStyle = freqColouring(subtypeStats[sequence][subtype]);
          } else {
            context.fillStyle = HEATMAP_NULL_COLOR;
          }

          context.fillRect(
            xAxis(subtype),
            yPos,
            heatmapWidth - 3,
            heatmapHeight
          );

          if (seqSubtypes.hasOwnProperty(subtype)) {
            context.fillStyle = "black";
            context.globalAlpha =
              highlightedColumn !== null && highlightedColumn !== undefined
                ? highlightedColumn === subtype
                  ? 1
                  : 0.2
                : 1;
            const freq = subtypeStats[sequence][subtype];
            const freqTextX =
              freq > 9
                ? xAxis(subtype) + (heatmapWidth - 15) / 2
                : xAxis(subtype) + (heatmapWidth - 6) / 2;
            context.fillText(freq, freqTextX, yPos + heatmapHeight / 2 + 5);
          }
        });
      });
  }
  return (
    <div>
      <div
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative",
        }}
      >
        <div
          id="heatmap"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex",
          }}
        >
          <canvas ref={ref} />
        </div>
      </div>
    </div>
  );
};
export default DataWrapper;
