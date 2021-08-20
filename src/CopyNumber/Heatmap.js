import React, { useState } from "react";
import * as d3 from "d3";
import Grid from "@material-ui/core/Grid";

import { useCanvas } from "../utils/useCanvas";
import { useD3 } from "../utils/useD3";
import { colorScale } from "./utils";
import Legend from "../Legend/Horizontal";

const MINIMAP_WIDTH = 150;
const PADDING = 10;

export const HEATMAP_SPACING = MINIMAP_WIDTH + PADDING;

const MIN_HEATMAP_ROW_SIZE = 7;
const MIN_MINIMAP_ROW_SIZE = 2;

const CopyNumberHeatmap = ({
  data,
  width,
  height,
  categories,
  chromosomes,
  onChange,
}) => {
  const [start, setStart] = useState(0);
  const [hoveredRow, setHoveredRow] = useState(null);
  const bpTotal = chromosomes.reduce((currSum, chr) => chr.length + currSum, 0);
  const heatmapWidth = width - MINIMAP_WIDTH - PADDING;

  const copyNumbers = colorScale.domain();
  const legendLabels = copyNumbers.map((value) => ({
    value,
    label: value === copyNumbers.length - 1 ? `≥${value}` : value,
  }));

  const numRows = Math.min(
    Math.floor(height / MIN_HEATMAP_ROW_SIZE),
    data.length
  );

  const shownData = data.slice(start, start + numRows + 1);
  const yScale = d3
    .scaleBand()
    .domain(shownData.map((row) => row["id"]))
    .range([0, height]);
  const rowHeight = yScale.bandwidth();

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawHeatmap(context, shownData, bpTotal, heatmapWidth, height);

      if (hoveredRow !== null) {
        context.globalAlpha = 1;
        context.fillStyle = "#000000";

        const hoveredY = hoveredRow * rowHeight;
        context.fillRect(0, hoveredY + rowHeight - 1, heatmapWidth, 1);
        context.fillRect(0, hoveredY, 1, rowHeight);
      }
    },
    heatmapWidth,
    height,
    [start, hoveredRow]
  );

  const indicatorRef = useD3(
    (svg) => {
      svg
        .on("mousemove", function () {
          const coordinates = d3.mouse(this);

          const rowIndex = Math.max(0, Math.floor(coordinates[1] / rowHeight));
          if (hoveredRow !== rowIndex) {
            setHoveredRow(rowIndex);
            onChange(shownData[rowIndex]);
          }
        })
        .on("mouseleave", function () {
          setHoveredRow(null);
          onChange(null);
        });
    },
    heatmapWidth,
    height,
    []
  );

  return (
    <Grid container direction="column" style={{ width }}>
      <Grid item style={{ textAlign: "right" }}>
        <Legend
          width={width}
          title={"Copy Number"}
          ticks={legendLabels}
          colorScale={colorScale}
        />
      </Grid>
      <Grid item container direction="row">
        <Grid
          item
          style={{
            position: "relative",
            height: height + 1,
          }}
        >
          <canvas ref={ref} style={{ position: "absolute" }} />
          <svg ref={indicatorRef} style={{ position: "absolute" }} />
        </Grid>
        <Grid item style={{ paddingLeft: PADDING + heatmapWidth }}>
          <Minimap
            data={data}
            width={MINIMAP_WIDTH}
            height={height}
            bpTotal={bpTotal}
            start={start}
            onChange={setStart}
            rowsHighlighted={numRows}
          />
        </Grid>
      </Grid>
    </Grid>
  );
};

const Labels = ({ data, labels, width, height }) => {};

const Minimap = ({
  data,
  width,
  height,
  bpTotal,
  start,
  onChange,
  rowsHighlighted,
}) => {
  const ratio = height / data.length;

  const nth = Math.ceil(MIN_MINIMAP_ROW_SIZE / ratio);

  const miniData =
    ratio < MIN_MINIMAP_ROW_SIZE
      ? data.filter((_, index) => index % nth === 0)
      : data;

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawHeatmap(context, miniData, bpTotal, width, height);
    },
    width,
    height,
    []
  );
  const brushSVGRef = useD3(
    (svg) => {
      const brush = d3.brushY().on("end", brushed);
      const boxSize = rowsHighlighted * ratio;
      const startingPoint = Math.floor(ratio * start);

      function brushed() {
        const selection = d3.event.selection;
        const newStart = Math.floor(selection[0] / ratio);
        if (d3.event.sourceEvent) {
          onChange(Math.max(0, newStart));
        }
      }

      const gBrush = svg.append("g");
      gBrush.call(brush);

      gBrush.selectAll(".handle").remove();
      gBrush.select(".overlay").remove();

      brush.move(gBrush, [startingPoint, startingPoint + boxSize]);
    },
    width,
    height,
    [start]
  );
  return (
    <div>
      <canvas ref={canvasRef} style={{ position: "absolute" }} />
      <svg ref={brushSVGRef} style={{ position: "absolute" }} />
    </div>
  );
};

const drawHeatmap = (context, data, bpTotal, width, height) => {
  const xScale = d3.scaleLinear().domain([0, bpTotal]).range([0, width]);
  const bandScale = d3.scaleLinear().domain([0, bpTotal]).range([0, width]);
  const yScale = d3
    .scaleBand()
    .domain(data.map((row) => row["id"]))
    .range([0, height]);
  const rowHeight = yScale.bandwidth();
  data.forEach((row) => {
    const y = yScale(row["id"]);
    const segs = row["segs"];

    segs.forEach((seg) => {
      context.fillStyle = colorScale(seg["state"]);
      context.fillRect(
        xScale(seg["genomeStart"]),
        y,
        bandScale(seg["length"]),
        rowHeight
      );
    });
  });
};

export default CopyNumberHeatmap;
