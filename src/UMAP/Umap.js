import React, { useEffect, useState } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";
import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis, drawMiniAxis } from "../DrawingUtils/utils.js";

import _ from "lodash";

import { useCanvas } from "../components/utils/useCanvas";
import { useD3 } from "../components/utils/useD3";

const radiusMax = 20;

const PADDING = 10;
const TITLE_HEIGHT = 30;
const AXIS_SPACE = 20;
const LEGEND_WIDTH = 180;

const LINE_GRAPH_SPACE = 50;

const NULL_POINT_COLOR = "#d2d7d3";
const UNHIGHLIGHTED_POINT_COLOR = "grey";
const POINT_RADIUS = 2;

const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";

const DataWrapper = ({
  chartName,
  data,
  chartDim,
  clonotypeLabels,
  selectedClonotype,
  hoveredClonotype,
  setSelectedClonotype,
}) => {
  const [{ xParam, yParam, clonotypeParam }] = useDashboardState();

  const setHighlighted = (event, value) => {
    if (event === "mouseenter") {
      setSelectedClonotype({ hover: value });
    } else if (event === "mousedown") {
      setSelectedClonotype({
        hover: null,
        selected: value,
      });
    } else if (event === "mouseout") {
      setSelectedClonotype({
        hover: null,
      });
    }
  };
  return (
    <UMAP
      chartDim={chartDim}
      chartName={chartName}
      data={data}
      highlighted={hoveredClonotype || selectedClonotype}
      xParam={xParam}
      yParam={yParam}
      subsetParam={clonotypeParam}
      subsetLabels={clonotypeLabels}
      setHighlighted={setHighlighted}
    />
  );
};

const UMAP = ({
  chartDim,
  chartName,
  data,
  highlighted,
  xParam,
  yParam,
  subsetParam,
  subsetLabels,
  setHighlighted,
}) => {
  const canvasWidth = chartDim["width"] - LEGEND_WIDTH - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;

  const chartWidth =
    canvasWidth - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;
  const chartHeight =
    canvasHeight - AXIS_SPACE - PADDING - PADDING - LINE_GRAPH_SPACE;

  const chartX = PADDING;
  const chartY = PADDING + LINE_GRAPH_SPACE;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([chartX, chartX + chartWidth]);
  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([chartY, chartY + chartHeight]);

  const subsetValues = subsetLabels.map(({ value }) => value);
  const subsetColors = subsetLabels.map(({ color }) => color);

  const colorScale = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(subsetColors);

  // filter data by top 10 (should probably be done in data wrapper tbh)
  // bin the data

  // for each clonotype
  // bin by x, y

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawUMAPAxis(
        context,
        canvasHeight - AXIS_SPACE - PADDING,
        xParam,
        yParam
      );
      drawPoints(
        context,
        data,
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        highlighted,
        subsetValues,
        colorScale
      );
    },
    canvasWidth,
    canvasHeight,
    [highlighted]
  );
  const svgRef = useD3((svg) => {}, LEGEND_WIDTH, chartHeight, [highlighted]);

  return (
    <Paper
      style={{
        margin: 10,
        padding: "10px 0px",
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
            paddingRight: PADDING,
            marginBottom: 10,
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
        <Grid container direction="row" style={{ padding: 0 }}>
          <Grid item>
            <canvas ref={canvasRef} />
          </Grid>
          <Grid item>
            <svg ref={svgRef} />
          </Grid>
        </Grid>
      </Grid>
    </Paper>
  );
};

const drawPoints = (
  context,
  data,
  xScale,
  yScale,
  xParam,
  yParam,
  subsetParam,
  highlighted,
  subsetLabels,
  colorScale
) => {
  context.lineWidth = 1;
  context.globalAlpha = 1;

  const [subsetData, backgroundData] = _.partition(
    data,
    (datum) => subsetLabels.indexOf(datum[subsetParam]) !== -1
  );

  context.fillStyle = NULL_POINT_COLOR;
  backgroundData.forEach((point) => {
    context.beginPath();
    context.arc(
      xScale(point[xParam]),
      yScale(point[yParam]),
      1.5,
      0,
      Math.PI * 2,
      true
    );
    context.fill();
  });

  subsetData.forEach((point) => {
    context.fillStyle = isHighlighted(point[subsetParam], highlighted)
      ? colorScale(point[subsetParam])
      : UNHIGHLIGHTED_POINT_COLOR;

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

const drawUMAPAxis = (context, startY, xParam, yParam) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  const START_X = AXIS_SPACE / 2;
  const START_Y = startY + AXIS_SPACE / 2;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
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

const isHighlighted = (datumValue, highlighted) =>
  highlighted === null || datumValue === highlighted;

export const clearAll = (context, chartDim) =>
  context.clearRect(0, 0, chartDim["chart"].x2 + 50, chartDim["chart"].y2 + 80);
export function drawPoint(
  context,
  point,
  fill,
  isCountInsignificant,
  x,
  y,
  xParam,
  yParam
) {
  const radius = isCountInsignificant ? 2 : point["radius"];

  context.beginPath();
  context.arc(x(point[xParam]), y(point[yParam]), radius, 0, Math.PI * 2, true);

  context.fillStyle = fill;
  context.fill();
}

const Umap2 = ({
  chartName,
  data,
  chartDim,
  selectedClonotype,
  hoveredClonotype,
  setSelectedClonotype,
}) => {
  const [
    {
      xParam,
      yParam,
      cellIdParam,
      clonotypeParam,
      fontSize,
      topTen,
      colors,
      topTenNumbering,
    },
  ] = useDashboardState();
  const [radiusAdjust, setRadius] = useState(10);
  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = data.filter((row) =>
    sampleTen.hasOwnProperty(row[clonotypeParam])
  );

  const dim = chartDim["chart"];
  // X axis
  var x = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([dim.x1, dim.x2]);

  // Y axis
  var y = d3
    .scaleLinear()
    .domain([yMin, yMax])
    .range([dim.y2, dim.y1]);

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      const selection =
        hoveredClonotype !== null
          ? hoveredClonotype
          : selectedClonotype !== null
          ? selectedClonotype
          : null;
      reDraw(
        context,
        x,
        y,
        dim,
        sampleData,
        sampleTen,
        colors,
        topTenNumbering,
        selection
      );
    },
    chartDim.width,
    chartDim.height,
    [selectedClonotype, hoveredClonotype, radiusAdjust]
  );
  const svgRef = useD3(
    (svg) => {
      drawLegend(svg);
    },
    250,
    250,
    []
  );

  // useEffect(() => {
  //   if (data.length > 0 && colors) {
  //     init(data, chartDim);
  //   }
  // }, [colors]);

  // useEffect(() => {
  //   if (context) {
  //     clearAll(context, chartDim);
  //     reDraw(
  //       context,
  //       x,
  //       y,
  //       dim,
  //       sampleData,
  //       sampleTen,
  //       colors,
  //       topTenNumbering
  //     );
  //   }
  // }, [radiusAdjust]);

  // useEffect(() => {
  //   if (context) {
  //     const selection =
  //       hoveredClonotype !== null
  //         ? hoveredClonotype
  //         : selectedClonotype !== null
  //         ? selectedClonotype
  //         : null;

  //     if (selection !== null) {
  //       clearAll(context, chartDim);
  //       context.beginPath();
  //       reDraw(
  //         context,
  //         x,
  //         y,
  //         dim,
  //         sampleData,
  //         sampleTen,
  //         colors,
  //         topTenNumbering,
  //         selection
  //       );
  //     } else {
  //       clearAll(context, chartDim);
  //       context.beginPath();
  //       reDraw(
  //         context,
  //         x,
  //         y,
  //         dim,
  //         sampleData,
  //         sampleTen,
  //         colors,
  //         topTenNumbering
  //       );
  //     }
  //   }
  // }, [selectedClonotype, hoveredClonotype, context]);
  function drawOutline(context, x, y, data, colors) {
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    context.globalAlpha = 0.5;
    data.forEach((point) => {
      context.beginPath();
      context.arc(
        x(point[xParam]),
        y(point[yParam]),
        1.5,
        0,
        Math.PI * 2,
        true
      );
      context.fillStyle = "#d2d7d3";

      context.fill();
    });
    context.globalAlpha = 1;
  }

  function drawLineGraph(
    context,
    x,
    y,
    dimensions,
    data,
    topTen,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    const maxValue = Math.max(...Object.entries(topTen).map((row) => row[1]));
    const minValue = Math.min(...Object.entries(topTen).map((row) => row[1]));

    const lineXFreq = d3
      .scaleLinear()
      .domain([minValue, maxValue])
      .range([dimensions.x2 + 30, dimensions.x2 + 65]);

    const lineYFreq = d3
      .scaleLinear()
      .domain([maxValue, minValue])
      .range([dimensions.y1 - 65, dimensions.y1 - 30]);

    const lineXaxis = d3
      .line()
      .x(function(d) {
        return x((d.x0 + d.x1) / 2);
      })
      .y(function(d) {
        const freq = Object.entries(d).length;
        return lineYFreq(freq);
      })
      .curve(d3.curveCatmullRom.alpha(0.5))
      .context(context);

    const lineYaxis = d3
      .line()
      .y(function(d) {
        return y((d.x0 + d.x1) / 2);
      })
      .x(function(d) {
        const freq = Object.entries(d).length;
        return lineXFreq(freq);
      })
      .curve(d3.curveCatmullRom.alpha(0.5))
      .context(context);

    var nestedSamples = Array.from(
      d3Array.group(data, (d) => d[clonotypeParam]),
      ([key, value]) => ({ key, value })
    );

    //if selected, move to end so it's drawn last
    if (selectedClonotype) {
      const keys = nestedSamples.map((row) => row["key"]);

      nestedSamples = [
        ...nestedSamples.filter((row) => row["key"] !== selectedClonotype),
        nestedSamples[keys.indexOf(selectedClonotype)],
      ];
    }
    nestedSamples.reduce((final, clonotype) => {
      const xBins = d3Array
        .bin()
        .value((d) => d[xParam])
        .domain(x.domain())
        .thresholds(x.ticks(10))(clonotype["value"]);

      context.beginPath();
      lineXaxis(xBins);

      const isSelected =
        selectedClonotype && selectedClonotype === clonotype["key"];

      context.lineWidth = isSelected ? 2.5 : 1.5;

      context.strokeStyle = selectedClonotype
        ? isSelected
          ? colors(clonotype["key"])
          : "#e8e8e8"
        : colors(clonotype["key"]);

      context.stroke();

      const yBins = d3Array
        .bin()
        .value((d) => d[yParam])
        .domain(y.domain())
        .thresholds(y.ticks(10))(clonotype["value"]);

      context.beginPath();
      lineYaxis(yBins);

      context.lineWidth = isSelected ? 2.5 : 1.5;

      context.strokeStyle = selectedClonotype
        ? isSelected
          ? colors(clonotype["key"])
          : "#e8e8e8"
        : colors(clonotype["key"]);

      context.stroke();
      final[clonotype["key"]] = { x: xBins, y: yBins };
      return final;
    }, []);
  }

  function drawPoints(
    context,
    x,
    y,
    dimensions,
    data,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    var nestedSamples = Array.from(
      d3Array.group(data, (d) => d[clonotypeParam]),
      ([key, value]) => ({ key, value })
    );
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    // for each clonotype
    const merge = nestedSamples.map((clonotype, i) => {
      // bin the data across the axis
      const xBins = d3Array
        .bin()
        .value((d) => d[xParam])
        .domain(x.domain())
        .thresholds(x.ticks(8))(clonotype["value"]);

      const yBins = d3Array
        .bin()
        .value((d) => d[yParam])
        .domain(y.domain())
        .thresholds(y.ticks(8))(clonotype["value"]);

      // for each bin,
      const xBinned = xBins.reduce((final, xBin) => {
        // this just checks to see if there are items in the bin
        const freq = Object.entries(xBin).length - 2;

        if (freq > 0) {
          const rows = Object.entries(xBin).filter((row, index) => {
            if (row[0] !== "x1" && row[1] !== "x0") {
              return true;
            } else {
              return false;
            }
          });

          // what?
          // ii think it's a map and reduce. it adds the freq to each record, and then adds it to a map by cell_id
          final = {
            ...final,
            ...rows.reduce((finalRow, row) => {
              finalRow[row[1][cellIdParam]] = { ...row[1], xRadius: freq };
              return finalRow;
            }, {}),
          };
        }
        return final;
      }, {});

      const yBinned = yBins.reduce((final, yBin) => {
        const freq = Object.entries(yBin).length - 2;

        if (freq > 0) {
          const rows = Object.entries(yBin).filter((row, index) => {
            if (row[0] !== "x1" && row[1] !== "x0") {
              return true;
            } else {
              return false;
            }
          });

          final = {
            ...final,
            ...rows.reduce((finalRow, row) => {
              finalRow[row[1][cellIdParam]] = { ...row[1], yRadius: freq };
              return finalRow;
            }, {}),
          };
        }
        return final;
      }, {});

      const merged = Object.keys(yBinned).map((item, i) =>
        Object.assign({}, yBinned[item], xBinned[item])
      );
      return merged;
    });

    const sortedmerge = merge
      .flat(1)
      .filter((point) => point.hasOwnProperty(cellIdParam))
      .map((point) => ({
        ...point,
        radius: (point["xRadius"] + point["yRadius"]) / radiusAdjust,
      }))
      .sort((a, b) => b.radius - a.radius);

    const isCountInsignificant =
      Math.max(...sortedmerge.map((point) => point["radius"])) < 1
        ? true
        : radiusAdjust == radiusMax
        ? true
        : false;

    sortedmerge.map((point) => {
      const fill = selectedClonotype ? "grey" : colors(point[clonotypeParam]);
      context.globalAlpha = selectedClonotype ? 0.5 : 1;
      drawPoint(
        context,
        point,
        fill,
        isCountInsignificant,
        x,
        y,
        xParam,
        yParam
      );
    });

    context.globalAlpha = 1;
    //if selected, move to end so it's drawn last
    if (selectedClonotype) {
      sortedmerge
        .filter((row) => row[clonotypeParam] === selectedClonotype)
        .map((point) => {
          const fill = colors(point[clonotypeParam]);
          drawPoint(
            context,
            point,
            fill,
            isCountInsignificant,
            x,
            y,
            xParam,
            yParam
          );
        });
    }
  }

  function reDraw(
    context,
    x,
    y,
    dim,
    sampleData,
    sampleTen,
    colors,
    topTenNumbering,
    selectedClonotype
  ) {
    drawMiniAxis(context, x, y, dim, xParam, yParam);
    //  drawAxisLabels(context, x, y, chartDim);
    drawOutline(context, x, y, data, colors);
    drawLineGraph(
      context,
      x,
      y,
      dim,
      sampleData,
      sampleTen,
      colors,
      topTenNumbering,
      selectedClonotype
    );
    drawPoints(
      context,
      x,
      y,
      dim,
      sampleData,
      colors,
      topTenNumbering,
      selectedClonotype
    );
  }
  function drawLegend(svg) {
    const mouseInteractions = (element) =>
      element
        .on("mouseenter", function(d) {
          setSelectedClonotype({
            hover: d[0],
          });
        })
        .on("mousedown", function(d, i) {
          d3.event.stopPropagation();
          setSelectedClonotype({
            hover: null,
            selected: d[0],
          });
        })
        .on("mouseout", function(event, d) {
          setSelectedClonotype({
            hover: null,
          });
        });

    const subsets = svg
      .selectAll("g")
      .data(topTen)
      .enter()
      .append("g")
      .attr("cursor", "pointer");

    subsets
      .append("rect")
      .attr("width", fontSize.legendSquare)
      .attr("height", fontSize.legendSquare)
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 5;
      })
      .attr("y", function(d, i) {
        return i * 20 + chartDim["legend"].y1 - 5;
      })
      .attr("fill", function(d) {
        return colors(d[0]);
      });

    subsets
      .append("text")
      .attr("x", function(d) {
        return chartDim["legend"].x1 + 20;
      })
      .attr("y", function(d, i) {
        return i * 20 + chartDim["legend"].y1;
      })
      .attr("dy", ".35em")
      .text(function(d) {
        return topTenNumbering[d[0]] + " - " + d[0] + " - " + d[1];
      })
      .attr("font-weight", "700")
      .attr("font-size", fontSize.legendFontSize + "px")
      .attr("fill", function(d) {
        return colors(d[0]);
      });

    subsets.call(mouseInteractions);

    // legend
    //   .selectAll("text")
    //   .on("mouseenter", null)
    //   .on("mousedown", null)
    //   .on("mouseout", null);
    // legend.selectAll("*").remove();

    // const legendRect = legend.selectAll("rect").data(topTen);

    // const legendRectEnter = legendRect
    //   .append("rect")
    //   .attr("width", fontSize.legendSquare)
    //   .attr("height", fontSize.legendSquare)
    //   .attr("x", function(d) {
    //     return chartDim["legend"].x1 + 5;
    //   })
    //   .attr("y", function(d, i) {
    //     return i * 20 + chartDim["legend"].y1 - 5;
    //   })
    //   .attr("fill", function(d) {
    //     return colors(d[0]);
    //   });

    // legendRect
    //   .enter()
    //   .append("rect")
    //   .attr("width", fontSize.legendSquare)
    //   .attr("height", fontSize.legendSquare)
    //   .attr("x", function(d) {
    //     return chartDim["legend"].x1 + 5;
    //   })
    //   .attr("y", function(d, i) {
    //     return i * 20 + chartDim["legend"].y1 - 5;
    //   })
    //   .attr("fill", function(d) {
    //     return colors(d[0]);
    //   })
    //   .call(mouseInteractions);

    // const legendText = legend.selectAll("text").data(topTen);

    // legendText
    //   .append("text")
    //   .attr("x", function(d) {
    //     return chartDim["legend"].x1 + 20;
    //   })
    //   .attr("y", function(d, i) {
    //     return i * 20 + chartDim["legend"].y1;
    //   })
    //   .attr("dy", ".35em")
    //   .text(function(d) {
    //     return topTenNumbering[d[0]] + " - " + d[0] + " - " + d[1];
    //   })
    //   .attr("font-weight", "700")
    //   .attr("font-size", fontSize.legendFontSize + "px")
    //   .attr("fill", function(d) {
    //     return colors(d[0]);
    //   })
    //   .attr("cursor", "pointer")
    //   .on("mouseenter", null)
    //   .on("mousedown", null)
    //   .on("mouseout", null);

    // legendText
    //   .enter()
    //   .append("text")
    //   .attr("x", function(d) {
    //     return chartDim["legend"].x1 + 20;
    //   })
    //   .attr("y", function(d, i) {
    //     return i * 20 + chartDim["legend"].y1;
    //   })
    //   .attr("dy", ".35em")
    //   .text(function(d) {
    //     return topTenNumbering[d[0]] + " - " + d[0] + " - " + d[1];
    //   })
    //   .attr("font-weight", "700")
    //   .attr("font-size", fontSize.legendFontSize + "px")
    //   .attr("fill", function(d) {
    //     return colors(d[0]);
    //   })
    //   .attr("cursor", "pointer")
    //   .call(mouseInteractions);
  }
  // function init(data, chartDim, selectedClonotype) {
  //   var canvas = d3.select("#umapCanvas");
  //   var currContext = canvasInit(canvas, chartDim.width, chartDim.height);

  //   currContext.fillStyle = "white";
  //   currContext.fillRect(0, 0, chartDim.width, chartDim.height);

  //   saveContext(currContext);

  //   reDraw(
  //     currContext,
  //     x,
  //     y,
  //     dim,
  //     sampleData,
  //     sampleTen,
  //     colors,
  //     topTenNumbering
  //   );
  // }

  return (
    <Paper style={{ margin: 10 }}>
      <Grid
        container
        direction="row"
        justify="flex-start"
        alignItems="flex-start"
        style={{
          width: chartDim["width"] + 250,
          height: chartDim["height"],
          position: "relative",
        }}
      >
        <Grid
          item
          xs={18}
          sm={9}
          id="scatterplot"
          style={{
            pointerEvents: "all",
            display: "flex",
            paddingRight: 0,
          }}
        >
          <canvas ref={canvasRef} id="umapCanvas" />
        </Grid>
        <Grid
          item
          xs={6}
          sm={3}
          style={{ paddingLeft: 0 }}
          container
          direction="column"
          justify="flex-start"
          alignItems="flex-start"
        >
          <Grid
            item
            style={{
              marginTop: chartDim["chart"]["x1"],
              width: "100%",
              height: 80,
              paddingTop: 40,
              marginLeft: -38,
              textAlign: "left",
            }}
          >
            {infoText[chartName]["title"] + "    "}

            <Info name={chartName} direction="s" />
          </Grid>
          <Grid item style={{ marginLeft: -50, height: 250 }}>
            <svg ref={svgRef} id="umapLegend" />
          </Grid>
          <Grid item style={{ marginLeft: -38 }}>
            <label
              style={{ fontSize: 12, marginTop: -20 }}
              for="customRange2"
              class="form-label"
            >
              Radius Adjustment
            </label>
          </Grid>
          <Grid style={{ marginLeft: -38 }}>
            <input
              type="range"
              min="4"
              max={radiusMax}
              step="0.5"
              value={radiusAdjust}
              onChange={(event) => {
                setRadius(event.target.value);
              }}
              style={{ direction: "rtl" }}
              class="form-range"
              id="customRange2"
              disabled={selectedClonotype !== null}
            />
          </Grid>
        </Grid>
      </Grid>
    </Paper>
  );
};

export default Umap2;
