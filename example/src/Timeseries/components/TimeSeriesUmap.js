import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import { quantileSorted } from "d3-array";
import _ from "lodash";
import infoText from "../InfoText.js";
import { InfoBar, useCanvas, useD3 } from "@shahlab/planetarium";
import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import Card from "@material-ui/core/Card";
import CardActions from "@material-ui/core/CardActions";
import CardContent from "@material-ui/core/CardContent";
import Button from "@material-ui/core/Button";
import ButtonGroup from "@material-ui/core/ButtonGroup";
import Typography from "@material-ui/core/Typography";
import Summary from "./Summary.js";

import { makeStyles } from "@material-ui/core/styles";
import { CONSTANTS } from "../config";

const PADDING = 10;
const TITLE_HEIGHT = 30;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const AXIS_COLOR = "#000000";

const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];

const NULL_POINT_COLOR = "#e8e8e8";
const POINT_RADIUS = 2;
const PERCENTILE_RANGE = [0.25, 0.75];
const LEGEND_SQUARE_LENGTH = 10;
const LEGEND_SQUARE_SPACING = 8;
const useStyles = makeStyles({
  root: {
    minWidth: 275,
    margin: "10px 10px 0px 10px",
    padding: "0px 0px",
    float: "right",
  },
  button: { float: "right" },
  buttonGroup: { float: "left", marginRight: 15 },
});

const DataWrapper = ({
  chartName,
  data,
  genes,
  chartDim,
  selected,
  hovered,
  setSelected,
}) => {
  const { xParam, yParam, subtypeParam, timepoint } = CONSTANTS;

  const setHighlighted = (event, value, type) => {
    if (event === "mouseenter") {
      setSelected({
        hover: data
          .filter((row) => isHighlighted(row[type], value))
          .reduce((final, curr) => {
            final[curr["cell_id"]] = true;
            return final;
          }, {}),
      });
    } else if (event === "mousedown") {
      setSelected({
        hover: null,
        selected: data
          .filter((row) => isHighlighted(row[type], value))
          .reduce((final, curr) => {
            final[curr["cell_id"]] = true;
            return final;
          }, {}),
      });
    } else if (event === "mouseout") {
      setSelected({
        hover: null,
        selected: null,
      });
    }
  };

  const groupedTimepoints = _.groupBy(data, (d) => d[timepoint]);

  const [currentTimepoint, setTimepoint] = useState(
    Object.keys(groupedTimepoints)[0]
  );

  const modifiedData = groupedTimepoints[currentTimepoint];

  return (
    <UMAP
      chartDim={chartDim}
      chartName={chartName}
      data={modifiedData}
      genes={genes}
      highlighted={hovered || selected}
      xParam={xParam}
      yParam={yParam}
      timepoint={timepoint}
      currentTimepoint={currentTimepoint}
      setTimepoint={setTimepoint}
      allTimepoints={Object.keys(groupedTimepoints)}
      subsetParam={subtypeParam}
      setSelected={setSelected}
      setHighlighted={setHighlighted}
    />
  );
};

const UMAP = ({
  chartDim,
  chartName,
  data,
  genes,
  highlighted,
  xParam,
  yParam,
  subsetParam,
  setHighlighted,
  setSelected,
  timepoint,
  currentTimepoint,
  allTimepoints,
  setTimepoint,
}) => {
  var mousePos = { x: 0, y: 0 };
  var lassPath = "";
  var polyList = [];

  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - LEGEND_WIDTH - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;

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

  const subsetGroups = _.groupBy(data, subsetParam);
  const subsetValues = Object.keys(subsetGroups).sort();

  const subsetColors = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(
      COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
    );

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      saveContext(context);
      drawPoints(
        context,
        data.filter((row) => row[timepoint] === currentTimepoint),
        xScale,
        yScale,
        xParam,
        yParam,
        subsetParam,
        highlighted,
        subsetColors
      );

      drawSubsetLabels(
        context,
        subsetGroups,
        xScale,
        yScale,
        xParam,
        yParam,
        highlighted
      );
      drawUMAPAxis(context, chartHeight, xParam, yParam);
      appendEventListenersToCanvas(context, data, setSelected);
    },
    canvasWidth,
    canvasHeight,
    [highlighted, data]
  );
  //taken from
  //https://gist.github.com/maxogden/574870
  function isPointInPoly(poly, pt) {
    for (var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
      ((poly[i][1] <= pt.y && pt.y < poly[j][1]) ||
        (poly[j][1] <= pt.y && pt.y < poly[i][1])) &&
        pt.x <
          ((poly[j][0] - poly[i][0]) * (pt.y - poly[i][1])) /
            (poly[j][1] - poly[i][1]) +
            poly[i][0] &&
        (c = !c);
    return c;
  }
  function setMousePosition(e, boundingRect) {
    mousePos.x = e.clientX - boundingRect.left;
    mousePos.y = e.clientY - boundingRect.top;
  }
  function setD3MousePosition(event, boundingRect) {
    mousePos.x = event.pageX - boundingRect.left;
    mousePos.y = event.pageY - boundingRect.top;
  }
  function drawLasso(context, boundingRect, x, y) {
    // mouse left button must be pressed
    if (d3.event.buttons !== 1) return;

    context.lineCap = "round";
    context.strokeStyle = "#c0392b";
    context.lineWidth = 3;

    context.moveTo(mousePos.x, mousePos.y);
    polyList = [...polyList, [mousePos.x, mousePos.y]];
    lassPath +=
      lassPath === ""
        ? "M " + mousePos.x + " " + mousePos.y + " "
        : "L " + mousePos.x + " " + mousePos.y + " ";

    //setMousePosition(e, boundingRect);
    setD3MousePosition(d3.event, boundingRect);
    context.lineTo(mousePos.x, mousePos.y);
    lassPath += "L " + mousePos.x + " " + mousePos.y + " ";
    polyList = [...polyList, [mousePos.x, mousePos.y]];

    context.stroke();
  }

  const appendEventListenersToCanvas = (context, data, setSelected) => {
    const scatterSelection = d3.select("#umapSelection");

    d3.select("#umapCanvas")
      .on("mousemove", function mousemove(e) {
        drawLasso(
          context,
          this.getBoundingClientRect(),
          d3.event.pageX,
          d3.event.pageY
        );
      })
      .on("mousedown", function mousedown() {
        context.restore();
        context.beginPath();
        polyList = [];

        setD3MousePosition(d3.event, this.getBoundingClientRect());
        //setMousePosition(e, this.getBoundingClientRect());
      })
      .on("mouseenter", function mouseenter(e) {
        setD3MousePosition(d3.event, this.getBoundingClientRect());
        //setMousePosition(e, this.getBoundingClientRect());
      })
      .on("mouseup", function mouseup(e) {
        var selectedNodes = [];

        data.map((point) => {
          const cords = {
            x: xScale(parseFloat(point[xParam])),
            y: yScale(parseFloat(point[yParam])),
          };
          if (isPointInPoly(polyList, cords)) {
            selectedNodes = [...selectedNodes, point];
          }
        });
        const lassoIDs = selectedNodes.reduce((final, curr) => {
          final[curr["cell_id"]] = true;
          return final;
        }, {});

        context.restore();

        if (selectedNodes.length > 0) {
          setSelected({ selected: lassoIDs });
        }

        context.save();

        context.strokeWidth = 1;
        context.globalAlpha = 0.2;
        context.strokeStyle = "purple";
        context.fillStyle = "black";

        context.fill();
        var lassoSvgPath = new Path2D(lassPath);
        lassoSvgPath.closePath();
        context.fill(lassoSvgPath, "evenodd");
        context.fill();
        context.closePath();
        context.stroke();
        lassPath = "";
      });
  };

  const svgRef = useD3((svg) => {
    drawLegend(
      svg,
      subsetValues,
      subsetColors,
      canvasHeight,
      highlighted,
      setHighlighted,
      subsetParam
    );
  }, []);

  const classes = useStyles();
  return (
    <Grid
      container
      direction="row"
      justify="flex-start"
      alignItems="stretch"
      style={{ marginTop: 15 }}
    >
      <Grid direction="column" justify="flex-start" alignItems="stretch">
        <Card
          className={classes.root}
          style={{
            width: chartDim["width"],
            marginBottom: "20px",
          }}
        >
          <CardActions disableSpacing>
            <Grid
              container
              direction="row"
              justify="space-between"
              alignItems="flex-start"
            >
              <ButtonGroup
                size="small"
                aria-label="small outlined button group"
                className={classes.buttonGroup}
              >
                {allTimepoints.map((timepoint, index) => (
                  <Button
                    onClick={() => {
                      setTimepoint(timepoint);
                    }}
                    style={{
                      opacity: currentTimepoint === timepoint ? 1 : 0.5,
                      color: "black",
                    }}
                  >
                    {timepoint}
                  </Button>
                ))}
              </ButtonGroup>
              <Typography variant="h6">{chartName}</Typography>
              <Button
                className={classes.button}
                variant="outlined"
                size="small"
                disabled={highlighted === null}
                onClick={() => {
                  setSelected(null);
                }}
              >
                Clear Selection
              </Button>
            </Grid>
          </CardActions>
        </Card>
        <Paper
          style={{
            width: chartDim["width"],
            margin: "2px 10px 10px 10px",
            padding: "10px 0px 10px 15px",
          }}
        >
          <Grid
            container
            direction="column"
            justify="flex-start"
            alignItems="stretch"
          >
            <Grid container direction="row" style={{ padding: 0 }}>
              <Grid item>
                <canvas ref={canvasRef} id="umapCanvas" />
              </Grid>
              <Grid item>
                <svg
                  ref={svgRef}
                  id="umapSelection"
                  style={{ yOverflow: "scroll" }}
                />
              </Grid>
            </Grid>
          </Grid>
        </Paper>
      </Grid>
      {highlighted && (
        <Summary
          data={Object.keys(genes).reduce((final, key) => {
            if (highlighted[key]) {
              final[key] = genes[key];
            }
            return final;
          }, [])}
          totalCellCount={data.length}
          colors={subsetColors}
          metadata={data.filter((cell) => highlighted[cell["cell_id"]])}
          chartDim={{
            width: chartDim["width"] / 2,
          }}
        />
      )}
    </Grid>
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
  colorScale,
  highlightedIds
) => {
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;
  data.forEach((point) => {
    if (highlighted) {
      context.fillStyle = highlighted[point["cell_id"]]
        ? colorScale(point[subsetParam])
        : NULL_POINT_COLOR;
    } else {
      context.fillStyle = isHighlighted(point[subsetParam], highlighted)
        ? colorScale(point[subsetParam])
        : NULL_POINT_COLOR;
    }

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
  if (highlighted || highlightedIds) {
    data
      .filter((row) => highlighted && highlighted[row["cell_id"]])
      .forEach((point) => {
        context.fillStyle = colorScale(point[subsetParam]);
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
  }
};

const drawUMAPAxis = (context, chartHeight, xParam, yParam) => {
  context.beginPath();

  const START_X = AXIS_SPACE / 2;
  const START_Y = chartHeight + AXIS_SPACE / 2;

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

const drawSubsetLabels = (
  context,
  subsetGroups,
  xScale,
  yScale,
  xParam,
  yParam,
  highlighted
) => {
  const subsetValues = Object.keys(subsetGroups);

  subsetValues.forEach((subset) => {
    const subsetData = subsetGroups[subset];

    const filteredData = filterOutliers(subsetData, xParam, yParam);
    const { xMin, xMax, yMin, yMax } = getBoxBounds(
      filteredData,
      xParam,
      yParam
    );

    const [x1, x2, y1, y2] = [
      xScale(xMin),
      xScale(xMax),
      yScale(yMin),
      yScale(yMax),
    ];
    const width = Math.abs(x2 - x1);
    const height = Math.abs(y2 - y1);

    context.font = "500 12px Helvetica";
    const textWidth = context.measureText(subset).width;
    context.globalAlpha = isHighlighted(subset, highlighted) ? 0.8 : 0.2;
    context.fillStyle = "white";
    context.fillRect(
      x1 + (width - textWidth) / 2 - 1,
      y1 - height / 2 - 7,
      textWidth + 2,
      14
    );

    context.globalAlpha = isHighlighted(subset, highlighted) ? 1 : 0.2;
    context.textAlign = "center";
    context.textBaseline = "middle";
    context.fillStyle = "black";
    context.fillText(subset, x1 + width / 2, y1 - height / 2);
  });
};

const filterOutliers = (data, xParam, yParam) => {
  const xValues = data
    .map((datum) => parseFloat(datum[xParam]))
    .sort((a, b) => a - b);
  const yValues = data
    .map((datum) => parseFloat(datum[yParam]))
    .sort((a, b) => a - b);

  const [xMin, xMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(xValues, range)
  );
  const [yMin, yMax] = PERCENTILE_RANGE.map((range) =>
    quantileSorted(yValues, range)
  );

  const xiqr = Math.abs(xMax - xMin) * 1.5;
  const yiqr = Math.abs(yMax - yMin) * 1.5;

  return data.filter((datum) => {
    const x = datum[xParam];
    const y = datum[yParam];

    return (
      xMin - xiqr <= x &&
      x <= xMax + xiqr &&
      yMin - yiqr <= y &&
      y <= yMax + yiqr
    );
  });
};

const getBoxBounds = (data, xParam, yParam) => {
  const xValues = data.map((datum) => datum[xParam]);
  const yValues = data.map((datum) => datum[yParam]);

  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const yMin = Math.min(...yValues);
  const yMax = Math.max(...yValues);

  return { xMin, xMax, yMin, yMax };
};

const drawLegend = (
  svg,
  subsetValues,
  colors,
  chartHeight,
  highlighted,
  setHighlighted,
  type
) => {
  const mouseEvents = (element) =>
    element.call((element) =>
      element
        .on("mouseenter", function (d) {
          //          d3.event.stopPropagation();
          setHighlighted("mouseenter", d, type);
        })
        .on("mousedown", function (d, i) {
          //      d3.event.stopPropagation();
          setHighlighted("mousedown", d, type);
        })
        .on("mouseout", function (d, i) {
          //  d3.event.stopPropagation();
          setHighlighted("mouseout", d);
        })
    );
  svg.attr("width", LEGEND_WIDTH).attr("height", chartHeight / 3);
  svg.append("text").text(type);
  const subsets = svg.selectAll("g").data(subsetValues);

  subsets
    .enter()
    .append("rect")
    .merge(subsets)
    .attr("width", LEGEND_SQUARE_LENGTH)
    .attr("height", LEGEND_SQUARE_LENGTH)
    .attr("x", 5)
    .attr("y", (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5)
    .attr("fill", (d) => colors(d))
    .call(mouseEvents);

  subsets
    .enter()
    .append("text")
    .merge(subsets)
    .attr("alignment-baseline", "hanging")
    .attr("text-align", "left")
    .attr("font", "Helvetica")
    .attr("font-weight", "500")
    .attr("font-size", "12px")
    .attr("fill", "#000000")
    .attr("x", LEGEND_SQUARE_LENGTH + 10)
    .attr("y", (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5)
    .text((d) => d)
    .call(mouseEvents);
};

const isHighlighted = (datumValue, highlighted) =>
  highlighted === null || datumValue === highlighted;

export default DataWrapper;
