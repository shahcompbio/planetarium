import React, { useState, useEffect } from "react";

import Grid from "@material-ui/core/Grid";
import Switch from "@material-ui/core/Switch";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";

import {
  useCanvas,
  useLasso,
  VerticalLegend,
  useGL,
  Layout,
} from "@shahlab/planetarium";
import Reg from "./Reg.js";
import Menu from "./Menu.js";

import { CssBaseline } from "@material-ui/core";

import "./Test.css";
import _ from "lodash";

import * as d3 from "d3";
const PADDING = 10;
const AXIS_SPACE = 20;
const NUM_LEGEND_WIDTH = 70;
const AXIS_LENGTH = 50;
const AXIS_FONT = "Helvetica";
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
];
const DataWrapper = ({ data, filters }) => {
  return <Test data={data} filters={filters} />;
};

export const drawAxis = ({ context, xPos, yPos, xLabel, yLabel }) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.lineWidth = 1;
  context.lineCap = "butt";
  context.moveTo(xPos, yPos);
  context.lineTo(xPos, yPos - AXIS_LENGTH);
  context.stroke();

  context.beginPath();
  context.moveTo(xPos, yPos);
  context.lineTo(xPos + AXIS_LENGTH, yPos);
  context.stroke();

  context.textAlign = "left";
  context.textBaseline = "middle";
  context.fillText(xLabel, xPos + AXIS_LENGTH + 2, yPos);
  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText(yLabel, -(yPos - AXIS_LENGTH - 2), xPos);
  context.restore();
};

const getColorScale = ({ data, subsetParam, isCategorical }) => {
  if (isCategorical) {
    const subsetGroups = _.groupBy(data, subsetParam);
    const subsetValues = Object.keys(subsetGroups).sort();
    return d3
      .scaleOrdinal()
      .domain(subsetValues)
      .range(
        COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
      );
  } else {
    const subsetData = data
      .filter((d) => d.hasOwnProperty(subsetParam))
      .map((d) => parseFloat(d[subsetParam]));

    const subsetMax = Math.max(...subsetData);
    return d3
      .scaleSequential(d3.interpolateViridis)
      .domain([0, subsetMax])
      .nice();
  }
};

const Test = ({
  data,
  filters,
  onLasso = (data) => {},
  disable = false,
  idParam = "cell_id",
  xParam = "UMAP_1",
  yParam = "UMAP_2",
  subsetParam = "cell_type",
  highlightIDs = null,
  onLegendHover = (value) => {},
  onLegendClick = (value) => {},
  labels = (value) => value,
  MoreInfoComponent = () => null,
  layerNames = ["umapCanvas1", "umapCanvas2"],
}) => {
  const [layerFilters, setLayerFilters] = useState({
    0: { filters: null },
    1: { filters: null },
  });

  const [activeLayer, setActiveLayer] = useState(0);
  const [canvas, setCanvas] = useState(null);

  useEffect(() => {
    if (layerFilters[activeLayer]["filters"] !== null) {
      const filters = layerFilters[activeLayer]["filters"];
      console.log(filters);
    }
  }, [layerFilters]);

  const isCategorical =
    typeof data.filter((datum) => datum.hasOwnProperty(subsetParam))[0][
      subsetParam
    ] === "string";

  const legendFilter = isCategorical
    ? (value, datum) => datum[subsetParam] === value
    : (value, datum) =>
        datum.hasOwnProperty(subsetParam) && datum[subsetParam] >= value[0];

  const canvasWidth = 700;
  const canvasHeight = 600;
  const [wrapperRef] = useGL(canvasWidth, canvasHeight, layerNames);

  const legendWidth = NUM_LEGEND_WIDTH;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const yData = data.map((d) => parseFloat(d[yParam]));
  const xData = data.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);
  const subsetColors = getColorScale({ data, subsetParam, isCategorical });

  const xScale = d3
    .scaleLinear()
    .domain([xMin, xMax])
    .range([PADDING, PADDING + chartWidth]);

  const yScale = d3
    .scaleLinear()
    .domain([yMax, yMin])
    .range([PADDING, PADDING + chartHeight]);

  const [hoveredLegend, setHoveredLegend] = useState(null);
  const [clickedLegend, setClickedLegend] = useState(null);

  const legendIDs = hoveredLegend || clickedLegend;

  useEffect(() => {
    const layers = layerNames.map((layer) => {
      return d3.select("#" + layer).node();
    });
    setCanvas(layers);
  }, [wrapperRef]);

  const [lassoData, drawLasso, addLassoHandler, resetLasso] = useLasso(
    data,
    xScale,
    yScale,
    xParam,
    yParam
  );
  const subsettedIDs = [];

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      canvas.id = "lasso";

      drawAxis({
        context,
        xPos: AXIS_SPACE / 2,
        yPos: canvasHeight - AXIS_SPACE / 2,
        xLabel: xParam,
        yLabel: yParam,
      });

      drawLasso(context);
      const disableLasso = disable || legendIDs !== null;
      addLassoHandler(canvas, disableLasso, onLasso);
    },
    canvasWidth,
    canvasHeight,
    [data, disable, subsetParam, subsettedIDs]
  );

  const getLegendData = (value) => {
    return value === null
      ? value
      : data
          .filter((datum) => legendFilter(value, datum))
          .map((datum) => datum[idParam]);
  };

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid style={{ height: canvasHeight, width: canvasWidth }}>
        <Switch
          checked={Boolean(activeLayer)}
          onChange={(event) => {
            setActiveLayer(!activeLayer);
          }}
        />

        <Layout title={""} infoText={""}>
          <Grid
            container
            direction="row"
            style={{ padding: 0, position: "relative", height: canvasHeight }}
          >
            <Grid item style={{ position: "relative" }} ref={wrapperRef}>
              {canvas && (
                <Reg
                  canvasRef={canvas[activeLayer]}
                  pointSize={activeLayer + 2}
                  data={data}
                  width={700}
                  height={600}
                  xParam={xParam}
                  yParam={yParam}
                  xScale={xScale}
                  yScale={yScale}
                  subsetParam={subsetParam}
                />
              )}
              <canvas ref={canvasRef} style={{ zIndex: 100 }} />
            </Grid>
            <Grid
              item
              style={{
                paddingLeft: "40px",
                float: "right",
                position: "absolute",
                right: 0,
              }}
            >
              <VerticalLegend
                width={legendWidth}
                height={canvasHeight / 2}
                colorScale={subsetColors}
                ticks={
                  isCategorical
                    ? subsetColors
                        .domain()
                        .sort()
                        .map((value) => ({ value, label: labels(value) }))
                    : 10
                }
                onHover={(value) => {
                  const legendData = getLegendData(value);
                  setHoveredLegend(legendData);
                  onLegendHover(value);
                }}
                onClick={(value) => {
                  const legendData = getLegendData(value);
                  setClickedLegend(legendData);
                  onLegendClick(value);
                }}
                disable={disable || lassoData !== null}
                reset={highlightIDs !== null}
              />
              <MoreInfoComponent />
            </Grid>
          </Grid>
        </Layout>
        <Menu
          filters={filters}
          activeLayer={activeLayer}
          layerFilters={layerFilters}
          setLayerFilters={(filters) => setLayerFilters(filters)}
        />
      </Grid>
    </MuiThemeProvider>
  );
};

export default DataWrapper;
