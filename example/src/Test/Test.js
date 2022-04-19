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
import { matchSorter } from "match-sorter";

import { CssBaseline } from "@material-ui/core";

import "./Test.css";
import _ from "lodash";

import * as d3 from "d3";
import { drawAxis } from "./util.js";

const PADDING = 10;
const AXIS_SPACE = 20;
const NUM_LEGEND_WIDTH = 70;
const AXIS_LENGTH = 50;
const AXIS_FONT = "Helvetica";
const AXIS_COLOR = "#000000";
const canvasWidth = 700;
const canvasHeight = 600;
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
const testFilters = {
  patient: ["K", "Q", "O", "P"],
  response: ["PR/SD"],
  cell_type: ["CD4", "B Cell", "CD8"],
};
const DataWrapper = ({ data, filters, idParam = "cell_id" }) => {
  const [filteredData, setFilteredData] = useState([...data]);
  const [layerFilters, setLayerFilter] = useState({
    0: { filters: { ...testFilters }, lasso: null },
    1: { filters: null, lasso: null },
  });
  const [activeLayer, setActiveLayer] = useState(0);

  const changeLayer = (event) => {
    const newActive = activeLayer ? 0 : 1;
    setActiveLayer(newActive);
    const filter = layerFilters[newActive]["filters"];
    const newFilteredData =
      filter === null
        ? data
        : getDataFromFilters(layerFilters[newActive]["filters"], data, idParam);
    setFilteredData([...newFilteredData]);
  };

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid style={{ height: canvasHeight, width: canvasWidth }}>
        <Switch
          checked={Boolean(activeLayer)}
          onChange={(event) => changeLayer(event)}
        />
        {filteredData.length > 0 ? (
          <Test
            data={data}
            filteredData={filteredData}
            setFilteredData={setFilteredData}
            activeLayer={activeLayer}
            layerFilters={layerFilters}
            filters={filters}
            onLasso={(lassoData) => {
              //console.log(lassoData);
              //console.log(layerFilters);

              var newLayerFilter = layerFilters;
              var newFilteredData = lassoData;
              if (lassoData !== null) {
                newLayerFilter[activeLayer]["lasso"] = {
                  selectionLength: lassoData.length,
                  selectedFrom: filteredData.length,
                };
              } else {
                //removed lasso from plot
                newLayerFilter[activeLayer]["lasso"] = null;
                newFilteredData = getDataFromFilters(
                  layerFilters[activeLayer]["filters"],
                  data,
                  idParam
                );
              }

              setLayerFilter({ ...newLayerFilter });
              setFilteredData([...newFilteredData]);
            }}
          />
        ) : (
          <div>No results</div>
        )}
        <Menu
          filters={filters}
          activeLayer={activeLayer}
          layerFilters={layerFilters}
          setLayerFilters={(filters) => {
            console.log(filters);
            return setLayerFilter({ ...filters });
          }}
        />
      </Grid>
    </MuiThemeProvider>
  );
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
const getDataFromFilters = (filters, data, idParam) => {
  //  console.log(filters);
  return Object.keys(filters).reduce(
    (dataLeft, filter) => {
      const interMatches = filters[filter].reduce((final, curr) => {
        /*or filter*/
        //  console.log(dataLeft);
        //  console.log(curr);
        const matches = matchSorter(dataLeft, curr, {
          maxRanking: matchSorter.rankings.EQUALS,
          keys: [(item) => item[filter]],
        });
        //console.log(matches);
        final = [...final, ...matches];
        return final;
      }, []);
      //  console.log(interMatches);
      /* get all the unique matches by cell id*/
      return [
        ...new Map(interMatches.map((item) => [item[idParam], item])).values(),
      ];
    },
    [...data]
  );
};

const Test = ({
  data,
  filteredData,
  setFilteredData,
  filters,
  activeLayer,
  layerFilters,
  disable = false,
  onLasso = (data) => {},
  idParam = "cell_id",
  xParam = "UMAP_1",
  yParam = "UMAP_2",
  subsetParam = "cell_type",
  highlightIDs = null,
  labels = (value) => value,
  MoreInfoComponent = () => null,
  layerNames = ["umapCanvas1", "umapCanvas2"],
}) => {
  const [canvas, setCanvas] = useState(null);

  useEffect(() => {
    console.log("edddditf");
    console.log(layerFilters[activeLayer]["filters"]);
    if (layerFilters[activeLayer]["filters"] !== null) {
      const newFilters = layerFilters[activeLayer]["filters"];
      console.log(newFilters);
      console.log(data);
      const newFilteredData = getDataFromFilters(newFilters, data, idParam);
      console.log(newFilteredData);
      setFilteredData([...newFilteredData]);
    }
  }, [layerFilters]);

  const isCategorical =
    typeof filteredData.filter((datum) => datum.hasOwnProperty(subsetParam))[0][
      subsetParam
    ] === "string";

  const legendFilter = isCategorical
    ? (value, datum) => datum[subsetParam] === value
    : (value, datum) =>
        datum.hasOwnProperty(subsetParam) && datum[subsetParam] >= value[0];

  const [wrapperRef] = useGL(canvasWidth, canvasHeight, layerNames);

  const legendWidth = NUM_LEGEND_WIDTH;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const yData = filteredData.map((d) => parseFloat(d[yParam]));
  const xData = filteredData.map((d) => parseFloat(d[xParam]));

  const yMin = Math.min(...yData);
  const yMax = Math.max(...yData);
  const xMin = Math.min(...xData);
  const xMax = Math.max(...xData);

  const subsetColors = getColorScale({
    data: filteredData,
    subsetParam,
    isCategorical,
  });

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
    filteredData,
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
    [filteredData, disable, subsetParam, subsettedIDs]
  );

  const getLegendData = (value) => {
    return value === null
      ? value
      : filteredData
          .filter((datum) => legendFilter(value, datum))
          .map((datum) => datum[idParam]);
  };

  return (
    <Layout title={""} infoText={""}>
      <Grid
        container
        direction="row"
        style={{ padding: 0, position: "relative", height: canvasHeight }}
      >
        <Grid
          item
          style={{ position: "relative", padding: 0 }}
          ref={wrapperRef}
        >
          {canvas && (
            <Reg
              canvasRef={canvas[activeLayer]}
              pointSize={activeLayer + 4}
              data={filteredData}
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
              //  onLegendHover(value);
            }}
            onClick={(value) => {
              const legendData = getLegendData(value);
              setClickedLegend(legendData);
              //  onLegendClick(value);
            }}
            disable={disable || lassoData !== null}
            reset={highlightIDs !== null}
          />
          <MoreInfoComponent />
        </Grid>
      </Grid>
    </Layout>
  );
};

export default DataWrapper;
