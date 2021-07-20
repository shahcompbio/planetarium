import React, { useState, useEffect } from "react";
import _ from "lodash";
import "./index.css";

import * as d3 from "d3";

import { UMAP, Fishtail, Layout, sortAlphanumeric } from "@shahlab/planetarium";
import TopGenesPanel from "./components/TopGenesPanel";
import Doughnut from "./components/Doughnut";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import Typography from "@material-ui/core/Typography";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";

const SUBSET_PARAM = "clone";
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

const DataWrapper = ({ dashboardID, api, data }) => {
  return (
    <TimeSeries
      metadata={data["metadata"]}
      genes={data["genes"]}
      dashboardID={dashboardID}
      api={api}
    />
  );
};

const Block = ({ children }) => (
  <Grid item>
    <Paper
      style={{
        margin: 10,
        padding: 10,
      }}
    >
      {children}
    </Paper>
  </Grid>
);

export const TimeSeries = ({ metadata, dashboardID, api }) => {
  const [hoveredTimepoint, setHoveredTimepoint] = useState(null);
  const [timepoint, setTimepoint] = useState(metadata[0]["timepoint"]);
  const [hoveredClone, setHoveredClone] = useState(null);
  const [clone, setClone] = useState(null);
  const [activeGraph, setActiveGraph] = useState(null);

  const subsetValues = _.uniq(
    metadata.map((datum) => datum[SUBSET_PARAM])
  ).sort(sortAlphanumeric);

  const colorScale = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(
      COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
    );

  const data = metadata.filter(
    (datum) => datum["timepoint"] === (hoveredTimepoint || timepoint)
  );

  const [highlightedCells, setHighlightedCells] = useState(data);

  useEffect(() => {
    document.title = `Timeseries - ${dashboardID}`;
  });

  useEffect(() => {
    setHighlightedCells(data);
  }, [hoveredTimepoint, timepoint]);

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="column">
        <Grid container direction="row">
          <Grid item container direction="column" xs={2}>
            <Block>
              <Typography color="textSecondary" gutterBottom>
                {dashboardID} ({hoveredTimepoint || timepoint})
              </Typography>
              <div>
                <Typography display="inline" variant="h5" gutterBottom>
                  {highlightedCells.length}
                </Typography>
                <Typography
                  display="inline"
                  variant="h6"
                  color="textSecondary"
                  gutterBottom
                >
                  /{data.length} Cells
                </Typography>
              </div>
            </Block>
            <Block>
              <TopGenesPanel
                ids={highlightedCells.map((datum) => datum["cell_id"])}
                api={api}
                dashboardID={dashboardID}
              />
            </Block>
            <Block>
              <Typography color="textSecondary" gutterBottom>
                Clone Breakdown
              </Typography>
              <Doughnut
                data={highlightedCells}
                width={250}
                height={250}
                subsetParam={SUBSET_PARAM}
                colors={colorScale}
              />
            </Block>
          </Grid>
          <Grid item>
            <Layout title="UMAP" infoText="UMAP">
              <UMAP
                data={data}
                width={750}
                height={650}
                xParam="umap_1"
                yParam="umap_2"
                subsetParam={SUBSET_PARAM}
                subset={hoveredClone || clone}
                idParam="cell_id"
                colorScale={colorScale}
                onLegendHover={setHoveredClone}
                onLegendClick={(value) => {
                  setClone(value);
                  setHighlightedCells(
                    value === null
                      ? data
                      : data.filter((datum) => datum[SUBSET_PARAM] === value)
                  );
                  setActiveGraph(value === null ? null : "umap");
                }}
                onLasso={(value) => {
                  setHighlightedCells(value);
                  setActiveGraph(value.length === data.length ? null : "umap");
                }}
                disable={activeGraph !== null && activeGraph !== "umap"}
              />
            </Layout>
          </Grid>
        </Grid>
        <Grid container direction="row">
          <Layout
            title="Timeseries expansion"
            infoText="Clonal expansion over time"
          >
            <Fishtail
              data={metadata}
              width={1100}
              height={400}
              subsetParam={SUBSET_PARAM}
              timepoint={timepoint}
              subset={hoveredClone || clone}
              colorScale={colorScale}
              onTimepointClick={setTimepoint}
              onTimepointHover={setHoveredTimepoint}
              onLegendHover={setHoveredClone}
              onLegendClick={(value) => {
                setClone(value);
                setActiveGraph(value === null ? null : "fishtail");
              }}
              disable={activeGraph !== null && activeGraph !== "fishtail"}
            />
          </Layout>
        </Grid>
      </Grid>
    </MuiThemeProvider>
  );
};

export default DataWrapper;
