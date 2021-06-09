import React, { useState } from "react";
import _ from "lodash";
import "./index.css";
import * as d3Array from "d3-array";

import TimeSeriesUMAP from "./components/TimeSeriesUmap";
import Summary from "./components/Summary";

import { makeStyles } from "@material-ui/core/styles";
import Popper from "@material-ui/core/Popper";
import Typography from "@material-ui/core/Typography";

import Grid from "@material-ui/core/Grid";

import Button from "@material-ui/core/Button";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardHeader from "@material-ui/core/CardHeader";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";

import { CONSTANTS, CLONOTYPE_COLORS } from "./config";

const NULL_SELECTED = {
  hover: null,
  selected: null,
};

const DataWrapper = ({ data }) => (
  <TimeSeries metadata={data["metadata"]} genes={data["genes"]} />
);

export const TimeSeries = ({ metadata, genes }) => {
  const [selected, setSelected] = useState(NULL_SELECTED);
  const { clonotypeParam, subtypeParam, logProbParam, timepoint } = CONSTANTS;

  // Remove none
  const clonotypeCounts = _.countBy(
    metadata.filter((datum) => datum[clonotypeParam] !== "None"),
    clonotypeParam
  );

  const clonotypeLabels = Object.keys(clonotypeCounts)
    .sort((a, b) => clonotypeCounts[b] - clonotypeCounts[a])
    .slice(0, 10)
    .map((value, index) => ({
      value,
      label: `SEQ${index + 1} - ${value}`,
      color: CLONOTYPE_COLORS[index],
    }));

  const subtypeTotals = _.countBy(metadata, subtypeParam);

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="flex-start"
      >
        <Grid
          item
          container
          direction="row"
          justify="flex-start"
          alignItems="flex-start"
        >
          <TimeSeriesUMAP
            chartName={"SUBTYPEUMAP"}
            data={metadata}
            genes={genes}
            selected={selected["selected"]}
            hovered={selected["hover"]}
            setSelected={(selection) => {
              if (selection === null) {
                setSelected(NULL_SELECTED);
              } else {
                setSelected({ ...selection });
              }
            }}
            chartDim={{
              width: 750,
              height: 600,
            }}
          />
        </Grid>
      </Grid>
    </MuiThemeProvider>
  );
};

export default DataWrapper;
