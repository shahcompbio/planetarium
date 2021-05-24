import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";

import Network from "./Chord.js";

import { makeStyles } from "@material-ui/core/styles";

import Typography from "@material-ui/core/Typography";
import Grid from "@material-ui/core/Grid";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";

const Circos = ({ data }) => {
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <DashboardProvider
        initialState={{
          ...initialState,
        }}
        reducer={dashboardReducer}
      >
        <Grid
          container
          direction="column"
          justify="flex-start"
          alignItems="flex-start"
        >
          <Grid style={{ margin: 15 }}>
            <Network
              chartName={"CIRCOS"}
              data={data}
              chartDim={{
                height: 800,
                width: 950,
              }}
            />
          </Grid>
        </Grid>
      </DashboardProvider>
    </MuiThemeProvider>
  );
};

export default Circos;
