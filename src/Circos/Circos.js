import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";

import Network from "./Chord.js";

import makeStyles from '@mui/styles/makeStyles';

import Typography from "@mui/material/Typography";
import Grid from "@mui/material/Grid";

import { theme } from "../theme/theme.js";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

const Circos = ({ data }) => {
  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
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
            justifyContent="flex-start"
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
      </ThemeProvider>
    </StyledEngineProvider>
  );
};

export default Circos;
