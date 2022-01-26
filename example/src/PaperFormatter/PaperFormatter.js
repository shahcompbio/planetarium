import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";

import Formatter from "./components/Formatter";

const DataWrapper = ({}) => {
  return <PaperFormatter x={8} y={11} type={"inches"} />;
};
const PaperFormatter = ({ x, y, type }) => {
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          <Formatter x={x} y={y} type={type} />
        </Block>
      </Grid>
    </MuiThemeProvider>
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

export default DataWrapper;
