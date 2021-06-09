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
import Typography from "@material-ui/core/Typography";

import { makeStyles } from "@material-ui/core/styles";
import { CONSTANTS } from "../config";

const PADDING = 10;
const TITLE_HEIGHT = 30;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const AXIS_COLOR = "#000000";

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
});

const DataWrapper = ({
  chartName,
  chartDim,
  data,
  isLassoSelected,
  setIsLassoSelected,
}) => {
  return (
    <Summary
      data={data}
      chartDim={chartDim}
      isLassoSelected={isLassoSelected}
      setIsLassoSelected={setIsLassoSelected}
    />
  );
};

const Summary = ({ data, isLassoSelected, setIsLassoSelected, chartDim }) => {
  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - LEGEND_WIDTH - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;
  console.log(data);
  const classes = useStyles();
  return (
    <Grid
      container
      direction="column"
      justify="flex-start"
      alignItems="stretch"
    >
      <Paper
        style={{
          margin: "2px 10px 10px 10px",
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
          <Grid container direction="row" style={{ padding: 0 }}>
            <div />
          </Grid>
        </Grid>
      </Paper>
    </Grid>
  );
};

export default DataWrapper;
