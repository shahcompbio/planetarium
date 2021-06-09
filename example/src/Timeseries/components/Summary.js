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

const Summary = ({ data, chartDim }) => {
  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;
  console.log(data);
  const sumData = Object.keys(data).reduce((final, cellId) => {
    const genesPerCell = Object.keys(data[cellId]);
    genesPerCell.map((gene) => {
      final[gene] = final[gene]
        ? final[gene] + data[cellId][gene]
        : data[cellId][gene];
    });
    return final;
  }, {});
  console.log(sumData);
  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;
  const average = 12;
  const classes = useStyles();
  return (
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
          <Card className={classes.root} variant="outlined">
            <CardContent>
              <Typography
                className={classes.title}
                color="textSecondary"
                gutterBottom
              >
                Mean:
              </Typography>
            </CardContent>
          </Card>
          <Card className={classes.root} variant="outlined">
            <CardContent>
              <Typography
                className={classes.title}
                color="textSecondary"
                gutterBottom
              >
                Average:
              </Typography>
            </CardContent>
          </Card>
        </Grid>
      </Grid>
    </Paper>
  );
};

export default Summary;
