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
import Doughnut from "./Doughnut.js";

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
    textAlign: "center",
  },
  cellCountContent: {
    margin: "0px",
    padding: "0px 0px",
    paddingBottom: "0px !important",
    textAlign: "center",
  },
  title: { marginBottom: 0 },
  button: { float: "right" },
});

const Summary = ({ metadata, data, chartDim, colors }) => {
  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;
  const [topNum, setTopNum] = useState(10);

  console.log(metadata);
  const sumData = Object.keys(data).reduce((final, cellId) => {
    const genesPerCell = Object.keys(data[cellId]);
    genesPerCell.map((gene) => {
      if (!final.hasOwnProperty(gene)) {
        final[gene] = { value: 0, count: 0 };
      }
      final[gene]["value"] = final[gene]["value"] + data[cellId][gene];
      final[gene]["count"] = final[gene]["count"] + 1;
    });
    return final;
  }, {});
  const average = Object.keys(sumData).reduce((final, gene) => {
    const avg = sumData[gene]["value"] / sumData[gene]["count"];
    final[gene] = avg;
    return final;
  }, {});

  const topTenAverge = Object.keys(average)
    .sort((a, b) => average[b] - average[a])
    .slice(0, topNum);

  const cloneBreakdown = _.groupBy(metadata, (d) => d["clone"]);
  const cellCount = metadata.length;

  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const classes = useStyles();
  return (
    <Paper
      style={{
        margin: "5px 10px 10px 10px",
        padding: "10px 10px",
        //width: chartDim["width"],
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="stretch"
      >
        <Card className={classes.root} variant="outlined">
          <CardContent className={classes.cellCountContent}>
            <Typography
              display="inline"
              className={classes.title}
              variant="h6"
              gutterBottom
            >
              Cell Count: {cellCount}
            </Typography>{" "}
          </CardContent>
        </Card>
        <Card className={classes.root} variant="outlined">
          <CardContent>
            <Typography
              className={classes.title}
              color="textSecondary"
              gutterBottom
            >
              Top {topNum} expressed genes:{" "}
            </Typography>
            {topTenAverge.map((gene) => (
              <div>
                {gene} - {average[gene]}
              </div>
            ))}
          </CardContent>
        </Card>
        <Card className={classes.root} variant="outlined">
          <CardContent>
            <Typography
              className={classes.title}
              color="textSecondary"
              gutterBottom
            >
              Clone Breakdown:
            </Typography>
            <Doughnut
              data={Object.keys(cloneBreakdown)
                .sort(function (a, b) {
                  return cloneBreakdown[b].length - cloneBreakdown[a].length;
                })
                .reduce((final, clone) => {
                  final = [
                    { key: clone, value: cloneBreakdown[clone] },
                    ...final,
                  ];
                  return final;
                }, [])}
              colors={colors}
              height={300}
              width={chartDim["width"]}
              totalCount={cellCount}
            />
          </CardContent>
        </Card>
      </Grid>
    </Paper>
  );
};

export default Summary;
