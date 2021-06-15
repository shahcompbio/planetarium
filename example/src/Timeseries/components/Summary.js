import React, { useState } from "react";
import * as d3 from "d3";
import _ from "lodash";
import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import Typography from "@material-ui/core/Typography";
import Doughnut from "./Doughnut.js";

import { makeStyles } from "@material-ui/core/styles";

const PADDING = 10;
const TITLE_HEIGHT = 30;

const useStyles = makeStyles({
  root: {
    minWidth: 275,
    margin: "10px 10px 0px 10px",
    padding: "0px 0px 10px 0px",
    float: "right",
    textAlign: "left",
  },
  doughnutCard: {
    margin: "10px 10px 0px 10px",
    float: "right",
    textAlign: "left",
  },
  doughnutContent: {
    paddingBottom: "0 !important",
  },
  geneCard: {
    minWidth: 275,
    margin: "10px 10px 0px 10px",
    padding: "0px 0px 0px 0px",
    float: "right",
    textAlign: "left",
  },
  cellCountTitle: {
    fontWeight: "bold",
    margin: 0,
  },
  cellCountTotalTitle: {
    color: "grey",
  },
  cellCountText: { margin: 0 },
  cellCountContent: {
    margin: "0px",
    padding: "10px 10px",
    paddingBottom: "0px !important",
    textAlign: "left",
  },
  geneTitle: {},
  title: { marginBottom: 0 },
  button: { float: "right" },
});

const Summary = ({ metadata, data, chartDim, colors, totalCellCount }) => {
  const [topNum, setTopNum] = useState(10);
  console.log(data);
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
            <div>
              <Typography
                display="inline"
                className={classes.cellCountTitle}
                variant="h3"
                gutterBottom
              >
                {cellCount}
              </Typography>
              <Typography
                display="inline"
                className={classes.cellCountTotalTitle}
                variant="h4"
                gutterBottom
              >
                /{totalCellCount}
              </Typography>
            </div>
            <div>
              <Typography
                display="inline"
                className={classes.cellCountText}
                variant="h6"
                gutterBottom
              >
                Cell Count
              </Typography>
            </div>
          </CardContent>
        </Card>
        <Card className={classes.geneCard} variant="outlined">
          <CardContent>
            <Typography
              className={classes.geneTitle}
              color="textSecondary"
              gutterBottom
            >
              Top {topNum} expressed genes:{" "}
            </Typography>
            {topTenAverge.map((gene) => (
              <div>
                <span>{d3.format("2s")(average[gene])}</span>
                <span style={{ width: 40, display: "inline-block" }} />
                <span style={{ fontWeight: "bold" }}>{gene}</span>
              </div>
            ))}
          </CardContent>
        </Card>
        <Card className={classes.doughnutCard} variant="outlined">
          <CardContent className={classes.doughnutContent}>
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
              height={180}
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
