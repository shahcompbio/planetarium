import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import * as d3Slider from "d3-simple-slider";
import { quantileSorted } from "d3-array";
import _ from "lodash";
import { InfoBar, useCanvas, useD3 } from "@shahlab/planetarium";
import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import Card from "@material-ui/core/Card";
import CardActions from "@material-ui/core/CardActions";
import CardContent from "@material-ui/core/CardContent";
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";

import { makeStyles } from "@material-ui/core/styles";
import { CONSTANTS } from "./config";

const PADDING = 10;
const TITLE_HEIGHT = 30;

const LEGEND_WIDTH = 180;
const AXIS_SPACE = 20;

const AXIS_COLOR = "#000000";

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
});

const Slider = ({
  chartDim,
  chartName,
  data,
  currentTimepoint,
  setTimepoint,
}) => {
  const [context, saveContext] = useState();
  const canvasWidth = chartDim["width"] - LEGEND_WIDTH - PADDING - PADDING;
  const canvasHeight = chartDim["height"] - TITLE_HEIGHT;
  const { timepoint } = CONSTANTS;
  const chartWidth = canvasWidth - AXIS_SPACE;
  const chartHeight = canvasHeight - AXIS_SPACE - PADDING - PADDING;

  const svgRef = useD3((svg) => {
    drawSlider(svg, canvasHeight, timepoint, setTimepoint);
  }, []);
  const drawSlider = () => {
    // Step
    const scale = d3.scaleLinear()
    .domain([0,1])
    .range(Object.keys(data));

    var sliderStep = d3Slider
      .sliderBottom()
      .min(obdata, (d) => d[timepoint]))
      .max(d3.max(data, (d) => d[timepoint]))
      .width(300)
      .tickFormat(d3.format(".2%"))
      .ticks(5)
      .step(0.005)
      .default(0.015)
      .on("onchange", (val) => {
        d3.select("p#value-step").text(d3.format(".2%")(val));
      });

    var gStep = d3
      .select("#timeSeriesSlider")
      .append("g")
      .attr("transform", "translate(30,30)");

    gStep.call(sliderStep);

    d3.select("p#value-step").text(d3.format(".2%")(sliderStep.value()));
  };
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
            <Grid item>
              <svg
                ref={svgRef}
                id="timeSeriesSlider"
                style={{ yOverflow: "scroll" }}
              />
            </Grid>
          </Grid>
        </Grid>
      </Paper>
    </Grid>
  );
};
export default Slider;
