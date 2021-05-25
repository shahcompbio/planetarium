import React, { useState, useRef } from "react";

import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";

import IconButton from "@material-ui/core/IconButton";
import CloseIcon from "@material-ui/icons/Close";
import { makeStyles } from "@material-ui/core/styles";
import Popper from "@material-ui/core/Popper";
import Typography from "@material-ui/core/Typography";

import Grid from "@material-ui/core/Grid";

import Button from "@material-ui/core/Button";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardHeader from "@material-ui/core/CardHeader";
import CardActions from "@material-ui/core/CardActions";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";
import SubTypeUmap from "../components/UMAP/subTypeUmap.js";

const App = ({ data }) => {
  const [selectedSubtype, setSelectedSubtype] = useState(
    initialState["defaultSelectedObject"]
  );
  const [selectedClonotype, setSelectedClonotype] = useState(
    initialState["defaultSelectedObject"]
  );

  const { metadata, probabilities, degs } = data;

  const filteredMetadata = metadata.filter(
    row => row[initialState["clonotypeParam"]] !== "None"
  );
  const topTen = Object.entries(
    _.countBy(filteredMetadata.map(row => row[initialState["clonotypeParam"]]))
  )
    .sort(([, a], [, b]) => b - a)
    .slice(0, 10);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = metadata.filter(
    row =>
      sampleTen.hasOwnProperty(row[initialState["clonotypeParam"]]) &&
      row[initialState["clonotypeParam"]] !== "None"
  );
  const topTenNumbering = Object.keys(sampleTen).reduce((final, seq, index) => {
    final[seq] = "SEQ" + (index + 1);
    return final;
  }, {});
  const clonotypes = _.groupBy(sampleData, initialState["clonotypeParam"]);

  const types = Object.keys(clonotypes);

  const colourList = [
    "#674172",
    "#098dde",
    "#fa832f",
    "#0e5702",
    "#c20c1e",
    "#911eb4",
    "#fc97bc",
    "#469990",
    "#b5762a",
    "#5aebed",
    "#8f8f3f",
    "#ed1a1a"
  ];
  var colors = d3
    .scaleOrdinal()
    .domain([...types])
    .range([...colourList]);

  const clonotypeLabels = Object.keys(sampleTen)
    .sort(([, a], [, b]) => b - a)
    .map((clonotype, index) => ({
      value: clonotype,
      label: `SEQ${index + 1} - ${clonotype}`,
      color: colors(clonotype)
    }));

  const subtypeTotals = _.countBy(metadata, initialState["subtypeParam"]);

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <DashboardProvider
        initialState={{
          ...initialState,
          sampleTen: sampleTen,
          sampleData: sampleData,
          topTenNumbering: topTenNumbering,
          topTen: topTen,
          colors: colors,
          clonotypes: clonotypes
        }}
        reducer={dashboardReducer}
      >
        <Grid
          container
          direction="column"
          justify="flex-start"
          alignItems="flex-start"
        >
          <SubTypeUmap
            data={metadata}
            selectedSubtype={selectedSubtype["selected"]}
            hoveredSubtype={selectedSubtype["hover"]}
            setSelectedSubtype={subtype => {
              if (subtype["selected"]) {
                setSelectedClonotype(initialState["defaultSelectedObject"]);
              }
              setSelectedSubtype(prevState => ({
                ...prevState,
                ...subtype
              }));
            }}
            chartDim={{ width: 750, height: 600 }}
          />
        </Grid>
      </DashboardProvider>
    </MuiThemeProvider>
  );
};

export default App;
