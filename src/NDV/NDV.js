import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";
import Heatmap from "../components/Heatmap/Heatmap";
import ClonotypeExpansion from "./ClonotypeExpansion";
import ProbabilityHistogram from "../components/Bar/ProbabilityHistogram";
import DEGTable from "./DEGTable";

import IconButton from "@material-ui/core/IconButton";
import CloseIcon from "@material-ui/icons/Close";
import { makeStyles } from "@material-ui/core/styles";
import Popper from "@material-ui/core/Popper";
import Typography from "@material-ui/core/Typography";

import Grid from "@material-ui/core/Grid";

import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardHeader from "@material-ui/core/CardHeader";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";

const NDV = ({ data }) => {
  const [selectedSubtype, setSelectedSubtype] = useState(
    initialState["defaultSelectedObject"]
  );
  const [selectedClonotype, setSelectedClonotype] = useState(
    initialState["defaultSelectedObject"]
  );

  const { metadata, probabilities, degs } = data;

  const filteredMetadata = metadata.filter(
    (row) => row[initialState["clonotypeParam"]] !== "None"
  );
  const topTen = Object.entries(
    _.countBy(
      filteredMetadata.map((row) => row[initialState["clonotypeParam"]])
    )
  )
    .sort(([, a], [, b]) => b - a)
    .slice(0, 10);

  const sampleTen = topTen.reduce((final, curr) => {
    final[curr[0]] = curr[1];
    return final;
  }, {});

  const sampleData = metadata.filter(
    (row) =>
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
    "#ed1a1a",
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
      color: colors(clonotype),
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
          clonotypes: clonotypes,
        }}
        reducer={dashboardReducer}
      >
        {(selectedClonotype["selected"] || selectedSubtype["selected"]) && (
          <Popup
            selected={
              selectedClonotype["selected"] || selectedSubtype["selected"]
            }
            setSelected={() => {
              setSelectedClonotype(initialState["defaultSelectedObject"]);
              setSelectedSubtype(initialState["defaultSelectedObject"]);
            }}
            type={selectedClonotype["selected"] ? "Clonotype" : "Subtype"}
          />
        )}
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
            <Layout
              chartName={"UMAP"}
              data={metadata}
              selectedClonotype={selectedClonotype["selected"]}
              hoveredClonotype={selectedClonotype["hover"]}
              setSelectedClonotype={(clonotype) => {
                if (clonotype["selected"]) {
                  setSelectedSubtype(initialState["defaultSelectedObject"]);
                }
                setSelectedClonotype({ ...clonotype });
              }}
            />
            <Layout
              chartName={"SUBTYPEUMAP"}
              data={metadata}
              selectedSubtype={selectedSubtype["selected"]}
              hoveredSubtype={selectedSubtype["hover"]}
              setSelectedSubtype={(subtype) => {
                if (subtype["selected"]) {
                  setSelectedClonotype(initialState["defaultSelectedObject"]);
                }
                setSelectedSubtype({ ...subtype });
              }}
            />
          </Grid>
          <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-start"
          >
            <Heatmap
              data={probabilities}
              chartDim={{
                height: 550,
                width: 750,
              }}
              column={initialState["clonotypeParam"]}
              row={initialState["subtypeParam"]}
              highlightedColumn={
                selectedSubtype["selected"] || selectedSubtype["hover"]
              }
              highlightedRow={
                selectedClonotype["selected"] || selectedClonotype["hover"]
              }
              columnLabels={clonotypeLabels}
              rowTotal={subtypeTotals}
            />
            <DEGTable
              chartName={"TABLE"}
              data={degs}
              selectedSubtype={
                selectedSubtype["selected"]
                  ? selectedSubtype["selected"]
                  : selectedSubtype["hover"]
              }
              chartDim={{
                height: 500,
                width: 750,
              }}
            />
          </Grid>
          <Grid
            item
            container
            direction="row"
            justify="flex-start"
            alignItems="flex-start"
          >
            <ClonotypeExpansion
              chartName={"BARPLOT"}
              data={probabilities}
              dim={{
                height: 475,
                width: 750,
              }}
              selectedSubtype={selectedSubtype}
              selectedClonotype={selectedClonotype}
              setSelectedSubtype={(subtype) => setSelectedSubtype(subtype)}
              setSelectedClonotype={(clonotype) =>
                setSelectedClonotype(clonotype)
              }
            />
            <ProbabilityHistogram
              chartName={"HISTOGRAM"}
              data={probabilities}
              chartDim={{
                chart: {
                  x1: 100,
                  y1: 50,
                  x2: 600,
                  y2: 400,
                },
                height: 475,
                width: 750,
              }}
              binParam={initialState["logXParam"]}
              lineParam={initialState["subtypeParam"]}
              highlightBarParam={initialState["clonotypeParam"]}
              highlightedBar={
                selectedClonotype["hover"] || selectedClonotype["selected"]
              }
              highlightedLine={
                selectedSubtype["hover"] || selectedSubtype["selected"]
              }
            />
          </Grid>
        </Grid>
      </DashboardProvider>
    </MuiThemeProvider>
  );
};
const useStyles = makeStyles({
  root: {
    minWidth: 275,
  },
  poppr: {
    width: 150,
    float: "right",
    right: "100px",
    top: "10px",
    left: "auto",
  },
});
const Popup = ({ selected, setSelected, type }) => {
  const classes = useStyles();
  return (
    <Popper
      open={true}
      placement={"bottom"}
      transition
      className={classes.popper}
    >
      <Card className={classes.root} variant="outlined">
        <CardHeader
          action={
            <IconButton aria-label="settings">
              <CloseIcon onClick={setSelected} />
            </IconButton>
          }
          title={"Selected " + type}
        />
        <CardContent>
          <Typography variant="body">{selected}</Typography>
        </CardContent>
      </Card>
    </Popper>
  );
};
export default NDV;
