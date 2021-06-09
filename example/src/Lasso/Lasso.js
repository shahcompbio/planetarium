import React, { useState } from "react";
import _ from "lodash";
import "./index.css";

import SubtypeUMAP from "./components/subTypeUmap";
import Summary from "./components/Summary";

import { makeStyles } from "@material-ui/core/styles";
import Popper from "@material-ui/core/Popper";
import Typography from "@material-ui/core/Typography";

import Grid from "@material-ui/core/Grid";

import Button from "@material-ui/core/Button";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardHeader from "@material-ui/core/CardHeader";

import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import CssBaseline from "@material-ui/core/CssBaseline";

import { CONSTANTS, CLONOTYPE_COLORS } from "./config";

const NULL_SELECTED = {
  hover: null,
  selected: null,
};

const DataWrapper = ({ data }) => (
  <Lasso
    metadata={data["metadata"]}
    probabilities={data["probabilities"]}
    degs={data["degs"]}
    genes={data["genes"]}
  />
);

export const Lasso = ({ metadata, probabilities, degs, genes }) => {
  const [isLassoSelected, setIsLassoSelected] = useState(false);
  const [selectedSubtype, setSelectedSubtype] = useState(NULL_SELECTED);
  const [selectedClonotype, setSelectedClonotype] = useState(NULL_SELECTED);

  const { clonotypeParam, subtypeParam, logProbParam } = CONSTANTS;

  // Remove none
  const clonotypeCounts = _.countBy(
    metadata.filter((datum) => datum[clonotypeParam] !== "None"),
    clonotypeParam
  );

  const clonotypeLabels = Object.keys(clonotypeCounts)
    .sort((a, b) => clonotypeCounts[b] - clonotypeCounts[a])
    .slice(0, 10)
    .map((value, index) => ({
      value,
      label: `SEQ${index + 1} - ${value}`,
      color: CLONOTYPE_COLORS[index],
    }));

  const subtypeTotals = _.countBy(metadata, subtypeParam);

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      {(selectedClonotype["selected"] || selectedSubtype["selected"]) && (
        <Popup
          selected={
            selectedClonotype["selected"] || selectedSubtype["selected"]
          }
          setSelected={() => {
            setSelectedClonotype(NULL_SELECTED);
            setSelectedSubtype(NULL_SELECTED);
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
          <Grid item>
            <SubtypeUMAP
              chartName={"SUBTYPEUMAP"}
              data={metadata}
              selectedSubtype={selectedSubtype["selected"]}
              hoveredSubtype={selectedSubtype["hover"]}
              setSelectedSubtype={(subtype) => {
                if (subtype["selected"]) {
                  setSelectedClonotype(NULL_SELECTED);
                }
                setSelectedSubtype((prevState) => ({
                  ...prevState,
                  ...subtype,
                }));
              }}
              chartDim={{
                width: 750,
                height: 600,
              }}
            />
          </Grid>
          <Grid item>
            <Summary
              data={genes}
              chartDim={{
                width: 450,
                height: 600,
              }}
            />
          </Grid>
        </Grid>
      </Grid>
    </MuiThemeProvider>
  );
};
const useStyles = makeStyles({
  root: {
    minWidth: 275,
  },
  header: { padding: 10, paddingBottom: 0 },
  body: {
    padding: 5,
    paddingLeft: 20,
  },
  button: { margin: 5, float: "right" },
  poppr: {
    width: 150,
    float: "right",
    right: "100px",
    top: "10px",
    left: "auto",
    margin: 10,
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
        <CardHeader className={classes.header} title={"Selected " + type} />
        <CardContent className={classes.body}>
          <Typography variant="body">{selected}</Typography>
        </CardContent>
        <Button
          color="primary"
          size="small"
          variant="outlined"
          className={classes.button}
          onClick={setSelected}
        >
          Clear
        </Button>
      </Card>
    </Popper>
  );
};
export default DataWrapper;
