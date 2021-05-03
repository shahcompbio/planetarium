import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";
import TextField from "@material-ui/core/TextField";
import Autocomplete from "@material-ui/lab/Autocomplete";

import Network from "../components/Network/Network";
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

const Networks = ({ data }) => {
  /*  const [selectedAliquot, setSelectedAliquot] = useState(
    "SPECTRUM-OV-007_PELVIC_PERITONEUM"
  );
  const [selectedCelltype, setSelectedCelltype] = useState(
    "Fibroblast-Endothelial.cell"
  );*/
  const [selectedAliquot, setSelectedAliquot] = useState(null);
  const [selectedCelltype, setSelectedCelltype] = useState(null);
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <DashboardProvider
        initialState={{
          ...initialState
        }}
        reducer={dashboardReducer}
      >
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
            <Grid style={{ margin: 15 }}>
              <Search
                data={[...new Set(data.map(row => row["aliquot"]).flat(1))]}
                selectOption={option => setSelectedAliquot(option)}
              />
              {selectedAliquot && (
                <Search
                  data={[
                    ...new Set(
                      data
                        .filter(row => row["aliquot"] === selectedAliquot)
                        .map(row => row["celltypes"])
                        .flat(1)
                    )
                  ]}
                  selectOption={option => setSelectedCelltype(option)}
                />
              )}
            </Grid>
            {selectedCelltype && selectedAliquot && (
              <Network
                chartName={"NETWORK"}
                data={data
                  .filter(row => row["aliquot"] === selectedAliquot)
                  .filter(row => row["celltypes"] === selectedCelltype)}
                chartDim={{
                  height: 800,
                  width: 950
                }}
              />
            )}
          </Grid>
        </Grid>
      </DashboardProvider>
    </MuiThemeProvider>
  );
};
const useStyles = makeStyles(theme => ({
  inputRoot: {
    marginBottom: 15,
    "& .MuiAutocomplete-popupIndicator": { color: "black" },
    color: "black",
    "& .MuiOutlinedInput-notchedOutline": {
      borderColor: "black"
    },
    "&:hover .MuiOutlinedInput-notchedOutline": {
      borderColor: "black"
    },
    "&.Mui-focused .MuiOutlinedInput-notchedOutline": {
      borderColor: "black"
    },
    "& .MuiInputLabel-formControl": {
      color: "black"
    }
  }
}));
const Search = ({ data, selectOption }) => {
  const classes = useStyles();
  return (
    <Autocomplete
      classes={classes}
      options={data}
      getOptionLabel={option => option}
      style={{ width: 300 }}
      onChange={(event, option) => {
        selectOption(option);
      }}
      renderInput={params => (
        <TextField
          {...params}
          InputLabelProps={{
            style: { color: "#black" }
          }}
          label="Search"
          variant="outlined"
        />
      )}
    />
  );
};
export default Networks;
