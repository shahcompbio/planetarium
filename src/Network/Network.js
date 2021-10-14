import React, { useState, useRef } from "react";
import Layout from "../Layout/Layout";
import _ from "lodash";
import * as d3 from "d3";
import "./App.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";
import TextField from "@mui/material/TextField";
import Autocomplete from '@mui/material/Autocomplete';

import Network from "../components/Network/Network";
import IconButton from "@mui/material/IconButton";
import CloseIcon from "@mui/icons-material/Close";
import makeStyles from '@mui/styles/makeStyles';
import Popper from "@mui/material/Popper";
import Typography from "@mui/material/Typography";

import Grid from "@mui/material/Grid";

import Button from "@mui/material/Button";
import Card from "@mui/material/Card";
import CardContent from "@mui/material/CardContent";
import CardHeader from "@mui/material/CardHeader";
import CardActions from "@mui/material/CardActions";

import { theme } from "../theme/theme.js";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

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
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <CssBaseline />
        <DashboardProvider
          initialState={{
            ...initialState,
          }}
          reducer={dashboardReducer}
        >
          <Grid
            container
            direction="column"
            justifyContent="flex-start"
            alignItems="flex-start"
          >
            <Grid
              item
              container
              direction="row"
              justifyContent="flex-start"
              alignItems="flex-start"
            >
              <Grid style={{ margin: 15 }}>
                <Search
                  data={[...new Set(data.map((row) => row["aliquot"]).flat(1))]}
                  selectOption={(option) => setSelectedAliquot(option)}
                />
                {selectedAliquot && (
                  <Search
                    data={[
                      ...new Set(
                        data
                          .filter((row) => row["aliquot"] === selectedAliquot)
                          .map((row) => row["celltypes"])
                          .flat(1)
                      ),
                    ]}
                    selectOption={(option) => setSelectedCelltype(option)}
                  />
                )}
              </Grid>
              {selectedCelltype && selectedAliquot && (
                <Network
                  chartName={"NETWORK"}
                  data={data
                    .filter((row) => row["aliquot"] === selectedAliquot)
                    .filter((row) => row["celltypes"] === selectedCelltype)}
                  chartDim={{
                    height: 800,
                    width: 950,
                  }}
                />
              )}
            </Grid>
          </Grid>
        </DashboardProvider>
      </ThemeProvider>
    </StyledEngineProvider>
  );
};
const useStyles = makeStyles((theme) => ({
  inputRoot: {
    marginBottom: 15,
    "& .MuiAutocomplete-popupIndicator": { color: "black" },
    color: "black",
    "& .MuiOutlinedInput-notchedOutline": {
      borderColor: "black",
    },
    "&:hover .MuiOutlinedInput-notchedOutline": {
      borderColor: "black",
    },
    "&.Mui-focused .MuiOutlinedInput-notchedOutline": {
      borderColor: "black",
    },
    "& .MuiInputLabel-formControl": {
      color: "black",
    },
  },
}));
const Search = ({ data, selectOption }) => {
  const classes = useStyles();
  return (
    <Autocomplete
      classes={classes}
      options={data}
      getOptionLabel={(option) => option}
      style={{ width: 300 }}
      onChange={(event, option) => {
        selectOption(option);
      }}
      renderInput={(params) => (
        <TextField
          {...params}
          InputLabelProps={{
            style: { color: "#black" },
          }}
          label="Search"
          variant="outlined"
        />
      )}
    />
  );
};
export default Networks;
