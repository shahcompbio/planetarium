import React, { useState } from "react";

import "./index.css";
import dashboardReducer, { initialState } from "../PlotState/dashboardReducer";
import { DashboardProvider } from "../PlotState/dashboardState";
import TextField from "@mui/material/TextField";
import Autocomplete from '@mui/material/Autocomplete';

import Radio from "@mui/material/Radio";
import RadioGroup from "@mui/material/RadioGroup";
import FormControlLabel from "@mui/material/FormControlLabel";
import FormControl from "@mui/material/FormControl";
import FormLabel from "@mui/material/FormLabel";

import makeStyles from '@mui/styles/makeStyles';
import Typography from "@mui/material/Typography";

import Grid from "@mui/material/Grid";
import PackingCircles from "./PackingCircles";

import { theme } from "../theme/theme.js";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

import { matchSorter } from "match-sorter";

const getDataByKey = (data, key) => [
  ...new Set(data.map(row => row[key]).flat(1))
];

const App = ({ data }) => {
  const [selected, setSelected] = useState({});
  const [modifiedData, setModifiedData] = useState([...data]);

  const handleFilterChange = (data, value, type) => {
    setSelected({ ...selected, [type]: value });
  };
  const filter = (data, keys, { inputValue }) => {
    return matchSorter(data, inputValue, {
      keys: [...keys]
    });
  };

  const filterOptions = (options, params, key) => {
    if (params.inputValue === null) {
      const newSelected = { ...selected, [key]: null };
      setSelected({ ...newSelected });
      const searchParams = Object.keys(newSelected)
        .map(key => newSelected[key])
        .filter(key => key !== null);
      if (searchParams.length === 0) {
        setModifiedData([...data]);
      } else {
        const filtered = filter(data, Object.keys(newSelected), searchParams);
        if (filtered.length !== modifiedData.length) {
          setModifiedData([...filtered]);
        }
      }
    } else if (params.inputValue !== "") {
      const newSelected = { ...selected, [key]: params.inputValue };
      setSelected({ ...newSelected });

      const filtered = filter(modifiedData, [key], params);

      if (filtered.length !== modifiedData.length) {
        setModifiedData([...filtered]);
      }
      return getDataByKey(filtered, key);
    } else {
      return getDataByKey(modifiedData, key);
    }
  };

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
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
                <Typography variant="h5" component="h2">
                  Filter:
                </Typography>

                <Search
                  data={[...getDataByKey(modifiedData, "jira_ticket")]}
                  selectedOption={selected["jira_ticket"] || null}
                  filterOptions={filterOptions}
                  type="jira_ticket"
                  title="Analysis Ticket"
                  selectOption={option =>
                    handleFilterChange(modifiedData, option, "jira_ticket")
                  }
                />
                <Search
                  data={[...getDataByKey(modifiedData, "pathology_disease_name")]}
                  filterOptions={filterOptions}
                  type="pathology_disease_name"
                  title="Tumour Type"
                  selectedOption={selected["pathology_disease_name"] || null}
                  selectOption={option =>
                    handleFilterChange(
                      modifiedData,
                      option,
                      "pathology_disease_name"
                    )
                  }
                />
                <Search
                  data={[...getDataByKey(modifiedData, "pool_id")]}
                  filterOptions={filterOptions}
                  type="pool_id"
                  title="Library"
                  selectedOption={selected["pool_id"] || null}
                  selectOption={option =>
                    handleFilterChange(modifiedData, option, "pool_id")
                  }
                />
                <Search
                  data={[
                    ...getDataByKey(modifiedData, "additional_pathology_info")
                  ]}
                  filterOptions={filterOptions}
                  type="additional_pathology_info"
                  title="Subtype"
                  selectedOption={selected["additional_pathology_info"] || null}
                  selectOption={option =>
                    handleFilterChange(
                      modifiedData,
                      option,
                      "additional_pathology_info"
                    )
                  }
                />
                <RadioOptions
                  options={[...getDataByKey(data, "taxonomy_id")]}
                  filterOptions={filterOptions}
                  type="taxonomy_id"
                  title="Taxonomy"
                  selectOption={option => {
                    return filterOptions(
                      modifiedData,
                      { inputValue: option },
                      "taxonomy_id"
                    );
                  }}
                />
              </Grid>
              <PackingCircles
                modifiedData={modifiedData}
                chartDim={{
                  height: 800,
                  width: 950
                }}
              />
            </Grid>
          </Grid>
        </DashboardProvider>
      </ThemeProvider>
    </StyledEngineProvider>
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
const RadioOptions = ({ options, title, selectOption }) => {
  const [value, setValue] = React.useState(null);

  const handleChange = event => {
    if (event.target.value === value) {
      selectOption(null);
      setValue("");
    } else {
      selectOption(event.target.value);
      setValue(event.target.value);
    }
  };
  return (
    <FormControl component="fieldset">
      <FormLabel component="legend">{title}</FormLabel>
      <RadioGroup key={title + "-radio"} value={value}>
        {options.map(option => (
          <FormControlLabel
            value={option}
            control={<Radio onClick={handleChange} />}
            label={option}
          />
        ))}
      </RadioGroup>
    </FormControl>
  );
};
const Search = ({
  data,
  selectOption,
  selectedOption,
  title,
  type,
  filterOptions
}) => {
  const classes = useStyles();

  return (
    <Autocomplete
      classes={classes}
      options={data}
      value={selectedOption}
      getOptionLabel={option => option}
      style={{ width: 300 }}
      renderOption={option => option}
      onChange={(event, option) => {
        filterOptions(data, { inputValue: option }, type);
      }}
      filterOptions={(options, params) => filterOptions(options, params, type)}
      renderInput={params => (
        <TextField
          {...params}
          InputLabelProps={{
            style: { color: "#black" }
          }}
          label={title}
          variant="outlined"
        />
      )}
    />
  );
};
export default App;
