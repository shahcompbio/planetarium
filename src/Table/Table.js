import React, { useState, useMemo } from "react";
import * as d3 from "d3";
import { useDashboardState } from "../PlotState/dashboardState";
import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";

import Button from "@material-ui/core/Button";
import IconButton from "@material-ui/core/IconButton";
import ClearIcon from "@material-ui/icons/Clear";
import TextField from "@material-ui/core/TextField";
import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import GetAppIcon from "@material-ui/icons/GetApp";
import * as d3Dsv from "d3-dsv";
import DataTable from "react-data-table-component";
const formatCols = ["adj_pval", "log_fc"];
const formatDecimal = d3.format(",.4f");
const Table = ({ chartName, data, chartDim, selectedSubtype }) => {
  const [filterText, setFilterText] = useState("");
  const [{ subtypeParam }] = useDashboardState();

  const columns = Object.keys(data[0]);
  const dataSource = selectedSubtype
    ? data.filter(row => row[subtypeParam] === selectedSubtype)
    : data;
  const filteredItems = dataSource.filter(
    item =>
      item["gene"] &&
      item["gene"].toLowerCase().includes(filterText.toLowerCase())
  );

  const subHeaderComponentMemo = useMemo(() => {
    const handleClear = () => {
      if (filterText) {
        setFilterText("");
      }
    };

    return (
      <FilterComponent
        onFilter={e => setFilterText(e.target.value)}
        onClear={handleClear}
        filterText={filterText}
        data={data}
      />
    );
  }, [filterText]);
  return (
    <Paper
      style={{
        margin: 10,
        height: chartDim["height"],
        width: chartDim["width"],
        padding: 15
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="flex-start"
        style={{
          width: "100%"
        }}
      >
        <Grid
          item
          style={{
            width: "100%",
            textAlign: "left",
            paddingTop: 15
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
        <Grid
          item
          style={{
            overflowY: "hidden",
            height: chartDim["height"] - 70,
            width: "100%"
          }}
        >
          <DataTable
            subHeader
            fixedHeader
            dense
            noHeader
            defaultSortAsc
            overflowY
            subHeaderComponent={subHeaderComponentMemo}
            compact
            columns={columns.map(col => {
              return formatCols.indexOf(col) !== -1
                ? {
                    name: col,
                    selector: col,
                    sortable: true,
                    right: true,
                    cell: row => (
                      <span>{formatDecimal(parseFloat(row[col]))}</span>
                    )
                  }
                : {
                    name: col,
                    selector: col,
                    sortable: true,
                    right: true
                  };
            })}
            data={filteredItems}
          />
        </Grid>
      </Grid>
    </Paper>
  );
};

const FilterComponent = ({ filterText, onFilter, onClear, data }) => (
  <div style={{ display: "flex" }}>
    <TextField
      color="primary"
      type="text"
      placeholder="Gene"
      aria-label="Gene"
      id="searchGenes"
      placeholder="Filter By Gene"
      aria-label="Search Input"
      value={filterText}
      onChange={onFilter}
      InputProps={{
        endAdornment: (
          <Button
            label="Clear"
            color="primary"
            variant="outlined"
            onClick={onClear}
            style={{ marginLeft: 15, marginBottom: 5 }}
          >
            <ClearIcon />
          </Button>
        )
      }}
    />
    <Button
      variant="outlined"
      label="Download"
      id="tsv-download"
      onClick={() => {
        const dataSource = new Blob([d3Dsv.tsvFormat(data)], {
          type: "text/tsv"
        });
        const tsvURL = window.URL.createObjectURL(dataSource);
        const tempLink = document.createElement("a");
        tempLink.href = tsvURL;
        tempLink.setAttribute("download", "filename.tsv");
        tempLink.click();
      }}
      color="secondary"
      style={{ marginLeft: 15 }}
    >
      <GetAppIcon />
    </Button>
  </div>
);
export default Table;
