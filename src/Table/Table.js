import React, { useState, useMemo } from "react";
import * as d3 from "d3";
import { useDashboardState } from "../PlotState/dashboardState";
import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";

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
      />
    );
  }, [filterText]);
  return (
    <div
      class="card"
      style={{
        margin: 10,
        height: chartDim["height"],
        width: chartDim["width"],
        padding: 15
      }}
    >
      <DataTable
        title={
          <div
            class="card-title"
            style={{
              width: "100%",
              height: 80,
              textAlign: "left"
            }}
          >
            <h6 class="card-title">
              {infoText[chartName]["title"] + "    "}
              <Info name={chartName} direction="s" />
            </h6>
          </div>
        }
        subHeader
        fixedHeader
        dense
        defaultSortAsc
        overflowY={false}
        subHeaderComponent={subHeaderComponentMemo}
        compact
        columns={columns.map(col => {
          return formatCols.indexOf(col) !== -1
            ? {
                name: col,
                selector: col,
                sortable: true,
                right: true,
                cell: row => <span>{formatDecimal(parseFloat(row[col]))}</span>
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
    </div>
  );
};
const Title = ({ chartName, chartDim }) => (
  <div
    class="card-title"
    style={{
      marginTop: chartDim["chart"]["x1"],
      width: "100%",
      height: 80,
      paddingTop: 40,
      paddingLeft: -50,
      textAlign: "left"
    }}
  >
    {infoText[chartName]["title"] + "    "}

    <Info name={chartName} direction="s" />
  </div>
);
const FilterComponent = ({ filterText, onFilter, onClear }) => (
  <div style={{ display: "flex" }}>
    <input
      type="text"
      class="form-control"
      placeholder="Gene"
      aria-label="Gene"
      aria-describedby="basic-addon1"
      id="searchGenes"
      type="text"
      placeholder="Filter By Gene"
      aria-label="Search Input"
      value={filterText}
      onChange={onFilter}
    />
    <button type="button" class="btn btn-primary" onClick={onClear}>
      X
    </button>
  </div>
);
export default Table;
