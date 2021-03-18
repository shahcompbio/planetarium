import React, { useState, useMemo } from "react";
import * as d3 from "d3";
import { useDashboardState } from "../PlotState/dashboardState";

import DataTable from "react-data-table-component";
const formatCols = ["adj_pval", "log_fc"];
const formatDecimal = d3.format(",.4f");
const Table = ({ data, chartDim, selectedSubtype }) => {
  const [filterText, setFilterText] = useState("");
  const [{ subtypeParam }] = useDashboardState();

  const columns = Object.keys(data[0]);
  const dataSource = selectedSubtype
    ? data.filter((row) => row[subtypeParam] === selectedSubtype)
    : data;
  const filteredItems = dataSource.filter(
    (item) =>
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
        onFilter={(e) => setFilterText(e.target.value)}
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
        overflow: "scroll",
        padding: 15,
      }}
    >
      <DataTable
        title=""
        subHeader
        fixedHeader
        dense
        defaultSortAsc
        overflowY={false}
        subHeaderComponent={subHeaderComponentMemo}
        compact
        columns={columns.map((col) => {
          return formatCols.indexOf(col) !== -1
            ? {
                name: col,
                selector: col,
                sortable: true,
                right: true,
                cell: (row) => (
                  <span>{formatDecimal(parseFloat(row[col]))}</span>
                ),
              }
            : {
                name: col,
                selector: col,
                sortable: true,
                right: true,
              };
        })}
        data={filteredItems}
      />
    </div>
  );
};
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
