import React from "react";

import { useDashboardState } from "../PlotState/dashboardState";

import DataTable from "react-data-table-component";

const Table = ({ data, chartDim, selectedSubtype }) => {
  const [{ subtypeParam }] = useDashboardState();

  const { columns } = data;
  const dataSource = selectedSubtype
    ? data.filter((row) => row[subtypeParam] === selectedSubtype)
    : data;

  return (
    <div
      style={{
        height: chartDim["height"],
        overflow: "auto",
        marginLeft: 50,
        marginTop: 70,
      }}
    >
      <DataTable
        title=""
        columns={columns.map((col) => ({
          name: col,
          selector: col,
          sortable: true,
          right: true,
        }))}
        data={dataSource}
      />
    </div>
  );
};
export default Table;
