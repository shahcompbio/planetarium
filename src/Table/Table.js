import React, { useState, useMemo, useCallback } from "react";

import * as d3 from "d3";
import * as d3Array from "d3-array";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";
//import Table from "react-bootstrap/Table";
import DataTable from "react-data-table-component";

const Table = ({ data, chartDim, selectedSubtype }) => {
  const [
    { xParam, yParam, cellIdParam, clonotypeParam, topTen, subtypeParam },
  ] = useDashboardState();

  console.log(data);
  const [selectedRows, setSelectedRows] = useState(null);
  // const { columns } = data;
  const columns = Object.keys(data[0]);
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
