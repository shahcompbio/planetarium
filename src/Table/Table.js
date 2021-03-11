import React, { useState, useMemo, useCallback } from "react";

import * as d3 from "d3";
import * as d3Array from "d3-array";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";

const Table = ({ data, chartDim, selectedSubtype }) => {
  const [
    { xParam, yParam, cellIdParam, clonotypeParam, topTen, subtypeParam }
  ] = useDashboardState();

  const [selectedRows, setSelectedRows] = useState(null);
  const { columns } = data;
  const dataSource = selectedSubtype
    ? data.filter(row => row[subtypeParam] === selectedSubtype)
    : data;
  console.log(selectedSubtype);
  return (
    <div style={{ height: chartDim["height"], overflow: "auto" }}>
      <table style={{ height: chartDim["height"] }}>
        <tr>
          {columns.map(column => (
            <th>{column}</th>
          ))}
        </tr>
        {dataSource.map(row => {
          return (
            <tr>
              {columns.map(column => {
                return <td>{row[column]}</td>;
              })}
            </tr>
          );
        })}
      </table>
    </div>
  );
};
export default Table;
