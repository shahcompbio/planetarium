import React from "react";

import * as d3 from "d3";

import {
  Box,
  CircularProgress,
  Typography,
  Table,
  TableBody,
  TableHead,
  TableCell,
  TableRow,
} from "@mui/material";

const TtestResults = ({ data, type, isLoading, count }) => {
  return count ? (
    <div style={{ height: 400, overflowY: "scroll", overflowX: "none" }}>
      <Typography color="textSecondary" gutterBottom>
        <b>{count}</b> cells selected
      </Typography>
      <Typography color="textSecondary" gutterBottom>
        Top genes
      </Typography>
      {count > 0 && data.length === 0 ? (
        <Box
          sx={{
            display: "flex",
            height: 250,
            width: 150,
            paddingTop: 25,
            justifyContent: "space-evenly",
          }}
        >
          <CircularProgress />
        </Box>
      ) : (
        <Table size="small" style={{ height: 300 }}>
          <TableHead>
            <TableRow>
              <TableCell style={{ fontSize: "0.7rem" }}>Gene</TableCell>
              <TableCell style={{ fontSize: "0.7rem" }} align="right">
                Adjusted P value
              </TableCell>
            </TableRow>
          </TableHead>
          <TableBody>
            {data.map((gene) => (
              <TableRow key={gene["gene"]} style={{ paddingBottom: "5px" }}>
                <TableCell style={{ borderBottom: "none" }}>
                  <b>{gene["gene"]}</b>
                </TableCell>
                <TableCell style={{ borderBottom: "none" }}>
                  {d3.format("0.3")(gene["p"])}
                </TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      )}
    </div>
  ) : (
    <div />
  );
};

export default TtestResults;
