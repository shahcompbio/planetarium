import React, { useEffect, useState } from "react";
import _ from "lodash";

import * as d3 from "d3";

import axios from "axios";

import Typography from "@material-ui/core/Typography";
import Table from "@material-ui/core/Table";
import TableBody from "@material-ui/core/TableBody";
import TableCell from "@material-ui/core/TableCell";
import TableRow from "@material-ui/core/TableRow";

const TOP_NUM = 10;

const TopGenesPanel = ({ ids, dashboardID, api }) => {
  // Top 10 genes of cells (specified by ids) by average expression
  const [genes, setGenes] = useState([]);
  useEffect(() => {
    const query = {
      size: 0,
      query: {
        bool: {
          filter: {
            terms: {
              cell_id: ids,
            },
          },
        },
      },
      aggs: {
        top_genes: {
          terms: {
            field: "gene",
            size: 10,
            order: { "expression_stats.avg": "desc" },
          },
          aggs: {
            expression_stats: { stats: { field: "expression" } },
          },
        },
      },
    };
    axios
      .post(`${api}/${dashboardID.toLowerCase()}/_search`, query)
      .then((res) => {
        const results = res.data.aggregations.top_genes.buckets;
        const data = results.map((record) => ({
          gene: record["key"],
          avg_expression: record["expression_stats"]["avg"],
        }));
        setGenes(data);
      });
  }, [ids.length, api]);

  return (
    <div>
      <Typography color="textSecondary" gutterBottom>
        Top {TOP_NUM} expressed (average) genes
      </Typography>
      <Table size="small">
        <TableBody>
          {genes.map((gene) => (
            <TableRow key={gene["gene"]} style={{ paddingBottom: "5px" }}>
              <TableCell style={{ borderBottom: "none" }}>
                <b>{gene["gene"]}</b>
              </TableCell>
              <TableCell style={{ borderBottom: "none" }}>
                {d3.format("2s")(gene["avg_expression"])}
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
};

export default TopGenesPanel;
