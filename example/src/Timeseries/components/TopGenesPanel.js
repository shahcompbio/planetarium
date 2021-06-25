import React from "react";
import _ from "lodash";

import * as d3 from "d3";

import Typography from "@material-ui/core/Typography";
import Table from "@material-ui/core/Table";
import TableBody from "@material-ui/core/TableBody";
import TableCell from "@material-ui/core/TableCell";
import TableRow from "@material-ui/core/TableRow";

const TOP_NUM = 10;

const TopGenesPanel = ({ ids, genes }) => {
  // Top 10 genes of cells (specified by ids) by average expression

  const filteredGenes = genes.filter((gene) => ids.includes(gene["cell_id"]));

  const groupedGenes = _.groupBy(filteredGenes, (datum) => datum["gene"]);
  const geneAverages = Object.keys(groupedGenes).map((gene) => ({
    gene,
    average: _.sumBy(groupedGenes[gene], "exp") / groupedGenes[gene].length,
  }));

  const topGenes = geneAverages
    .sort((a, b) => b["average"] - a["average"])
    .slice(0, TOP_NUM);

  return (
    <div>
      <Typography color="textSecondary" gutterBottom>
        Top {TOP_NUM} expressed (average) genes
      </Typography>
      <Table size="small">
        <TableBody>
          {topGenes.map((gene) => (
            <TableRow style={{ paddingBottom: "5px" }}>
              <TableCell style={{ borderBottom: "none" }}>
                <b>{gene["gene"]}</b>
              </TableCell>
              <TableCell style={{ borderBottom: "none" }}>
                {d3.format("2s")(gene["average"])}
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
};

export default TopGenesPanel;
