import React, { useState } from "react";
import Select from "./components/Select";
import UMAP from "./components/UMAP";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";

const DataWrapper = ({ data }) => (
  <CITESEQ
    cells={data["cells"]}
    genes={data["genes"]}
    proteins={data["proteins"]}
  />
);
const CITESEQ = ({ cells, genes, proteins }) => {
  const [gene, setGene] = useState(genes[0]);
  const [protein, setProtein] = useState(proteins[0]);
  const [highlighted, setHighlighted] = useState(null);
  const [activeGraph, setActiveGraph] = useState(null);
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          <Grid container direction="column" alignItems="center">
            <Grid item>
              <Select
                options={genes}
                title={"Gene"}
                value={gene}
                onSelect={setGene}
              />
            </Grid>
            <Grid item>
              <UMAP
                width={500}
                height={500}
                data={cells}
                xParam={"rna_UMAP_1"}
                yParam={"rna_UMAP_2"}
                subsetParam={gene}
                idParam={"cell_id"}
                onLasso={(value) => {
                  const subset = value.length === cells.length ? null : value;
                  setHighlighted(subset);
                  setActiveGraph(subset === null ? null : "gene_umap");
                }}
                highlighted={activeGraph !== "gene_umap" ? highlighted : null}
              />
            </Grid>
          </Grid>
        </Block>

        <Block>
          <Grid container direction="column" alignItems="center">
            <Grid item>
              <Select
                options={proteins}
                title={"Protein"}
                value={protein}
                onSelect={setProtein}
              />
            </Grid>
            <Grid item>
              <UMAP
                width={500}
                height={500}
                data={cells}
                xParam={"protein_UMAP_1"}
                yParam={"protein_UMAP_2"}
                subsetParam={protein}
                idParam={"cell_id"}
                onLasso={(value) => {
                  const subset = value.length === cells.length ? null : value;
                  setHighlighted(subset);
                  setActiveGraph(subset === null ? null : "protein_umap");
                }}
                highlighted={
                  activeGraph !== "protein_umap" ? highlighted : null
                }
              />
            </Grid>
          </Grid>
        </Block>
      </Grid>
    </MuiThemeProvider>
  );
};

const Block = ({ children }) => (
  <Grid item>
    <Paper
      style={{
        margin: 10,
        padding: 10,
      }}
    >
      {children}
    </Paper>
  </Grid>
);

export default DataWrapper;
