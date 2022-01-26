import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import Heatmap from "./components/Heatmap";

const DataWrapper = ({ data }) => {
  const [fileName, setFileName] = useState("renormed_protein.h5ad");
  const [defaultGenes, setDefaultGenes] = useState([
    "HLA-DRB5",
    "HLA-DPB1",
    "HLA-DPA1",
    "HLA-DRB1",
    "HLA-DRA",
    "HLA-DQB1",
    "HLA-DMA",
    "CD74",
    "CIITA",
    "HLA-DQA1",
    "HLA-DQA2",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-DOA",
    "HLA-DMB",
    "HLA-DMA",
    "HLA-G",
    "HLA-DQB1-AS1",
    "HLA-F",
    "B2M",
  ]);

  return (
    <div>
      <div style={{ with: 400 }}>
        <Paper style={{ width: 700, margin: 10, padding: 10 }}>
          <div>
            <div>
              <label style={{ width: "100%" }}>
                Input Genes (comma seperated)
              </label>
            </div>
            <textarea name="paragraph_text" cols="50" rows="5" id="genes-list">
              {defaultGenes.join(",")}
            </textarea>
          </div>
          <div>
            <label style={{ width: "100%" }}>Input File</label>
          </div>
          <input
            type="file"
            id="input"
            type="text"
            value={fileName}
            onChange={(event) => setFileName(event.target.value)}
          />
          <button
            style={{ float: "right", padding: 5, margin: 10, marginTop: -15 }}
            onClick={() => {
              var input = document.getElementById("input");
              if (input !== "") {
                const genes = document.getElementById("genes-list").value;
                fetch("http://localhost:5000/render/" + input.value + "/", {
                  method: "POST",
                  body: JSON.stringify({ data: genes }),
                })
                  .then(function (response) {
                    return response.json();
                  })
                  .then(function (data) {
                    console.log(data);
                  });
              }
            }}
          >
            Render
          </button>
        </Paper>
      </div>

      <StaticFigures
        patientsData={data["patients"]}
        geneData={data["geneData"]}
      />
    </div>
  );
};
const StaticFigures = ({ patientsData, geneData }) => {
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          {patientsData ? (
            <Heatmap data={geneData} patients={patientsData} />
          ) : (
            <span />
          )}
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
