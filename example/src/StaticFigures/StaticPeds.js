import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import HeatmapPeds from "./components/HeatmapPeds";
import _ from "lodash";

const DataWrapper = ({ data }) => {
  console.log(data);
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
  const order = [
    "CD244",
    "CLEC12A",
    "HLA-DRA",
    "GGT1",
    "IL3RA",
    "CD19",
    "CD24",
    "CD81",
    "CD48",
    "MS4A1",
    "TFRC",
    "C5AR1",
    "CD36",
    "DPP4",
    "CD82",
    "CLEC12A",
    "ITGAX",
    "ITGAL",
    "ITGB2",
    "PECAM1",
    "CD2",
    "CD48",
    "ITGAL",
    "B3GAT1",
    "CD3D",
    "CD99",
    "CD69",
    "CD40LG",
    "NECTIN2",
    "TFRC",
    "CD38",
    "CD48",
    "FCGR2A",
    "ENTPD1",
    "ICAM1",
    "CD2",
    "CD5",
    "CD48",
    "CD3D",
    "CD8A",
    "CLEC12A",
    "ITGAX",
    "ITGAL",
    "CD101",
    "ITGB2",
    "IL3RA",
    "CD4",
    "PTPRC",
    "CD36",
    "CD48",
  ];
  const geneOrder = [
    "PRSS57",
    "CYBA",
    "GAPDH",
    "C1QTNF4",
    "ACTG1",
    "CD74",
    "HLA-DRA",
    "CD79A",
    "MS4A1",
    "IGHM",
    "HBB",
    "PRDX2",
    "HBD",
    "BLVRB",
    "UROD",
    "LYZ",
    "S100A8",
    "S100A6",
    "S100A9",
    "FTL",
    "NKG7",
    "GNLY",
    "B2M",
    "KLRD1",
    "GZMA",
    "CDK6",
    "DNAJB6",
    "SOX4",
    "ANKRD28",
    "JMJD1C",
    "MZB1",
    "SSR4",
    "SEC11C",
    "IGKC",
    "IGHG1",
    "IL7R",
    "TNFAIP3",
    "CD3E",
    "CD69",
    "ETS1",
    "CD74",
    "HLA-DRB1",
    "HLA-DRA",
    "HLA-DPA1",
    "HLA-DRB5",
    "TCF4",
    "IRF8",
    "CD74",
    "UGCG",
    "JCHAIN",
  ];
  console.log(data);
  const modifiedData = data.reduce((final, curr) => {
    const geneList = Object.keys(curr).filter((d) => d !== "");
    const modCurr = {
      name: curr[""],
      genes: order.reduce((thisFinal, d) => {
        console.log(curr);
        console.log(curr[d]);
        return [...thisFinal, [d, curr[d]]];
      }, []),
    };
    console.log([...final, { ...modCurr }]);
    return [...final, modCurr];
  }, []);

  console.log(modifiedData);
  const groupedData = _.groupBy(modifiedData, (d) => d["name"]);

  const groups = Object.keys(groupedData);
  return (
    <div>
      <div style={{ with: 400 }}>
        <Paper style={{ width: 1700, margin: 10, padding: 10 }}>
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

      <StaticFigures data={groupedData} groups={groups} />
    </div>
  );
};
const StaticFigures = ({ data, groups }) => {
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          {groups ? <HeatmapPeds data={data} patients={groups} /> : <span />}
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
