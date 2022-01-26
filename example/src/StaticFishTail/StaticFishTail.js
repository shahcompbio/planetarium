import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import { Fishtail, Sankey } from "@shahlab/planetarium";

const DataWrapper = ({ data }) => {
  const [fileName, setFileName] = useState("prog_AE.h5ad");

  return (
    <div>
      <div style={{ with: 400 }}>
        <Paper style={{ width: 700, margin: 10, padding: 10 }}>
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
                fetch("http://localhost:5000/render/" + input.value + "/")
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

      <StaticFigures data={data["data"]} />
    </div>
  );
};
const StaticFigures = ({ data }) => {
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          {data ? (
            <Fishtail
              width={700}
              height={300}
              data={data}
              subsetParam={"clone"}
              cloneParam={"clone"}
              timepointOrder={["Pre", "Post"]}
              addTwoTimepointCurve={true}
              timepointParam={"timepoint"}
              //    subsetParam="clone"
              //    treatment="timepoint"
            />
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
