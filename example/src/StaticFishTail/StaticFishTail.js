import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import { Fishtail, Sankey } from "@shahlab/planetarium";
import _ from "lodash";
import * as d3 from "d3";

const DataWrapper = ({ data, cloneColor }) => {
  const [fileName, setFileName] = useState("steve_02.h5ad");
  const parseCloneColor = Object.keys(cloneColor).reduce((final, d) => {
    final[d] = parseFloat(cloneColor[d]);
    return final;
  }, {});
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
                fetch("http://localhost:5000/render2/" + input.value + "/")
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

      <StaticFigures data={data} cloneColor={parseCloneColor} />
    </div>
  );
};
const StaticFigures = ({ data, cloneColor }) => {
  const timepoints = [...new Set(data.map((item) => item.timepoint))];
  const subsetValues = ["AE_0", "AE_1", "AE_2", "AE_3", "AE_4", "AE_5"];

  const colorScale = d3
    .scaleSequential(d3.interpolateRdYlBu)
    .domain(
      d3.extent(Object.keys(cloneColor).map((d) => cloneColor[d])).reverse()
    );

  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          {data ? (
            <Fishtail
              width={600}
              height={600}
              data={data}
              subsetParam={"clone"}
              cloneParam={"clone"}
              timepointOrder={timepoints}
              colorScale={colorScale}
              cloneColor={cloneColor}
              legendSorting={timepoints}
              interpolateColor={true}
              timepointOrder={["EP-02-B", "EP-02-F"]}
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
