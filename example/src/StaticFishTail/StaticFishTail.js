import React, { useState } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import { Fishtail, Sankey } from "@shahlab/planetarium";
import _ from "lodash";
import * as d3 from "d3";
/*const cloneColor = {
  L: 0.294,
  O: 0.4468,
};*/
/*const cloneColor = {
  "AE_3-pre": 0.29674882613378417,
  "AE_3-post": 0.23518855684429804,
  "AE_1-pre": 0.1592361056273969,
  "AE_1-post": 0.16845969969976557,
  "AE_2-pre": -0.133290918779613,
  "AE_2-post": -0.05529150992486014,
  "AE_5-pre": -0.13487067178536521,
  "AE_5-post": -0.1748878242053734,
  "AE_0-pre": 0.0012924036369184586,
  "AE_0-post": -0.08845248144109227,
  "AE_4-pre": 0.22772281846699516,
  "AE_4-post": 0.27497373627812005,
};*/

/*const cloneColor = {
  M: 0.10615526,
  I: 0,
  Q: 0.63132566,
  G: 0.018968515,
  E: 0,
};*/
/*const cloneColor = {
  J: 0.057269838,
  I: 0,
  V: 0.03315382,
  F: 0,
  Y: 0.04308413,
  G: 0.0868441,
  P: 0.09321452,
  X: 0,
};*/
const DataWrapper = ({ data, cloneColor }) => {
  const [fileName, setFileName] = useState("steve_02.h5ad");
  console.log(cloneColor);
  console.log(data);
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
  console.log(data);
  const timepoints = [...new Set(data.map((item) => item.timepoint))];
  //const COLOR_ARRAY = ["rgb(50, 136, 189)", "#79b9e0"];
  const subsetValues = ["AE_0", "AE_1", "AE_2", "AE_3", "AE_4", "AE_5"];

  //  const subsetValues = ["I", "E", "M", "G", "Q"];
  //const COLOR_ARRAY = ["#8f0b03", "#a1054d", "#E1341E", "#fc479b", "#fc9b95"];

  //  const subsetValues = ["I", "F", "X", "V", "Y", "J", "G", "P"];
  /*  const COLOR_ARRAY = [
    "#5f6e03",
    "#156e03",
    "#074030",
    "#1bbf93",
    "#36bf3a",
    "#9ec965",
    "#e4fc77",
    "#a6ffc9",
  ];*/
  /*  const colorScale = d3
    .scaleOrdinal(d3.quantize(d3.interpolateViridis, subsetValues.length))
    .domain(subsetValues);
*/
  console.log(
    d3.extent(Object.keys(cloneColor).map((d) => cloneColor[d])).reverse()
  );
  const colorScale = d3
    .scaleSequential(d3.interpolateRdYlBu)
    //.scaleSequential(d3.interpolateViridis)
    .domain(
      d3.extent(Object.keys(cloneColor).map((d) => cloneColor[d])).reverse()
    );
  //  const colorScale = d3
  //    .scaleOrdinal(d3.interpolateViridis)
  //    .domain(subsetValues);

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
              timepointOrder={["EP-04-B", "EP-04-F"]}
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
