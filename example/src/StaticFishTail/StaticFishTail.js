import React, { useState, useEffect } from "react";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";
import { theme } from "../theme/theme.js";
import { MuiThemeProvider } from "@material-ui/core/styles";
import { CssBaseline } from "@material-ui/core";
import { Fishtail, Sankey } from "@shahlab/planetarium";
import _ from "lodash";
import * as d3 from "d3";
/*      {_.groupBy(data, "patient").map((d) => {
        console.log(d);
        return <StaticFigures data={d} cloneColor={parseCloneColor} timePointOrder={d.map(r=>r[""])}/>;
      })}*/
const DataWrapper = ({ data, cloneColor }) => {
  const [fileName, setFileName] = useState("steve_02.h5ad");
  console.log(cloneColor);

  const cColor = _.groupBy(Object.keys(cloneColor), function (o) {
    return o.split("_")[0];
  });
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
      <div id="tool"></div>
      <StaticFigures data={data} cloneColor={parseCloneColor} cColor={cColor} />
    </div>
  );
};
const StaticFigures = ({ data, cloneColor, cColor }) => {
  const timepoints = [...new Set(data.map((item) => item.timepoint))];
  const subsetValues = ["AE_0", "AE_1", "AE_2", "AE_3", "AE_4", "AE_5"];

  console.log(_.groupBy(data, "patient"));
  return (
    <MuiThemeProvider theme={theme}>
      <CssBaseline />
      <Grid container direction="row">
        <Block>
          {data ? (
            [
              ..._(data)
                .groupBy("patient")
                .map((d) => {
                  const patientNUm = d[0]["patient"];
                  console.log(cColor);
                  if (patientNUm === "P2") {
                    console.log(patientNUm);
                    const thisColorScale = cColor["02"];
                    console.log(thisColorScale);
                    const colorScale = d3
                      .scaleSequential(d3.interpolateRdYlBu)
                      .domain([1.42, -1.0527498]);
                    console.log("hello");
                    console.log(
                      d3
                        .extent(thisColorScale.map((d) => cloneColor[d]))
                        .reverse()
                    );
                    //[1.451, 0]
                    //
                    /*  d3
                          .extent(thisColorScale.map((d) => cloneColor[d]))
                          .reverse()*/
                    //2 - 0.92251694, -1.0527498
                    //3 - 1.6665053, -0.6125342
                    //4 - 1.4224248, -0.2644508
                    console.log(colorScale.domain());
                    const pickColorScale = _.pick(cloneColor, thisColorScale);
                    const ordering = Object.keys(pickColorScale)
                      .filter((c) => c.indexOf("B") !== -1)
                      .sort((a, b) => cloneColor[b] - cloneColor[a]);
                    console.log(ordering);
                    console.log(ordering.map((c) => cloneColor[c]));
                    console.log(cloneColor);
                    console.log(pickColorScale);
                    console.log(d);

                    return (
                      <Fishtail
                        key={patientNUm + "fish"}
                        width={600}
                        height={800}
                        data={d}
                        subsetParam={"clone"}
                        cloneParam={"clone"}
                        timepointOrder={timepoints}
                        colorScale={colorScale}
                        cloneColor={pickColorScale}
                        legendSorting={ordering}
                        interpolateColor={true}
                        normalize={false}
                        timepointOrder={["B", "F"]}
                        addTwoTimepointCurve={true}
                        timepointParam={"timepoint"}
                        //    subsetParam="clone"
                        //    treatment="timepoint"
                      />
                    );
                  } else {
                    return <span />;
                  }
                }),
            ]
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
