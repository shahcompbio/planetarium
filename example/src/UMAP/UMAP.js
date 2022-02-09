import React, { useState, useEffect } from "react";
import _ from "lodash";
import * as d3 from "d3";

import TTestUmap from "./TTestUmap";

import { Select } from "@shahlab/planetarium";

import Grid from "@mui/material/Grid";

import { theme } from "./theme";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import CssBaseline from "@mui/material/CssBaseline";

import Paper from "@mui/material/Paper";
import MetaData from "./MetaData";
import Filters from "./Filters";

import { CONSTANTS } from "./config";

const PHENOTYPE_COLORS = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];
const DataWrapper = ({ data }) => (
  <VDJ
    metadata={data["metadata"]}
    filters={data["filters"]}
    degs={data["degs"]}
  />
);
//const url = "https://spectrum-staging.shahlab.mskcc.org";
const url = process.env.HOST ? process.env.HOST : "http://localhost:5000";
export const VDJ = ({ metadata, degs, filters }) => {
  const { clonotypeParam, subtypeParam, xParam, yParam } = CONSTANTS;

  const [selectClone, setSelectClone] = useState(null);
  const [selectIDs, setSelectIDs] = useState(null);
  const [activeGraph, setActiveGraph] = useState(null);

  const [subset, setSubset] = useState(subtypeParam);
  const [selectSubset, setSelectSubset] = useState(null);

  const [selectFilters, setSelectFilters] = useState(null);

  const [tTestData, settTestData] = useState([]);
  const data =
    selectFilters === null
      ? metadata
      : metadata.filter(
          (datum) => datum[selectFilters[0]] === selectFilters[1]
        );

  useEffect(() => {
    if (selectIDs !== null) {
      if (selectIDs.length !== 0) {
        const param = selectIDs.join(",");
        fetch(url + "/ttestSantosh/", {
          method: "POST",
          credentials: "include",
          body: JSON.stringify({ data: param }),
        })
          .then((res) => res.json())
          .then((result) => {
            if (result.data) {
              settTestData(result.data);
            }
          });
      }
    } else if (selectSubset !== null) {
      const ids = data
        .filter((datum) => datum[subset] === selectSubset)
        .map((datum) => datum["cell_id"])
        .join(",");

      fetch(url + "/ttestSantosh/", {
        method: "POST",
        credentials: "include",
        body: JSON.stringify({ data: ids }),
      })
        .then((res) => res.json())
        .then((result) => {
          if (result.data) {
            settTestData(result.data);
          }
        });
    }
  }, [selectIDs, selectSubset]);

  const phenotypeValues = Object.keys(_.groupBy(data, subset)).sort();
  const phenotypeColorScale = d3
    .scaleOrdinal()
    .domain(phenotypeValues)
    .range(
      PHENOTYPE_COLORS.slice(
        0,
        Math.min(phenotypeValues.length, PHENOTYPE_COLORS.length)
      )
    )
    .unknown("#e8e8e8");

  const highlightData =
    selectIDs !== null
      ? data.filter((datum) => selectIDs.includes(datum["cell_id"]))
      : selectClone !== null
      ? data.filter((datum) => datum[clonotypeParam] === selectClone)
      : selectSubset !== null
      ? data.filter((datum) => datum[subset] === selectSubset)
      : null;

  const highlightIDs =
    highlightData === null
      ? null
      : highlightData.map((datum) => datum["cell_id"]);

  return (
    <StyledEngineProvider injectFirst>
      <ThemeProvider theme={theme}>
        <CssBaseline />
        <Grid
          container
          direction="column"
          justifyContent="flex-start"
          alignItems="flex-start"
          style={{
            minWidth: 1600,
            xOverflow: "scroll",
            padding: 15,
            marginBottom: 10,
          }}
        >
          <Grid
            item
            container
            direction="row"
            justifyContent="flex-start"
            alignItems="flex-start"
          >
            <Paper
              elevation={0}
              style={{
                background: "none",
                margin: 10,
                padding: 10,
              }}
            >
              <MetaData
                sample="Santosh"
                hasSelection={
                  highlightData || selectClone || selectSubset || selectFilters
                }
                setHighlight={() => {
                  setSelectIDs(null);
                  setSelectClone(null);
                  setSelectSubset(null);
                  setActiveGraph(null);
                  setSelectFilters(null);
                }}
                totalCount={metadata.length}
              />
              <Filters
                selected={selectFilters}
                filters={filters}
                setFilters={setSelectFilters}
              />
            </Paper>

            <TTestUmap
              width={700}
              height={600}
              data={data}
              xParam={xParam}
              yParam={yParam}
              subsetParam={subset}
              idParam="cell_id"
              colorScale={phenotypeColorScale}
              onLasso={(data) => {
                setSelectIDs(
                  data === null ? null : data.map((datum) => datum["cell_id"])
                );
                setActiveGraph(data === null ? null : "phenoUMAP");
              }}
              onLegendClick={(value) => {
                setSelectSubset(value);
                setActiveGraph(value === null ? null : "phenoUMAP");
              }}
              disable={activeGraph !== null && activeGraph !== "phenoUMAP"}
              highlightIDs={highlightIDs}
              Select={
                <span style={{ marginRight: 10 }}>
                  <Select
                    width={200}
                    title={"Color By"}
                    value={subset}
                    options={filters.map((datum) => datum["name"])}
                    onSelect={setSubset}
                  />
                </span>
              }
              tTestData={tTestData}
            />
          </Grid>
        </Grid>
      </ThemeProvider>
    </StyledEngineProvider>
  );
};

export default DataWrapper;
