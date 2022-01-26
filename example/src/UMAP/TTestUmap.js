import React, { useEffect, useState } from "react";
import _ from "lodash";

import { Layout, UMAP } from "@shahlab/planetarium";
import { CircularProgress } from "@mui/material";

import TtestResults from "./tTestResults";

const TTestUmap = ({
  data,
  xParam,
  yParam,
  subsetParam,
  idParam,
  colorScale,
  onLasso,
  onLegendClick,
  disable,
  highlightIDs,
  Select,
  loadingTest,
  tTestData,
}) => {
  return (
    <Layout title={""} infoText={""} addIcon={[Select]}>
      <UMAP
        width={700}
        height={600}
        data={data}
        xParam={xParam}
        yParam={yParam}
        subsetParam={subsetParam}
        idParam={idParam}
        colorScale={colorScale}
        onLasso={onLasso}
        onLegendClick={onLegendClick}
        disable={disable}
        highlightIDs={highlightIDs}
        MoreInfoComponent={() => (
          <TtestResults
            data={tTestData}
            count={highlightIDs ? highlightIDs.length : null}
          />
        )}
      />
    </Layout>
  );
};

export default TTestUmap;
