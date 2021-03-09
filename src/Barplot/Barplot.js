import React, { useEffect } from "react";
import * as d3 from "d3";
import * as d3Array from "d3-array";

import { useDashboardState } from "../PlotState/dashboardState";

import { canvasInit, drawAxis } from "../DrawingUtils/utils.js";

const Umap = ({ data, chartDim }) => {
  const [
    { xParam, yParam, cellIdParam, clonotypeParam, topTen, colors }
  ] = useDashboardState();

  useEffect(() => {
    drawAll(data, chartDim);
  }, [data]);
}
function drawAll(data, chartDim){

}

  return (
    <div>
      <div style={{ width: 600, height: 700, position: "relative" }}>
        <div
          id="barplot"
          style={{
            position: "absolute",
            pointerEvents: "all",
            display: "flex"
          }}
        >
          <canvas id="barplotCanvas" />
        </div>
      </div>
    </div>
  );
};
export default Umap;
