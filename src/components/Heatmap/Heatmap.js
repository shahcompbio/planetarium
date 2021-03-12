import React from "react";
import * as d3 from "d3";

const Heatmap = ({ data }) => {
  return (
    <div style={{ width: 600, height: 700, position: "relative" }}>
      <div
        id="heatmap"
        style={{
          position: "absolute",
          pointerEvents: "all",
          display: "flex",
        }}
      >
        <canvas id="heatmapCanvas" />
      </div>
    </div>
  );
};
