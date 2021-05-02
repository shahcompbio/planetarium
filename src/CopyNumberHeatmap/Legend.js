import React, { useCallback } from "react";

import * as d3 from "d3";
import { colorScale } from "./utils.js";
import { heatmapConfig } from "./config.js";

const Legend = ({ heatmapConfig, maxState }) => {
  function LegendCallback() {
    const setRef = useCallback(node => {
      if (node) {
        const legend = d3.select(node);

        legend
          .append("text")
          .attr("class", "legend title")
          .attr("x", 0)
          .attr("y", heatmapConfig.legend.titleHeight)
          .attr("text-anchor", "start")
          .attr("font-size","x-small")
          .text("Copy Number");

        Array.from(Array(maxState + 1).keys()).forEach((copyNum, i) => {
          var color = colorScale(i);
          const xCord =
            i *
              (heatmapConfig.legend.squareSize +
                heatmapConfig.legend.squareSpacing + 2) +
            heatmapConfig.legend.titleWidth;

          legend
            .append("rect")
            .attr("class", "legend")
            .attr("x", xCord+ 15)
            .attr("y", 6)
            .attr("width", heatmapConfig.legend.squareSize)
            .attr("height", heatmapConfig.legend.squareSize)
            .attr("fill", color);

          legend
            .append("text")
            .attr("class", "legend")
            .attr("x", function() {
              if (i === 10) {
                return xCord - 4 + 15;
              } else if (i === 11) {
                return xCord - 1 + 15;
              } else {
                return xCord + 1 + 15;
              }
            })
            .attr(
              "y",
              heatmapConfig.legend.textOffset + heatmapConfig.legend.squareSize + 4
            )
            .attr("font-size","xx-small")
            .attr("text-anchor", "start")
            .attr("dominant-baseline", "hanging")
            .text(d => (i === maxState ? copyNum + "+" : copyNum));
        });
      }
      ref.current = node;
    }, []);
    return [setRef];
  }
  const [ref] = LegendCallback();

  return (
    <svg
      id="legend"
      ref={ref}
      width={heatmapConfig.legend.width}
      height={heatmapConfig.legend.height}
      style={{ float: "right"}}
    />
  );
};
export default Legend;
