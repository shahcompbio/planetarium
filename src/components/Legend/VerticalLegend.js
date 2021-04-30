import React from "react";
import * as d3 from "d3";

import { useD3 } from "../utils/useD3";

const LEGEND_SQUARE_LENGTH = 10;
const LEGEND_SQUARE_SPACING = 8;

const VerticalLegend = ({ width, height, labels, setHighlighted }) => {
  const mouseEvents = (element) =>
    element
      .on("mouseenter", function(d) {
        d3.event.stopPropagation();
        setHighlighted("mouseenter", d["value"]);
      })
      .on("mousedown", function(d, i) {
        d3.event.stopPropagation();
        setHighlighted("mousedown", d["value"]);
      })
      .on("mouseout", function(d, i) {
        d3.event.stopPropagation();
        setHighlighted("mouseout", d["value"]);
      });

  const svgRef = useD3(
    (svg) => {
      const subsets = svg
        .selectAll("g")
        .data(labels)
        .enter()
        .append("g")
        .attr("cursor", "pointer")
        .call(mouseEvents);

      subsets
        .append("rect")
        .attr("width", LEGEND_SQUARE_LENGTH)
        .attr("height", LEGEND_SQUARE_LENGTH)
        .attr("x", 5)
        .attr(
          "y",
          (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5
        )
        .attr("fill", (d) => d["color"]);

      subsets
        .append("text")
        .attr("alignment-baseline", "hanging")
        .attr("dominant-baseline", "hanging")
        .attr("text-align", "left")
        .attr("font-family", "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr("x", LEGEND_SQUARE_LENGTH + 10)
        .attr(
          "y",
          (d, i) => i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) + 5
        )
        .text((d) => d["label"]);
    },
    width,
    height,
    []
  );

  return <svg ref={svgRef} />;
};

export default VerticalLegend;
