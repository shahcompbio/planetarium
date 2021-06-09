import React from "react";
import * as d3 from "d3";

import { useD3 } from "../utils/useD3";

const LEGEND_SQUARE_LENGTH = 10;
const LEGEND_SQUARE_SPACING = 8;

const STEP = LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING;

const PADDING = 20;

const VerticalLegend = ({ width, labels, title, setHighlighted }) => {
  const legendWidth = width;
  const legendHeight = PADDING * 2 + labels.length * STEP;

  const svgRef = useD3(
    (svg) => {
      let offset = 5;
      if (title !== undefined) {
        const titleText = svg
          .append("text")
          .attr("alignment-baseline", "center")
          .attr("dominant-baseline", "center")
          .attr("text-anchor", "middle")
          .attr("font-family", "Helvetica")
          .attr("font-weight", "500")
          .attr("font-size", "12px")
          .attr("fill", "#000000")
          .attr("x", legendWidth / 2)
          .attr("y", PADDING)
          .text(title);

        const titleHeight = PADDING + titleText.node().getBBox().height;
        svg.attr("height", legendHeight + titleHeight);
        offset += titleHeight;
      }

      const subsets = svg
        .selectAll("g")
        .data(labels)
        .enter()
        .append("g")
        .attr("cursor", "pointer");

      subsets
        .append("rect")
        .attr("width", LEGEND_SQUARE_LENGTH)
        .attr("height", LEGEND_SQUARE_LENGTH)
        .attr("x", 5)
        .attr("y", (d, i) => i * STEP + offset)
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
        .attr("y", (d, i) => i * STEP + offset)
        .text((d) => d["label"]);

      const findLabelValue = (mouseY) => {
        const index = Math.round((mouseY - offset) / STEP) - 1;

        if (0 <= index && index < labels.length) {
          return labels[index]["value"];
        }

        return null;
      };

      const mousemove = () => {
        const mouseY = d3.event.clientY;

        const label = findLabelValue(mouseY);

        setHighlighted(label === null ? "mouseout" : "mouseenter", label);
      };

      const mouseout = () => {
        setHighlighted("mouseout", null);
      };

      const click = () => {
        const mouseY = d3.event.clientY;

        const label = findLabelValue(mouseY);

        setHighlighted("mousedown", label);
      };

      svg
        .on("mousemove", mousemove)
        .on("mouseout", mouseout)
        .on("click", click);
    },
    legendWidth,
    legendHeight,
    []
  );

  return <svg ref={svgRef} />;
};

export default VerticalLegend;
