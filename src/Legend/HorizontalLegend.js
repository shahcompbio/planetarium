import React from "react";
import * as d3 from "d3";

import { useD3 } from "../utils/useD3";

const LEGEND_SQUARE_LENGTH = 12;
const LEGEND_SQUARE_SPACING = 8;
const PADDING_SIDE = 10;
const PADDING_TOP = 20;

const HorizontalLegend = ({ labels, title, fontFamily = null }) => {
  const legendWidth =
    labels.length * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) +
    PADDING_SIDE * 2;
  const legendHeight = PADDING_TOP + LEGEND_SQUARE_LENGTH * 2;

  //   const mouseEvents = (element) =>
  //     element
  //       .on("mouseenter", function(d) {
  //         d3.event.stopPropagation();
  //         setHighlighted("mouseenter", d["value"]);
  //       })
  //       .on("mousedown", function(d, i) {
  //         d3.event.stopPropagation();
  //         setHighlighted("mousedown", d["value"]);
  //       })
  //       .on("mouseout", function(d, i) {
  //         d3.event.stopPropagation();
  //         setHighlighted("mouseout", d["value"]);
  //       });

  const svgRef = useD3(
    (svg) => {
      const titleText = svg
        .append("text")
        .attr("alignment-baseline", "center")
        .attr("dominant-baseline", "center")
        .attr("text-anchor", "left")
        .attr("font-family", fontFamily ? fontFamily["regular"] : "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr("x", PADDING_SIDE)
        .attr("y", (1 * legendHeight) / 3)
        .text(title);

      const titleWidth = PADDING_SIDE + titleText.node().getBBox().width;
      svg.attr("width", legendWidth + titleWidth);
      // console.log(titleText.node().getBBox());

      const subsets = svg
        .selectAll("g")
        .data(labels)
        .enter()
        .append("g")
        .attr("cursor", "pointer");
      // .call(mouseEvents);

      subsets
        .append("rect")
        .attr("width", LEGEND_SQUARE_LENGTH)
        .attr("height", LEGEND_SQUARE_LENGTH)
        .attr(
          "x",
          (d, i) =>
            titleWidth +
            PADDING_SIDE +
            i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING)
        )
        .attr("y", 5)
        .attr("fill", (d) => d["color"]);

      subsets
        .append("text")
        .attr("alignment-baseline", "hanging")
        .attr("dominant-baseline", "hanging")
        .attr("text-anchor", "middle")
        .attr("font-family", fontFamily ? fontFamily["regular"] : "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr(
          "x",
          (d, i) =>
            titleWidth +
            PADDING_SIDE +
            i * (LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING) +
            LEGEND_SQUARE_LENGTH / 2
        )
        .attr("y", LEGEND_SQUARE_LENGTH + 8)
        .text((d) => d["label"]);
    },
    legendWidth,
    legendHeight,
    []
  );

  return <svg ref={svgRef} />;
};

export default HorizontalLegend;
