import React, { useState } from "react";
import * as d3 from "d3";

import { useD3 } from "../utils/useD3";

const LEGEND_SQUARE_LENGTH = 10;
const LEGEND_SQUARE_SPACING = 8;

const STEP = LEGEND_SQUARE_LENGTH + LEGEND_SQUARE_SPACING;

const PADDING = 20;

const VerticalLegend = ({
  width,
  labels,
  title,
  disable = false,
  onClick = (value) => {},
  onHover = (value) => value,
  fontFamily = null,
}) => {
  const [selected, setSelected] = useState(null);
  const legendWidth = width;
  const legendHeight = PADDING * 2 + labels.length * STEP;

  const svgRef = useD3(
    (svg) => {
      let offset = 6;
      const textOffset = fontFamily ? fontFamily["labelOffset"] : 0;
      if (title !== undefined) {
        const titleText = svg
          .append("text")
          .attr("alignment-baseline", "center")
          .attr("dominant-baseline", "center")
          .attr("text-anchor", "middle")
          .attr("font-family", fontFamily ? fontFamily["bold"] : "Helvetica")
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
        .attr("pointer-events", "none")
        .data(labels)
        .enter()
        .append("g")
        .attr("cursor", disable ? null : "pointer");

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
        .attr("font-family", fontFamily ? fontFamily["regular"] : "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr("x", LEGEND_SQUARE_LENGTH + 10)
        .attr("y", (d, i) => i * STEP + offset + textOffset)
        .text((d) => d["label"]);

      const findLabelValue = (mouseY) => {
        const index = Math.round((mouseY - offset) / STEP);

        if (0 <= index && index < labels.length) {
          return labels[index]["value"];
        }

        return null;
      };

      const mousemove = (d, i, e) => {
        const mouseY = d3.mouse(e[0])[1];

        const label = findLabelValue(mouseY);

        onHover(label);
      };

      const mouseout = () => {
        onHover(null);
      };

      const click = (d, i, e) => {
        const mouseY = d3.mouse(e[0])[1];

        const label = findLabelValue(mouseY);

        const selectedValue = label === selected ? null : label;
        setSelected(selectedValue);
        onClick(selectedValue);
      };

      svg
        .on("mousemove", (d, i, e) => {
          if (disable) {
            return;
          }
          mousemove(d, i, e);
        })
        .on("mouseout", () => {
          if (disable) {
            return;
          }
          mouseout();
        })
        .on("click", (d, i, e) => {
          if (disable) {
            return;
          }
          click(d, i, e);
        });
    },
    legendWidth,
    legendHeight,
    [selected, disable]
  );

  return <svg ref={svgRef} />;
};

export default VerticalLegend;
