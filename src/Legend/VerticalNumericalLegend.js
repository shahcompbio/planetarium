import React from "react";
import * as d3 from "d3";

import { useD3 } from "../utils/useD3";

const PADDING = 20;
const STOPS = 10;
const format = d3.format(".2f");

const VerticalNumericalLegend = ({
  width,
  height,
  title,
  colorScale,
  ticks = 2,
}) => {
  const ID = Math.random().toString();
  const gradWidth = width - 50;
  const svgRef = useD3(
    (svg) => {
      var gradHeight = height;
      var offset = 0;
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
          .attr("x", gradWidth / 2)
          .attr("y", PADDING)
          .text(title);

        const titleHeight = PADDING + titleText.node().getBBox().height;
        gradHeight -= titleHeight;
        offset = titleHeight;
      }
      const numTicks = Math.max(ticks, 2);
      const tickArray = Array.from(Array(numTicks).keys()).map((i) => i + 1);
      const tickToDomain = d3
        .scaleLinear()
        .range(colorScale.domain())
        .domain([1, numTicks]);
      const tickToY = d3
        .scaleLinear()
        .range([offset, offset + gradHeight - 10])
        .domain([numTicks, 1]);

      const stopArray = Array.from(Array(STOPS).keys()).map((i) => i + 1);
      const stopToDomain = d3
        .scaleLinear()
        .range(colorScale.domain())
        .domain([1, STOPS]);
      const stopToPer = d3.scaleLinear().range([0, 100]).domain([1, STOPS]);

      var grad = svg
        .append("defs")
        .append("linearGradient")
        .attr("id", ID)
        .attr("x1", "0%")
        .attr("x2", "0%")
        .attr("y1", "100%")
        .attr("y2", "0%");

      grad
        .selectAll("stop")
        .data(stopArray)
        .enter()
        .append("stop")
        .style("stop-color", (d) => colorScale(stopToDomain(d)))
        .attr("offset", (d) => `${stopToPer(d)}%`);

      svg
        .append("rect")
        .attr("x", 0)
        .attr("y", offset)
        .attr("width", gradWidth)
        .attr("height", gradHeight)
        .style("fill", `url(#${ID})`);

      svg
        .append("g")
        .selectAll("text")
        .data(tickArray)
        .enter()
        .append("text")
        .attr("alignment-baseline", "hanging")
        .attr("dominant-baseline", "hanging")
        .attr("text-align", "left")
        .attr("font-family", "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr("x", gradWidth + 10)
        .attr("y", (d) => tickToY(d))
        .text((d) => format(tickToDomain(d)));
    },
    width,
    height,
    [title, colorScale, ticks]
  );
  return <svg ref={svgRef} />;
};

export default VerticalNumericalLegend;
