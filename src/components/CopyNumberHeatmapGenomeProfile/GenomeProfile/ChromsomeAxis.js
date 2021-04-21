import React from "react";
import { Y_AXIS_WIDTH } from "./utils";
import {useD3} from "../../utils/useD3";

const ChromosomeAxis = ({ chromosomes, bpTotal, width, height }) => {
  const ref = useD3((svg) => {
    const backgroundColors = ["#fefefe", "#eee"];
    const bpRatio = width / bpTotal;

    svg.style("transform", `translate(${Y_AXIS_WIDTH}px, 0px)`);

    svg
      .append("g")
      .attr("class", "chromosome-axis-boxes")
      .selectAll(".chromosome-axis-box")
      .data(chromosomes)
      .enter()
      .append("rect")
      .attr("class", (d) => "chromosome-axis-box chrom-" + d.chr)
      .attr("x", (d) => d.genomeStart * bpRatio)
      .attr("y", 0)
      .attr("width", (d) => d.length * bpRatio)
      .attr("height", height)
      .attr("fill", (d, i) => backgroundColors[i % 2]);

    svg
      .append("g")
      .attr("class", "chromosome-axis-texts")
      .selectAll(".chromosome-axis-text")
      .data(chromosomes)
      .enter()
      .append("text")
      .attr("class", (d) => "chromosome-axis-text chrom-" + d.chr)
      .style("fill", "black")
      .attr("x", (d) => (d.genomeStart + d.length / 2) * bpRatio)
      .attr("y", height)
      .attr("text-anchor", "middle")
      .attr("font-size", "10px")
      .text((d) => d.chr);
  });
  return <svg id="chromosome-axis" ref={ref} width={width} height={height} />;
};

export default ChromosomeAxis;