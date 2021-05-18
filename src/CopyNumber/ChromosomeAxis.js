import React from "react";
import { useD3 } from "../utils/useD3";
import { BACKGROUND_COLORS } from "./utils";

export const HEIGHT = 14;

const ChromosomeAxis = ({ chromosomes, width }) => {
  const bpTotal = chromosomes.reduce((currSum, chr) => currSum + chr.length, 0);
  const ref = useD3(
    (svg) => {
      const bpRatio = width / bpTotal;

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
        .attr("height", HEIGHT)
        .attr(
          "fill",
          (d, i) => BACKGROUND_COLORS[i % BACKGROUND_COLORS.length]
        );

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
        .attr("y", HEIGHT / 2)
        .attr("alignment-baseline", "middle")
        .attr("text-anchor", "middle")
        .attr("font-size", "10px")
        .text((d) => d.chr);
    },
    width,
    HEIGHT,
    []
  );
  return <svg ref={ref} width={width} height={HEIGHT} />;
};

export default ChromosomeAxis;
