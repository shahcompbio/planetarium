import React from "react";
import { useD3 } from "../../utils/useD3";
import {
  getGenomeYScale,
  getCopyNumberArray,
  BACKGROUND_COLORS,
} from "../utils";

const ProfileBackground = ({ chromosomes, width, height, maxState }) => {
  const bpTotal = chromosomes.reduce((currSum, chr) => currSum + chr.length, 0);

  const ref = useD3(
    (svg) => {
      const bpRatio = width / bpTotal;
      const copyNumbers = getCopyNumberArray(maxState);
      const genomeYScale = getGenomeYScale(maxState, height);

      svg
        .append("g")
        .attr("class", "gw-background")
        .selectAll(".gw-background-box")
        .data(chromosomes)
        .enter()
        .append("rect")
        .attr("class", (d) => "gw-background-box chrom-" + d.chr)
        .attr("x", (d) => d.genomeStart * bpRatio)
        .attr("y", 0)
        .attr("width", (d) => d.length * bpRatio)
        .attr("height", height)
        .attr(
          "fill",
          (d, i) => BACKGROUND_COLORS[i % BACKGROUND_COLORS.length]
        );

      svg
        .append("g")
        .attr("class", "gw-background-lines")
        .selectAll(".gw-background-line")
        .data(copyNumbers)
        .enter()
        .append("line")
        .attr("x1", 0)
        .attr("x2", width)
        .attr("y1", (d) => genomeYScale(d))
        .attr("y2", (d) => genomeYScale(d))
        .attr("stroke", "#aaa")
        .attr("stroke-dasharray", "1, 2");
    },
    width,
    height,
    []
  );

  return (
    <svg
      id="profile-background"
      ref={ref}
      width={width}
      height={height}
      style={{ position: "absolute" }}
    />
  );
};

export default ProfileBackground;
