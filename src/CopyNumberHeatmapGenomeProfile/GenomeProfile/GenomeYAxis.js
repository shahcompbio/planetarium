import React from "react";
import {useD3} from "../../utils/useD3";
import * as d3 from "d3";
import { getGenomeYScale, TOP_PADDING } from "./utils";

const GenomeYAxis = ({ width, height, maxState }) => {
  const ref = useD3((svg) => {
    const genomeYScale = getGenomeYScale(maxState, height);

    const yAxisTicks = genomeYScale
      .ticks()
      .filter((tick) => Number.isInteger(tick));

    var yAxis = d3
      .axisLeft(genomeYScale)
      .tickValues(yAxisTicks)
      .tickFormat(d3.format("d"));

    svg
      .append("g")
      .attr("class", "genome-y-axis")
      .style("transform", `translate(${width - 2}px, ${TOP_PADDING}px)`)
      .call(yAxis);
  },width,height,[]);
  return <svg id="genome-axis" ref={ref} width={width} height={height} />;
};

export default GenomeYAxis;
