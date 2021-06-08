import React, { useState } from "react";
import _ from "lodash";

import * as d3 from "d3";
import { useD3 } from "../utils/useD3";
import { isValueHighlighted as isHighlighted } from "../utils/isHighlighted";

const AXIS_HEIGHT = 20;
const PADDING = 10;

const SLIDER_BAR_WIDTH = 6;
const SLIDER_BAR_COLOR = "#ebebeb";
const HIGHLIGHTED_BAR_COLOR = "#5e5e5e";
const SLIDE_BAR_STROKE_COLOR = "#ffcb3d";

const COLOR_ARRAY = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];

const Fishtail = ({ data, subsetParam, width, height }) => {
  const [highlightedTimepoint, setHighlightedTimepoint] = useState(null);
  const [selectedTimepoint, setSelectedTimepoint] = useState(null);
  const chartWidth = width - 2 * PADDING;
  const chartHeight = height - AXIS_HEIGHT;

  const timeValues = _.uniq(data.map((datum) => datum["timepoint"])).sort(
    (a, b) => a.substring(1) - b.substring(1)
  );

  const subsetValues = _.uniq(data.map((datum) => datum[subsetParam])).sort();

  const counts = timeValues.map((timepoint) => {
    const timeData = data.filter((datum) => datum["timepoint"] === timepoint);

    const subsetCounts = subsetValues.map(
      (subset) =>
        timeData.filter((datum) => datum[subsetParam] === subset).length
    );

    return {
      timepoint,
      ...subsetValues.reduce(
        (rsf, subset, index) => ({
          ...rsf,
          [subset]: subsetCounts[index],
        }),
        {}
      ),
    };
  });

  const timeScale = d3
    .scalePoint()
    .domain(timeValues)
    .range([PADDING, PADDING + chartWidth]);

  const color = d3
    .scaleOrdinal()
    .domain(subsetValues)
    .range(
      COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
    );

  const drawArea = (svg) => {
    const series = d3
      .stack()
      .keys(subsetValues)
      .offset(d3.stackOffsetWiggle)
      .order(d3.stackOrderInsideOut)(counts);

    const yScale = d3
      .scaleLinear()
      .domain([
        d3.min(series, (d) => d3.min(d, (d) => d[0])),
        d3.max(series, (d) => d3.max(d, (d) => d[1])),
      ])
      .range([0, height - AXIS_HEIGHT]);

    const area = d3
      .area()
      .curve(d3.curveBasis)
      .x((d) => timeScale(d.data.timepoint))
      .y0((d) => yScale(d[0]))
      .y1((d) => (Number.isNaN(d[1]) ? yScale(d[0]) : yScale(d[1])));

    svg
      .append("g")
      .selectAll("path")
      .data(series)
      .join("path")
      .attr("fill", ({ key }) => color(key))
      .attr("d", area)
      .append("title")
      .text(({ key }) => key);
  };

  const drawAxis = (svg) => {
    const xAxis = d3.axisBottom(timeScale);

    svg
      .append("g")
      .style("transform", `translate(0px, ${chartHeight}px)`)
      .call(xAxis);
  };

  const drawSlider = (svg, highlighted, selected) => {
    svg
      .append("g")
      .selectAll("rect")
      .data(timeValues)
      .enter()
      .append("rect")
      .attr("x", (d) => timeScale(d) - SLIDER_BAR_WIDTH / 2)
      .attr("y", 0)
      .attr("width", SLIDER_BAR_WIDTH)
      .attr("height", chartHeight)
      .attr("fill", (d) =>
        highlighted !== null && highlighted === d
          ? HIGHLIGHTED_BAR_COLOR
          : SLIDER_BAR_COLOR
      )
      .attr("fill-opacity", 0.5)
      .attr("stroke", (d) =>
        selected !== null && selected === d ? HIGHLIGHTED_BAR_COLOR : null
      );

    const mousemove = () => {
      const mouseX = d3.event.clientX;

      const timepointIndex = Math.round(mouseX / timeScale.step());

      setHighlightedTimepoint(timeValues[timepointIndex]);
    };

    const mousedown = () => {
      const mouseX = d3.event.clientX;

      const timepointIndex = Math.round(mouseX / timeScale.step());

      setSelectedTimepoint(timeValues[timepointIndex]);
    };
    svg.on("mousemove", mousemove).on("click", mousedown);
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg);
      drawAxis(svg);
      drawSlider(svg, highlightedTimepoint, selectedTimepoint);
    },
    width,
    height,
    [highlightedTimepoint, selectedTimepoint]
  );

  return <svg ref={ref} />;
};

export default Fishtail;
