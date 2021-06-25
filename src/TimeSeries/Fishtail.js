import React, { useState } from "react";
import _ from "lodash";

import * as d3 from "d3";
import { useD3 } from "../utils/useD3";
import { Grid } from "@material-ui/core";
import { isValueHighlighted as isHighlighted } from "../utils/isHighlighted";
import sortAlphanumeric from "../utils/sortAlphanumeric";
import VerticalLegend from "../Legend/VerticalLegend";

const AXIS_HEIGHT = 20;
const PADDING = 10;

const NULL_AREA_COLOR = "#d2d7d3";

const SLIDER_BAR_WIDTH = 6;
const SLIDER_BAR_COLOR = "#ebebeb";
const HIGHLIGHTED_BAR_COLOR = "#5e5e5e";

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

const Fishtail = ({
  data,
  subsetParam,
  timepointParam = "timepoint",
  width,
  height,
  timepoint = null,
  subset = null,
  disable = false,
  colorScale = null,
  onLegendHover = (value) => {},
  onLegendClick = (value) => {},
  onTimepointHover = (timepoint) => {},
  onTimepointClick = (timepoint) => {},
}) => {
  const [highlightedTimepoint, setHighlightedTimepoint] = useState(null);
  const [selectedTimepoint, setSelectedTimepoint] = useState(null);
  const [highlightedSubset, setHighlightedSubset] = useState(null);
  const [selectedSubset, setSelectedSubset] = useState(null);
  const chartWidth = width - 2 * PADDING;
  const chartHeight = height - AXIS_HEIGHT;

  const timeValues = _.uniq(data.map((datum) => datum[timepointParam])).sort(
    sortAlphanumeric
  );

  const subsetValues = _.uniq(data.map((datum) => datum[subsetParam])).sort(
    sortAlphanumeric
  );

  const subsetOverall = highlightedSubset || subset || selectedSubset;
  const timepointOverall =
    highlightedTimepoint || timepoint || selectedTimepoint;

  const counts = timeValues.map((timepoint) => {
    const timeData = data.filter(
      (datum) => datum[timepointParam] === timepoint
    );

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

  const color =
    colorScale ||
    d3
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
      .attr("pointer-events", "none")
      .selectAll("path")
      .data(series)
      .join("path")
      .attr("fill", ({ key }) =>
        isHighlighted(key, subsetOverall) ? color(key) : NULL_AREA_COLOR
      )
      .attr("d", area)
      .attr("stroke", "#FFFFFF")
      .attr("stroke-width", 1)
      .attr("stroke-opacity", 0.2)
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
      .attr("pointer-events", "none")
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

    const mousemove = (d, i, e) => {
      const mouseX = d3.mouse(e[0])[0];

      const timepointIndex = Math.round(mouseX / timeScale.step());
      const timepointValue = timeValues[timepointIndex];

      setHighlightedTimepoint(timepointValue);
      onTimepointHover(timepointValue);
    };

    const mousedown = (d, i, e) => {
      const mouseX = d3.mouse(e[0])[0];

      const timepointIndex = Math.round(mouseX / timeScale.step());
      const timepointValue = timeValues[timepointIndex];
      const selectedValue = timepointValue === selected ? null : timepointValue;

      setSelectedTimepoint(selectedValue);
      onTimepointClick(selectedValue);
    };
    svg
      .on("mousemove", (d, i, e) => {
        if (disable) {
          return;
        }
        mousemove(d, i, e);
      })
      .on("click", (d, i, e) => {
        if (disable) {
          return;
        }
        mousedown(d, i, e);
      })
      .on("mouseout", () => {
        if (disable) {
          return;
        }
        setHighlightedTimepoint(null);
        onTimepointHover(null);
      });
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg);
      drawAxis(svg);
      drawSlider(svg, highlightedTimepoint, timepoint || selectedTimepoint);
    },
    width,
    height,
    [
      data,
      highlightedTimepoint,
      selectedTimepoint,
      timepoint,
      subsetOverall,
      disable,
    ]
  );

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <svg ref={ref} />
      </Grid>
      <Grid item>
        <VerticalLegend
          width={100}
          title={subsetParam}
          labels={subsetValues.map((value) => ({
            value,
            label: value,
            color: color(value),
          }))}
          disable={disable}
          onHover={(value) => {
            setHighlightedSubset(value);
            onLegendHover(value);
          }}
          onClick={(value) => {
            setSelectedSubset(value);
            onLegendClick(value);
          }}
        />
      </Grid>
    </Grid>
  );
};

export default Fishtail;
