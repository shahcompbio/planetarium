import React, { useState } from "react";
import PropTypes from "prop-types";
import _ from "lodash";

import * as d3 from "d3";
import { useD3 } from "../utils/useD3";
import { Grid } from "@mui/material";
import isHighlighted from "../utils/isHighlighted";
import sortAlphanumeric from "../utils/sortAlphanumeric";
import VerticalLegend from "../Legend/Vertical";

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
const addFakeTimepoints = (data, timepointParam) =>
  data.reduce((final, d) => {
    var dup = Object.assign({}, d);
    dup[timepointParam] = dup[timepointParam] + ".5";

    final = [...final, d, dup];
    return final;
  }, []);

const addGradientToSvg = (svg, preColor, postColor, name) => {
  var lg = svg
    .append("defs")
    .append("linearGradient")
    .attr("id", "mygrad-" + name)
    .attr("x1", "0%")
    .attr("x2", "100%")
    .attr("y1", "0%")
    .attr("y2", "0%"); //since its a vertical linear gradient
  lg.append("stop")
    .attr("offset", "0%")
    .style("stop-color", preColor) //end in red
    .style("stop-opacity", 1);

  lg.append("stop")
    .attr("offset", "100%")
    .style("stop-color", postColor) //start in blue
    .style("stop-opacity", 1);
};
const Fishtail = ({
  data,
  subsetParam,
  timepointParam = "timepoint",
  width = 700,
  height = 400,
  timepoint = null,
  timepointOrder = null,
  legendSorting = null,
  subset = null,
  disable = false,
  colorScale = null,
  cloneColor = null,
  addTwoTimepointCurve = false,
  interpolateColor = false,
  normalize = true,
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

  const orgData =
    addTwoTimepointCurve && timepointOrder.length === 2
      ? addFakeTimepoints(data, timepointParam)
      : data;

  const timeValues = timepointOrder
    ? addTwoTimepointCurve
      ? timepointOrder.reduce(
          (final, t) => (final = [...final, t, t + ".5"]),
          []
        )
      : timepointOrder
    : _.uniq(orgData.map((datum) => datum[timepointParam])).sort(
        sortAlphanumeric
      );

  const subsetValues = legendSorting
    ? _.uniq(orgData.map((datum) => datum[subsetParam])).sort(
        (a, b) => legendSorting.indexOf(a) - legendSorting.indexOf(b)
      )
    : _.uniq(orgData.map((datum) => datum[subsetParam])).sort(sortAlphanumeric);

  const subsetOverall = highlightedSubset || subset || selectedSubset;
  const timepointOverall =
    highlightedTimepoint || timepoint || selectedTimepoint;

  const counts = timeValues.map((timepoint) => {
    const timeData = orgData.filter(
      (datum) => datum[timepointParam] === timepoint
    );
    //2 - 0.6667
    //3 - 0.3846
    //4 - 0.709

    const percent = timepoint.indexOf("F") !== -1 ? 0.6667 : 1;
    console.log(timepoint);
    console.log(percent);
    const subsetCounts = normalize
      ? subsetValues.map(
          (subset) =>
            timeData.filter((datum) => datum[subsetParam] === subset).length /
            timeData.length
        )
      : subsetValues.map(
          (subset) =>
            Math.round(
              timeData.filter((datum) => datum[subsetParam] === subset).length *
                percent
            ) / timeData.length
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
    .domain([...timeValues])
    .range([PADDING, PADDING + chartWidth]);

  const color =
    colorScale ||
    d3
      .scaleOrdinal()
      .domain(subsetValues)
      .range(
        COLOR_ARRAY.slice(0, Math.min(subsetValues.length, COLOR_ARRAY.length))
      );
  console.log(subsetValues);
  /*const order = [
    "C26",
    "C7",
    "C12",
    "C4",
    "C18",
    "C24",
    "C14",
    "C9",
    "C6",
    "C0",
    "C13",
    "C16",
    "C23",
    "C10",
    "C2",
    "C19",
    "C5",
    "C21",
    "C8",
    "C11",
    "C15",
    "C1",
    "C22",
    "C17",
    "C20",
    "C27",
    "C28",
    "C25",
    "C3",
  ].map((d) => "04_04_" + d);*/
  /*const order = [
    "C1",
    "C3",
    "C6",
    "C8",
    "C0",
    "C7",
    "C11",
    "C12",
    "C15",
    "C14",
    "C5",
    "C4",
    "C13",
    "C2",
    "C10",
    "C9",
  ].map((d) => "03_03_" + d);*/
  /*    "C6",
      "C17",
      "C16",
      "C5",
      "C1",
      "C8",
      "C4",
      "C9",
      "C13",
      "C3",
      "C19",
      "C12",
      "C11",
      "C20",
      "C0",
      "C2",
      "C10",
      "C7",
      "C14",
      "C15",
      "C18",*/

  /*const order = [
    "C6",
    "C5",
    "C17",
    "C16",
    "C4",
    "C1",
    "C9",
    "C8",
    "C13",
    "C3",
    "C19",
    "C12",
    "C11",
    "C20",
    "C0",
    "C10",
    "C7",
    "C2",
    "C14",
    "C15",
    "C13",
    "C2",
    "C18",
  ].map((d) => "01_" + d);*/
  const order = legendSorting.map((d) => d.split("-")[0]);
  console.log(counts);
  console.log(order);
  const drawArea = (svg) => {
    const series = d3
      .stack()
      .keys(order)
      .offset(d3.stackOffsetNone)
      .order(d3.stackOrderNone)(counts);
    //d3.stackOrderAscending

    const yScale = d3
      .scaleLinear()
      .domain([
        d3.min(series, (d) => d3.min(d, (d) => d[0])),
        d3.max(series, (d) => d3.max(d, (d) => d[1])),
      ])
      .range([0, height - AXIS_HEIGHT]);

    const area = d3
      .area()
      .x((d) => {
        return timeScale(d.data.timepoint);
      })
      .y0((d) => yScale(d[0]))
      .y1((d) => (Number.isNaN(d[1]) ? yScale(d[0]) : yScale(d[1])))
      .curve(d3.curveCardinal);

    svg
      .append("g")
      .attr("pointer-events", "none")
      .selectAll("path")
      .data(series)
      .join("path")
      .attr("fill", ({ key }) => {
        if (interpolateColor) {
          const pre = timepointOrder[0];
          const post = timepointOrder[1];
          const preNum = cloneColor[key + "-" + pre];
          const postNum = cloneColor[key + "-" + post];
          const preColor = color(preNum);
          const postColor = color(postNum);

          addGradientToSvg(svg, preColor, postColor, key);
          return "url(#mygrad-" + key + ")";
        } else {
          return isHighlighted(key, subsetOverall)
            ? cloneColor !== null
              ? color(cloneColor[key])
              : color(key)
            : NULL_AREA_COLOR;
        }
      })
      .attr("id", ({ key }) => key)
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
      orgData,
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
          ticks={subsetValues}
          colorScale={color}
          cloneColor={cloneColor}
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

Fishtail.propTypes = {
  /**
   * Data points to visualize
   */
  data: PropTypes.arrayOf(PropTypes.object).isRequired,

  /**
   * Key in data that specifies subgroup
   */
  subsetParam: PropTypes.string.isRequired,

  /**
   * Key in data that specifies timepoint
   */
  timepointParam: PropTypes.string,

  /**
   * Width of plot
   */
  width: PropTypes.number,

  /**
   * Height of plot
   */
  height: PropTypes.number,

  /**
   * Value of timepoint to highlight
   */
  timepoint: PropTypes.string,

  /**
   * Array consisting of order of timepoints
   **/
  timepointOrder: PropTypes.array,

  /**
   * Order of legend options to sort by
   **/
  legendSorting: PropTypes.array,

  /**
   * Value of subset to highlight
   */
  subset: PropTypes.string,

  /**
   * To disable interactions on plot
   */
  disable: PropTypes.bool,
  /**
   * Override default color scale of subsets
   */
  colorScale: PropTypes.func,

  /**
   * Color of gradient between two timepoints
   **/
  //cloneColor: PropTypes.arrayOf(PropTypes.object),

  /**
   * Option to add curve to plot by adding fake timepoints
   **/
  addTwoTimepointCurve: PropTypes.bool,

  /**
   * Option to interpolate color, must include cloneColor
   **/
  interpolateColor: PropTypes.bool,

  /**
   * Option to remove normalization of strands
   **/
  normalize: PropTypes.bool,

  /**
   * Callback when value on legend is hovered
   */
  onLegendHover: PropTypes.func,

  /**
   * Callback when value on legend is clicked
   */
  onLegendClick: PropTypes.func,

  /**
   * Callback when timepoint section is hovered
   */
  onTimepointHover: PropTypes.func,

  /**
   * Callback when timepoint section is clicked
   */
  onTimepointClick: PropTypes.func,
};

export default Fishtail;
