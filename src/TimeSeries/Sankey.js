import React, { useState } from "react";
import * as d3 from "d3";

import { curveBumpX } from "d3-shape";

import _ from "lodash";
import { useD3 } from "../utils/useD3";

import sortAlphanumeric from "../utils/sortAlphanumeric";

const PADDING = 15;

const NODE_WIDTH = 16;
const NODE_VERTICAL_PADDING = 10;
const NULL_COLOR = "#d2d7d3";

const AXIS_HEIGHT = 40;

const Sankey = ({
  data,
  subsetParam,
  cloneParam,
  width,
  height,
  colorScale = null,
  timepointParam = "timepoint",
  timepointOrder = null,
}) => {
  const chartWidth = width - 2 * PADDING;
  const chartHeight = height - AXIS_HEIGHT;

  const [hoveredNode, setHoveredNode] = useState(null);
  const [selectedNode, setSelectedNode] = useState(null);
  const [hoveredLink, setHoveredLink] = useState(null);
  const [selectedLink, setSelectedLink] = useState(null);
  const [hoveredSourceNode, setHoveredSourceNode] = useState(null);
  const [selectedSourceNode, setSelectedSourceNode] = useState(null);
  const [hoveredTargetNode, setHoveredTargetNode] = useState(null);
  const [selectedTargetNode, setSelectedTargetNode] = useState(null);

  // Data processing - this counts the number of cells in each subtype and clone

  const timepoints =
    timepointOrder ||
    _.uniq(data.map((datum) => datum[timepointParam])).sort(sortAlphanumeric);

  const subsets = _.uniq(data.map((datum) => datum[subsetParam])).sort(
    sortAlphanumeric
  );

  const groupedTime = _.groupBy(data, (datum) => datum[timepointParam]);

  const counts = timepoints.reduce((rsf, timepoint) => {
    const timeData = groupedTime[timepoint];

    const subsetValues = _.uniq(
      timeData.map((datum) => datum[subsetParam])
    ).sort(sortAlphanumeric);

    const subsetCounts = processData(timeData, subsetParam, cloneParam);

    const cloneCounts = processData(timeData, cloneParam, subsetParam);

    return {
      ...rsf,
      [timepoint]: {
        total: timeData.length,
        clones: cloneCounts,
        subsets: subsetCounts,
        subsetValues: subsetValues,
      },
    };
  }, {});

  // Visualization processing

  const timeScale = d3
    .scalePoint()
    .domain(timepoints)
    .range([PADDING + NODE_WIDTH / 2, PADDING + chartWidth - NODE_WIDTH / 2]);

  const subsetScales = timepoints.reduce((scaleMap, timepoint) => {
    const timeData = counts[timepoint];
    const subsetsInTime = timeData["subsetValues"];

    const heightScale = d3
      .scaleLinear()
      .domain([0, timeData["total"]])
      .range([
        0,
        chartHeight - (subsetsInTime.length - 1) * NODE_VERTICAL_PADDING,
      ]);

    let yPos = PADDING;
    let yScaleObj = {};

    subsetsInTime.forEach((subset) => {
      yScaleObj[subset] = yPos;

      yPos += heightScale(timeData["subsets"][subset]["total"]);
      yPos += NODE_VERTICAL_PADDING;
    });

    return {
      ...scaleMap,
      [timepoint]: { height: heightScale, y: (subset) => yScaleObj[subset] },
    };
  }, {});

  // Processing data to draw
  const color = colorScale || d3.scaleOrdinal(subsets, d3.schemeTableau10);

  const nodes = timepoints.reduce((rsf, timepoint) => {
    const timeData = counts[timepoint];

    const subsetsInTime = timeData["subsetValues"];

    const scales = subsetScales[timepoint];
    const records = subsetsInTime.map((subset) => ({
      id: `${timepoint}_${subset}`,
      x0: timeScale(timepoint),
      width: NODE_WIDTH,
      y0: scales.y(subset),
      height: scales.height(timeData["subsets"][subset]["total"]),
      color: color(subset),
      label: subset,
    }));

    return [...rsf, ...records];
  }, []);

  // More data processing, but this time to calculate the cell ratios for drawing, which needs information from the timepoints before and after
  const linkRatios = timepoints.reduce((rsf, timepoint, index) => {
    const sourceData = counts[timepoint];
    const subsetValues = counts[timepoint]["subsetValues"];

    const ratios = subsetValues.reduce((rsf, subset) => {
      let record = {};

      if (index > 0) {
        record["before"] = processPairSubsetData(
          sourceData["subsets"][subset],
          subset,
          counts[timepoints[index - 1]]
        );
      }

      if (index < timepoints.length - 1) {
        record["after"] = processPairSubsetData(
          sourceData["subsets"][subset],
          subset,
          counts[timepoints[index + 1]]
        );
      }

      return {
        ...rsf,
        [subset]: record,
      };
    }, {});

    return { ...rsf, [timepoint]: ratios };
  }, {});

  const timepointsButLast =
    timepoints.length > 2
      ? timepoints[(0, timepoints.length - 2)]
      : [timepoints[0]];

  // for each timepoint
  // for each subtype
  // generate links for each target (subtype in next timepoint)
  const links = timepointsButLast.reduce(
    (rsf, sourceTimepoint, index) => {
      const sourceData = counts[sourceTimepoint];
      const sourceLinkData = linkRatios[sourceTimepoint];
      const sourceScales = subsetScales[sourceTimepoint];

      const targetTimepoint = timepoints[index + 1];
      const targetData = counts[targetTimepoint];
      const targetScales = subsetScales[targetTimepoint];
      const targetLinkData = linkRatios[targetTimepoint];

      const x0 = timeScale(sourceTimepoint) + NODE_WIDTH;
      const x1 = timeScale(targetTimepoint);

      const sourceSubsets = sourceData["subsetValues"];

      const rsfRecord = sourceSubsets.reduce((rsf2, sourceSubset) => {
        const subsetLinkData = sourceLinkData[sourceSubset]["after"];

        const rsf2Record = targetData["subsetValues"].map((targetSubset) => {
          return {
            id: `${sourceTimepoint}_${sourceSubset}_${targetSubset}`,
            source: `${sourceTimepoint}_${sourceSubset}`,
            target: `${targetTimepoint}_${targetSubset}`,
            x0,
            x1,
            y0: calculateLinkYStart(
              sourceScales,
              subsetLinkData,
              sourceSubset,
              targetSubset
            ),
            h0: calculateLinkHeight(sourceScales, subsetLinkData, targetSubset),
            y1: calculateLinkYStart(
              targetScales,
              targetLinkData[targetSubset]["before"],
              targetSubset,
              sourceSubset
            ),
            h1: calculateLinkHeight(
              targetScales,
              targetLinkData[targetSubset]["before"],
              sourceSubset
            ),
            c0: color(sourceSubset),
            c1: color(targetSubset),
          };
        });

        return [...rsf2, ...rsf2Record];
      }, []);

      return [...rsf, ...rsfRecord];
    },

    []
  );

  const ref = useD3(
    (svg) => {
      const node = svg
        .append("g")
        .selectAll("rect")
        .data(nodes)
        .join("rect")
        .attr("x", (d) => d.x0)
        .attr("y", (d) => d.y0)
        .attr("height", (d) => d.height)
        .attr("width", (d) => d.width)
        .attr("fill", (d) =>
          isIDHighlighted(d.id, [
            hoveredNode,
            selectedNode,
            hoveredSourceNode,
            hoveredTargetNode,
            selectedSourceNode,
            selectedTargetNode,
          ])
            ? d.color
            : NULL_COLOR
        )
        .attr("stroke", (d) =>
          isIDHighlighted(d.id, [
            hoveredNode,
            selectedNode,
            hoveredSourceNode,
            hoveredTargetNode,
            selectedSourceNode,
            selectedTargetNode,
          ])
            ? selectedNode === d.id
              ? "currentColor"
              : d.color
            : NULL_COLOR
        )
        .on("mouseenter", (d, i, e) => {
          // if something is selected, do nothing
          if (selectedNode || selectedLink) {
            return;
          }
          setHoveredNode(d.id);
          setHoveredLink(null);
        })
        .on("mouseout", (d, i, e) => {
          // if something is selected, do nothing
          if (selectedNode || selectedLink) {
            return;
          }
          setHoveredNode(null);
        })
        .on("click", (d, i, e) => {
          if (selectedNode === d.id) {
            setSelectedNode(null);
          } else {
            setSelectedNode(d.id);
          }
          setHoveredNode(null);
          setHoveredLink(null);
          setSelectedLink(null);
          setSelectedSourceNode(null);
          setSelectedTargetNode(null);
        });

      const link = svg
        .append("g")
        .selectAll("g")
        .data(links)
        .join("g")
        .style("mix-blend-mode", "multiply");

      const uid = `O-${Math.random().toString(16).slice(2)}`;
      link
        .append("linearGradient")
        .attr("id", (d) => `${uid}-link-${d.id}`)
        .attr("gradientUnits", "userSpaceOnUse")
        .attr("x1", (d) => d.x0)
        .attr("x2", (d) => d.x1)
        .call((gradient) =>
          gradient
            .append("stop")
            .attr("offset", "0%")
            .attr("stop-color", (d) => d.c0)
        )
        .call((gradient) =>
          gradient
            .append("stop")
            .attr("offset", "100%")
            .attr("stop-color", (d) => d.c1)
        );

      const area = d3
        .area()
        .curve(curveBumpX)
        .x((d) => d[0])
        .y0((d) => d[1])
        .y1((d) => d[2]);

      link
        .append("path")
        .attr("d", (d) =>
          area([
            [d.x0, d.y0, d.y0 + d.h0],
            [d.x1, d.y1, d.y1 + d.h1],
          ])
        )
        .attr("fill-opacity", 0.8)
        .attr("fill", (d) =>
          [hoveredLink, selectedLink, hoveredNode, selectedNode].every(
            (value) => !value
          ) ||
          [hoveredLink, selectedLink].some((value) => d.id === value) ||
          [hoveredNode, selectedNode].some((value) => d.source === value) ||
          [hoveredNode, selectedNode].some((value) => d.target === value)
            ? `url(#${uid}-link-${d.id})`
            : NULL_COLOR
        )
        .attr("stroke-width", (d) => (selectedLink === d.id ? 2 : 0))
        .attr("stroke", (d) => `url(#${uid}-link-${d.id})`)
        .on("mouseenter", (d, i, e) => {
          // if something is selected, do nothing
          if (selectedNode || selectedLink) {
            return;
          }
          setHoveredNode(null);
          setHoveredLink(d.id);
          setHoveredSourceNode(d.source);
          setHoveredTargetNode(d.target);
        })
        .on("mouseout", (d, i, e) => {
          // if something is selected, do nothing
          if (selectedNode || selectedLink) {
            return;
          }
          setHoveredLink(null);
          setHoveredSourceNode(null);
          setHoveredTargetNode(null);
        })
        .on("click", (d, i, e) => {
          if (selectedLink === d.id) {
            setSelectedLink(null);
            setSelectedSourceNode(null);
            setSelectedTargetNode(null);
          } else {
            setSelectedLink(d.id);
            setSelectedSourceNode(d.source);
            setSelectedTargetNode(d.target);
          }
          setHoveredNode(null);
          setHoveredLink(null);
          setSelectedNode(null);
        });

      svg
        .append("g")
        .attr("font-family", "sans-serif")
        .attr("font-size", 11)
        .selectAll("text")
        .data(nodes)
        .join("text")
        .attr("x", (d) => (d.x0 < width / 2 ? d.x0 + d.width + 6 : d.x0 - 6))
        .attr("y", (d) => d.y0 + d.height / 2)
        .attr("dy", "0.35em")
        .attr("text-anchor", (d) => (d.x0 < width / 2 ? "start" : "end"))
        .attr("pointer-events", "none")
        .text((d) => d.label);

      const xAxis = d3.axisBottom(timeScale).tickSize(0);

      svg
        .append("g")
        .style(
          "transform",
          `translate(${NODE_WIDTH / 2}px, ${chartHeight + PADDING + 2}px)`
        )
        .call(xAxis)
        .call((g) => g.select(".domain").remove())
        .call((g) =>
          g
            .selectAll(".tick text")
            .attr("font-family", "sans-serif")
            .attr("font-size", 12)
            .attr("font-weight", 500)
        );
    },
    width,
    height,
    [links, width, height]
  );

  return <svg ref={ref} />;
};

// For initial counting - this iterates through the data and keeps track of (say groupParam = subtype, countParam = clone)
// all subtype values
// for each subtype, all clone values
// for each subtype, total number of cells
// for each subtype, total number of cells per clone
// for each subtype, summation of cells per clone (helpful for y placement)
const processData = (data, groupParam, countParam) => {
  const groupValues = _.uniq(data.map((datum) => datum[groupParam])).sort(
    sortAlphanumeric
  );

  const groupedData = _.groupBy(data, (datum) => datum[groupParam]);

  const counts = groupValues.reduce((rsf, value) => {
    const groupData = groupedData[value];

    const groupCounts = _.countBy(groupData, (datum) => datum[countParam]);

    const total = groupData.length;

    const countValues = Object.keys(groupCounts).sort(sortAlphanumeric);

    let ssf = 0;
    let cumulativeCount = {};

    countValues.forEach((countValue) => {
      cumulativeCount[countValue] = ssf;

      ssf += groupCounts[countValue];
    });

    return {
      ...rsf,
      [value]: {
        counts: groupCounts,
        sum: { ...cumulativeCount },
        values: countValues,
        total,
      },
    };
  }, {});

  return counts;
};

// For link ratios
// The timepoints only share clone data, so when calculating how much space to allocate to each link, you need to figure out how many cells to send from timepoint 1 to timepoint 2
// For example:
// timepoint 1, source subtype, clone 1 = 100 cells
// timepoint 2, target subtype 1, clone 1 = 20 cells
// timepoint 2, target subtype 2, clone 1 = 30 cells

// then from source subtype, you will want to send: (20/50) -> 40% of the cells (=40 cells) to target subtype 1, and (30/50) => 60% of the cells (=60 cells) to target subtype 2

// Then to figure out how many cells to send from source subtype to target subtype, you need to sum up this ratio across all clones shared by source and target subtypes

// That total/count scaled proportionally provides the height for the link
// The y placement would be the summation of heights of all target subtypes before it - this is provided by 'sums' in the final object

const processPairSubsetData = (sourceSubsetData, sourceSubset, targetData) => {
  const targetSubsets = targetData["subsetValues"];
  const sourceClones = sourceSubsetData["values"];

  let targetSum = 0;
  let sumObj = {};

  const record = targetSubsets.reduce((rsf2, targetSubset) => {
    const targetSubsetData = targetData["subsets"][targetSubset];
    const clones = targetSubsetData["values"].filter((clone) =>
      sourceClones.includes(clone)
    );

    let cellSum = 0;

    const rsf2Record = clones.reduce((rsf3, clone) => {
      const cellCount =
        (sourceSubsetData["counts"][clone] *
          targetData["clones"][clone]["counts"][targetSubset]) /
        targetData["clones"][clone]["total"];

      const rsf3Record = {
        count: cellCount,
        sum: cellSum,
      };

      cellSum += cellCount;

      return {
        ...rsf3,
        [clone]: rsf3Record,
      };
    }, {});

    sumObj[targetSubset] = targetSum;
    targetSum += cellSum;

    return {
      ...rsf2,
      [targetSubset]: {
        total: cellSum,
        clones: rsf2Record,
        cloneValues: clones,
      },
    };
  }, {});

  return { ...record, sums: sumObj };
};

const calculateLinkYStart = (
  scales,
  sourceData,
  sourceSubset,
  targetSubset
) => {
  // height of subset
  const subsetHeight = scales.y(sourceSubset);

  // offset of targetSubset start within subset
  const cloneStartHeight = scales.height(sourceData["sums"][targetSubset]);

  return subsetHeight + cloneStartHeight;
};

const calculateLinkHeight = (scales, sourceData, targetSubset) => {
  return scales.height(sourceData[targetSubset]["total"]);
};

const isIDHighlighted = (id, values) => {
  // if nothing is highlighted, return true
  if (values.every((value) => !value)) {
    return true;
  }

  return values.some((value) => value === id);
};

export default Sankey;
