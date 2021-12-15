import React, { useState } from "react";
import * as d3 from "d3";

import { curveBumpX } from "d3-shape";

import _ from "lodash";
import { useD3 } from "../utils/useD3";

import isHighlighted from "../utils/isHighlighted";
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

  const [highlightedNode, setHighlightedNode] = useState(null);

  // Data processing

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

  const timepointsButLast =
    timepoints.length > 2
      ? timepoints[(0, timepoints.length - 2)]
      : [timepoints[0]];

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

  // for each time point
  // for each subset
  // for each targetSubset
  // x0: sourceSubset
  // X1: targetSubset
  // y0: sourceSubsetStartY + offset of targetSubsets prior
  // h0: summation of sourceCloneCells * proportion of target clone cells in targetSubset
  // y1: targetSubsetStartY + offset of sourceSubsets prior
  // h0: summation of targetCloneCells * proportion of source clone cells in sourceSubset

  // subset startY is sum of cells up to subset (subset sum by clone)

  // height for given subset is summation of clones in given subset * proportion each clone in targetSubset

  //  offset Y for given subset is the summation of all heights before given subset (so all targetSubsets before)

  // data should have
  // - timepoint
  // - subset
  // - count by clone
  // - sum by clone
  // - clone
  // - count by subset
  // - sum by subset
  //

  // - timepointOrder

  console.log(linkRatios);
  const links = timepointsButLast.reduce(
    (rsf, sourceTimepoint, index) => {
      const sourceData = counts[sourceTimepoint];
      const sourceLinkData = linkRatios[sourceTimepoint];
      const sourceScales = subsetScales[sourceTimepoint];

      const targetTimepoint = timepoints[index + 1];
      const targetData = counts[targetTimepoint];
      const targetScales = subsetScales[targetTimepoint];
      const targetLinkData = linkRatios[targetTimepoint];
      // for each subset-clone pair, generate a list of links

      const x0 = timeScale(sourceTimepoint) + NODE_WIDTH;
      const x1 = timeScale(targetTimepoint);

      const sourceSubsets = sourceData["subsetValues"];

      const rsfRecord = sourceSubsets.reduce((rsf2, sourceSubset) => {
        const subsetLinkData = sourceLinkData[sourceSubset]["after"];

        const rsf2Record = targetData["subsetValues"].reduce(
          (rsf3, targetSubset) => {
            const clones = subsetLinkData[targetSubset]["cloneValues"];

            const rsf3Records = clones.map((clone) => ({
              id: `${sourceTimepoint}_${sourceSubset}_${targetSubset}_${clone}`,
              source: `${sourceTimepoint}_${sourceSubset}`,
              target: `${targetTimepoint}_${targetSubset}`,
              x0,
              x1,
              y0: calculateLinkYStart(
                clone,
                sourceScales,
                subsetLinkData,
                sourceSubset,
                targetSubset
              ),
              h0: calculateLinkHeight(
                clone,
                sourceScales,
                subsetLinkData,
                targetSubset
              ),
              y1: calculateLinkYStart(
                clone,
                targetScales,
                targetLinkData[targetSubset]["before"],
                targetSubset,
                sourceSubset
              ),
              h1: calculateLinkHeight(
                clone,
                targetScales,
                targetLinkData[targetSubset]["before"],
                sourceSubset
              ),
              c0: color(sourceSubset),
              c1: color(targetSubset),
            }));

            return [...rsf3, ...rsf3Records];
          },
          []
        );

        return [...rsf2, ...rsf2Record];
      }, []);

      // const sourceCloneLinks = sourceSubsets.reduce((subsetLinksSF, subset) => {
      //   const sourceClones = sourceData["subsets"][subset]["values"];

      //   const linkRecords = sourceClones.reduce((clonesRSF, clone) => {
      //     const x0 = timeScale(sourceTimepoint) + NODE_WIDTH;
      //     const x1 = timeScale(targetTimepoint);

      //     const targetSubsetClones = targetData["clones"].hasOwnProperty(clone)
      //       ? targetData["clones"][clone]["values"]
      //       : [];

      //     const records = targetSubsetClones.map((targetSubset) => {
      //       return {
      //         id: `${sourceTimepoint}_${subset}_${targetSubset}_${clone}`,
      //         source: `${sourceTimepoint}_${subset}`,
      //         target: `${targetTimepoint}_${targetSubset}`,
      //         x0,
      //         x1,
      //         y0: calculateLinkYStart(
      //           clone,
      //           sourceScales,
      //           sourceData,
      //           targetData,
      //           subset,
      //           targetSubset
      //         ),
      //         h0: calculateLinkHeight(
      //           clone,
      //           sourceScales,
      //           sourceData,
      //           targetData,
      //           subset,
      //           targetSubset
      //         ),
      //         y1: calculateLinkYStart(
      //           clone,
      //           targetScales,
      //           targetData,
      //           sourceData,
      //           targetSubset,
      //           subset
      //         ),
      //         h1: calculateLinkHeight(
      //           clone,
      //           targetScales,
      //           targetData,
      //           sourceData,
      //           targetSubset,
      //           subset
      //         ),
      //         c0: color(subset),
      //         c1: color(targetSubset),
      //       };
      //     });

      //     return [...clonesRSF, ...records];
      //   }, []);

      //   return [...subsetLinksSF, ...linkRecords];
      // }, []);

      return [...rsf, ...rsfRecord];
    },

    []
  );

  const ref = useD3(
    (svg) => {
      const node = svg
        .append("g")
        .attr("stroke", "currentColor")
        .selectAll("rect")
        .data(nodes)
        .join("rect")
        .attr("x", (d) => d.x0)
        .attr("y", (d) => d.y0)
        .attr("height", (d) => d.height)
        .attr("width", (d) => d.width)
        .attr("fill", (d) =>
          isHighlighted(d.id, highlightedNode) ? d.color : NULL_COLOR
        )
        .on("mouseenter", (d, i, e) => {
          setHighlightedNode(d.id);
        })
        .on("mouseout", (d, i, e) => {
          setHighlightedNode(null);
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
          isHighlighted(d.source, highlightedNode) ||
          isHighlighted(d.target, highlightedNode)
            ? `url(#${uid}-link-${d.id})`
            : NULL_COLOR
        );

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

const processData = (data, groupParam, countParam) => {
  // calculate counts per group
  // summation of counts per group
  // total per group
  // list of values in group
  // list of values in each count
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
  clone,
  scales,
  sourceData,
  sourceSubset,
  targetSubset
) => {
  // Subset

  // height of subset
  // offset of target subset start within subset (cumulative of ones before it)

  // Clone

  // height of subset
  const subsetHeight = scales.y(sourceSubset);

  // offset of clone start within subset
  const cloneStartHeight = scales.height(sourceData["sums"][targetSubset]);

  // offset of specific clone within clone, % of clones from target timepoint above target subset * total cells in clone in source subset
  const cloneLinkHeight = scales.height(
    sourceData[targetSubset]["clones"][clone]["sum"]
  );

  return subsetHeight + cloneStartHeight + cloneLinkHeight;
};

const calculateLinkHeight = (clone, scales, sourceData, targetSubset) => {
  // Subset
  // number of cells in source subset
  // proportion of clone cells in target subset across all clones in source

  // Clone
  // number of cells in source subset

  // proportion of clone cells in target subset (compared to total of clone)

  return scales.height(sourceData[targetSubset]["clones"][clone]["count"]);
};

export default Sankey;
