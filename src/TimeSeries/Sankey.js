// Filter data for clones in both pre / post

// Aggregate data by timepoint, clone, phenotype - we want a count

import React, { useState } from "react";
import * as d3 from "d3";

import { curveBumpX } from "d3-shape";

import _ from "lodash";
import { useD3 } from "../utils/useD3";
import sortAlphanumeric from "../utils/sortAlphanumeric";

const PADDING = 10;

const NODE_WIDTH = 15;
const NODE_VERTICAL_PADDING = 10;

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
  const chartHeight = height;

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
    .range([PADDING, PADDING + chartWidth]);

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
    }));

    return [...rsf, ...records];
  }, []);

  // links need to iterate over pairs
  // - matched by clones (then subsetted by subset)
  // - needs to be in adjacent timepoints

  const timepointsButLast =
    timepoints.length > 2
      ? timepoints[(0, timepoints.length - 2)]
      : [timepoints[0]];

  const links = timepointsButLast.reduce(
    (rsf, sourceTimepoint, index) => {
      const sourceData = counts[sourceTimepoint];
      const sourceScales = subsetScales[sourceTimepoint];

      const targetTimepoint = timepoints[index + 1];
      const targetData = counts[targetTimepoint];
      const targetScales = subsetScales[targetTimepoint];
      // for each subset-clone pair, generate a list of links

      const sourceSubsets = sourceData["subsetValues"];
      const sourceCloneLinks = sourceSubsets.reduce((subsetLinksSF, subset) => {
        const sourceClones = sourceData["subsets"][subset]["values"];

        const linkRecords = sourceClones.reduce((clonesRSF, clone) => {
          const x0 = timeScale(sourceTimepoint) + NODE_WIDTH;
          const x1 = timeScale(targetTimepoint) - NODE_WIDTH;

          // for a subset-clone pair, iterate over every instance of clone in target timepoint

          const targetSubsetClones = targetData["clones"].hasOwnProperty(clone)
            ? targetData["clones"][clone]["values"]
            : [];

          const records = targetSubsetClones.map((targetSubset) => {
            // !!! offset needs to account that it might not start at exactly subset y0
            // so i need to know both how clone count increases within a subset
            // but also something about knowing proportions in the opposite timepoint? oy veyv

            //  ineed to keep track of same clone different subset offset?
            return {
              id: `${sourceTimepoint}_${subset}_${targetSubset}_${clone}`,
              x0,
              x1,
              y0: calculateLinkYStart(
                clone,
                sourceScales,
                sourceData,
                targetData,
                subset,
                targetSubset
              ),
              h0: calculateLinkHeight(
                clone,
                sourceScales,
                sourceData,
                targetData,
                subset,
                targetSubset
              ),
              y1: calculateLinkYStart(
                clone,
                targetScales,
                targetData,
                sourceData,
                targetSubset,
                subset
              ),
              h1: calculateLinkHeight(
                clone,
                targetScales,
                targetData,
                sourceData,
                targetSubset,
                subset
              ),
              c0: color(subset),
              c1: color(targetSubset),
            };
          });

          return [...clonesRSF, ...records];
        }, []);

        return [...subsetLinksSF, ...linkRecords];
      }, []);

      return [...rsf, ...sourceCloneLinks];
    },

    []
  );
  // to draw links, needs to know:
  // - source timepoint (x0 position)
  // - target timepoint (x1 position)
  // - source subset, clone, cell count (y starting positions)
  // - target subset, clone, cell count (y ending positions)

  const ref = useD3(
    (svg) => {
      const color = colorScale || d3.scaleOrdinal(nodes, d3.schemeTableau10);

      const node = svg
        .append("g")
        .selectAll("rect")
        .data(nodes)
        .join("rect")
        .attr("x", (d) => d.x0)
        .attr("y", (d) => d.y0)
        .attr("height", (d) => d.height)
        .attr("width", (d) => d.width)
        .attr("fill", (d) => d.color);

      const link = svg
        .append("g")
        .attr("fill", "none")
        .attr("stroke-opacity", 0.5)
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
            [d.x0 + 5, d.y0, d.y0 + d.h0],
            [d.x1, d.y1, d.y1 + d.h1],
          ])
        )
        .attr("fill", (d) => `url(#${uid}-link-${d.id})`);
      // .attr("stroke", ({ index: i }) => `url(#${uid}-link-${i})`)
      // .attr("stroke-width", ({ width }) => Math.max(1, width));
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

const calculateLinkYStart = (
  clone,
  scales,
  sourceData,
  targetData,
  sourceSubset,
  targetSubset
) => {
  // height of subset
  const subsetHeight = scales.y(sourceSubset);

  // offset of clone start within subset
  const cloneStartHeight = scales.height(
    sourceData["subsets"][sourceSubset]["sum"][clone]
  );

  // offset of specific clone within clone, % of clones from target timepoint above target subset * total cells in clone in source subset
  const cloneLinkHeight = scales.height(
    (sourceData["subsets"][sourceSubset]["counts"][clone] *
      targetData["clones"][clone]["sum"][targetSubset]) /
      targetData["clones"][clone]["total"]
  );

  return subsetHeight + cloneStartHeight + cloneLinkHeight;
};

const calculateLinkHeight = (
  clone,
  scales,
  sourceData,
  targetData,
  sourceSubset,
  targetSubset
) => {
  // number of cells in source subset

  // proportion of clone cells in target subset (compared to total of clone)

  return scales.height(
    (sourceData["subsets"][sourceSubset]["counts"][clone] *
      targetData["clones"][clone]["counts"][targetSubset]) /
      targetData["clones"][clone]["total"]
  );
};

export default Sankey;
