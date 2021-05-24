import React, { useState } from "react";
import * as d3 from "d3";
import * as d3Chord from "d3-chord";
import * as d3Array from "d3-array";
import { getData, strokeToFill } from "gradient-path";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useSvg } from "../utils/useSvg";

const PROP_AXIS_FONT = "normal 10px Helvetica";
const CAT_LABEL_FONT = "normal 12px Helvetica";
const PADDING = 10;
const scale = d3.scaleOrdinal(d3.schemeCategory10);
const labelMapping = { l: "Left Tumor", r: "Right Tumor", s: "Spleen" };
const Network = ({ data, chartDim, chartName }) => {
  const canvasWidth = chartDim["width"];
  const canvasHeight = chartDim["height"];

  const chartWidth = canvasWidth - PADDING - PADDING;
  const chartHeight = canvasHeight - PADDING - PADDING;

  const outerRadius = Math.min(chartWidth, chartHeight) * 0.5 - 60;
  const innerRadius = outerRadius - 10;
  const ref = useSvg(
    svgRef => {
      const svg = d3
        .select("#canvasSelection")
        .attr("viewBox", [0, 0, chartWidth, chartHeight]);

      const svgTwo = d3
        .select("#canvasSelectionTwo")
        .attr("viewBox", [0, 0, chartWidth, chartHeight]);

      const modifiedData = data.filter(
        row => row["clonotype"] !== row["CDR3 AA"]
      );
      requestAnimationFrame(() => {
        drawChord(
          svg,
          modifiedData,
          innerRadius,
          outerRadius,
          canvasWidth,
          canvasHeight,
          "clonotype"
        );
        drawChord(
          svgTwo,
          modifiedData,
          innerRadius,
          outerRadius,
          canvasWidth,
          canvasHeight,
          "CDR3 AA"
        );
      });
    },
    canvasWidth,
    canvasHeight,
    [data]
  );

  return (
    <Grid container direction="row" justify="flex-start" alignItems="stretch">
      <Paper
        style={{
          margin: 10,
          padding: PADDING,
          height: chartDim["height"],
          width: chartDim["width"]
        }}
      >
        <Grid
          container
          direction="column"
          justify="flex-start"
          alignItems="stretch"
        >
          <Grid
            style={{
              width: chartWidth,
              height: chartHeight,
              position: "relative",
              background: "#323333"
            }}
            item
          >
            <svg
              id="canvasSelection"
              ref={ref}
              style={{
                width: chartWidth,
                height: chartHeight,
                position: "relative",
                background: "aliceblue"
              }}
            />
          </Grid>
        </Grid>
      </Paper>
      <Paper
        style={{
          margin: 10,
          padding: PADDING,
          height: chartDim["height"],
          width: chartDim["width"]
        }}
      >
        <Grid
          container
          direction="column"
          justify="flex-start"
          alignItems="stretch"
        >
          <Grid
            style={{
              width: chartWidth,
              height: chartHeight,
              position: "relative",
              background: "#323333"
            }}
            item
          >
            <svg
              id="canvasSelectionTwo"
              ref={ref}
              style={{
                width: chartWidth,
                height: chartHeight,
                position: "relative",
                background: "white"
              }}
            />
          </Grid>
        </Grid>
      </Paper>
    </Grid>
  );
};
const getFreq = list =>
  list.reduce((final, row) => {
    final[row["value"]] = final.hasOwnProperty(row["value"])
      ? final[row["value"]] + row["freq"]
      : row["freq"];
    return final;
  }, {});

const getFrequencySum = (matches, fullList) =>
  d3.sum([...matches].map(match => fullList[match]));
const getDifference = (a, b) =>
  d3Array.difference(Object.keys(a), Object.keys(b));

const getIntersection = (a, b) =>
  d3Array.intersection(Object.keys(a), Object.keys(b));
const getSetIntersection = (a, b) => d3Array.intersection([...a], [...b]);
const getCorrectSite = (siteAndClonotype, site) =>
  siteAndClonotype[0][site]
    ? siteAndClonotype[0][site]
    : siteAndClonotype[1][site]
    ? siteAndClonotype[1][site]
    : siteAndClonotype[2][site];

const drawChord = (
  svg,
  data,
  innerRadius,
  outerRadius,
  width,
  height,
  parameter
) => {
  const colors = [
    "#69b40f",
    "#ec1d25",
    "#008fc8",
    "#10218b",
    "#134b24",
    "#737373"
  ];
  console.log([...new Set(data.map(d => d.mouse))]);
  const mouse = "10";
  const modifiedData = data.filter(d => d.mouse === mouse);
  console.log(parameter);
  //45
  //27

  const conditions = Array.from(
    d3Array.group(modifiedData, d => d.condition, d => d.site)
  );
  console.log(conditions);
  const names = Array.from(conditions[0]);
  const allSites = Array.from(names[1]);

  const siteAndClonotype = allSites.map((site, index) => ({
    [site[0]]: site[1].map(row => ({
      value: row[parameter],
      freq: parseFloat(row["Frequency"])
    }))
  }));

  svg
    .append("g")
    .attr("class", "title")
    .append("text")
    .attr("font-size", 20)
    .attr("font-weight", "bold")
    .attr("x", 10)
    .attr("y", 25)
    .text(names[0] + " - " + mouse);
  console.log(siteAndClonotype);
  const s = getFreq(Array.from(getCorrectSite(siteAndClonotype, "s")));
  console.log(s);
  const r = getFreq(Array.from(getCorrectSite(siteAndClonotype, "r")));
  const l = getFreq(Array.from(getCorrectSite(siteAndClonotype, "l")));

  const sr = getFrequencySum(getIntersection(s, r), s);
  const rs = getFrequencySum(getIntersection(s, r), r);

  const sl = getFrequencySum(getIntersection(s, l), s);
  const ls = getFrequencySum(getIntersection(s, l), l);

  const rl = getFrequencySum(getIntersection(r, l), r);
  const lr = getFrequencySum(getIntersection(r, l), l);

  const sSize = getFrequencySum(
    getSetIntersection(getDifference(s, r), getDifference(s, l)),
    s
  );
  const rSize = getFrequencySum(
    getSetIntersection(getDifference(r, s), getDifference(r, l)),
    r
  );
  const lSize = getFrequencySum(
    getSetIntersection(getDifference(l, s), getDifference(l, r)),
    l
  );
  console.log(rSize);
  console.log(sr);
  console.log(rs);
  console.log(rl);
  //   s r l
  // s (1,1),(1,2),(1,3)
  // r (2,1),(2,2),(2,3)
  // l (3,1),(3,2),(3,3)
  const matrix = [[sSize, sr, sl], [rs, rSize, rl], [ls, lr, lSize]];

  //const chordIntersections = allSites.map((site,index)=>{
  //  site[1]
  //})

  const chord = d3Chord
    .chord()
    .padAngle(10 / innerRadius)
    .sortSubgroups(d3.ascending)
    .sortChords(d3.ascending);
  const arc = d3
    .arc()
    .innerRadius(innerRadius)
    .outerRadius(outerRadius);

  const ribbon = d3Chord
    .ribbon()
    .radius(innerRadius - 1)
    .padAngle(1 / innerRadius);

  const tickStep = d3.tickStep(
    0,
    d3.sum(data.map(row => row["Frequency"]).flat()),
    100
  );

  const formatValue = d3.format(".1~%");
  function ticks({ startAngle, endAngle, value }) {
    const k = (endAngle - startAngle) / value;

    return d3.range(0, value, tickStep).map(value => {
      return { value, angle: value * k + startAngle };
    });
  }

  const color = d3.scaleOrdinal(allSites.map(site => site[0]), colors);

  const chords = chord(matrix);
  //console.log(chords);
  const group = svg
    .append("g")
    .attr("class", "chordNodes")
    .attr("font-size", 10)
    .attr("font-family", "sans-serif")
    .selectAll("g")
    .data(chords.groups)
    .join("g");
  const name = ["s", "r", "l"];
  group
    .append("path")
    .attr("fill", "grey")
    //  .attr("fill", d => color(name[d.index]))
    .attr("d", arc);

  const groupTick = group
    .append("g")
    .selectAll("g")
    .data(ticks)
    .join("g")
    .attr("transform", d => {
      return `rotate(${(d.angle * 180) / Math.PI -
        90}) translate(${outerRadius},0)`;
    });

  groupTick
    .enter()
    .append("line")
    .attr("stroke", "currentColor")
    .attr("x2", 6);

  groupTick
    .append("text")
    .attr("x", 8)
    .attr("dy", "0.35em")
    .attr("fill", "black")
    .attr("transform", d =>
      d.angle > Math.PI ? "rotate(180) translate(-16)" : null
    )
    .attr("text-anchor", d => (d.angle > Math.PI ? "end" : null))
    .text(d => d.value);

  group
    .select("text")
    .attr("font-weight", "bold")
    .attr("font-size", 11)
    .text(function(d) {
      return this.getAttribute("text-anchor") === "end"
        ? `↑ ${labelMapping[name[d.index]]}`
        : `${labelMapping[name[d.index]]} ↓`;
    });
  const colored = {
    rl: "#69b40f",
    lr: "#69b40f",
    rs: "#ec1d25",
    sr: "#ec1d25",
    sl: "#008fc8",
    ls: "#008fc8"
  };

  svg
    .append("g")
    .attr("fill-opacity", 0.8)
    .attr("class", "chords")
    .selectAll("path")
    .data(chords)
    .join("path")
    .style("mix-blend-mode", "multiply")
    .attr("fill", d => {
      return colored[name[d.source.index] + name[d.target.index]];
    })
    //  .attr("fill", "grey")
    //.attr("fill", d => color(name[d.source.index]))
    .attr("d", ribbon)
    .attr("opacity", d => (d.source.index === d.target.index ? 0 : 0.8))
    .append("title")
    .text(d => {
      return `${d.source.value} ${name[d.target.index]} → ${
        name[d.source.index]
      }${
        d.source.index === d.target.index
          ? ""
          : `\n${d.target.value} ${name[d.source.index]} → ${
              name[d.target.index]
            }`
      }}`;
    });

  d3.selectAll(".chordNodes").attr(
    "transform",
    "translate(" + width / 2 + "," + height / 2 + ")"
  );
  d3.selectAll(".chords").attr(
    "transform",
    "translate(" + width / 2 + "," + height / 2 + ")"
  );
};

export default Network;
