import React, { useState } from "react";
import _ from "lodash";

import * as d3 from "d3";
import { useD3 } from "@shahlab/planetarium";
import { Grid } from "@material-ui/core";

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

const Doughnut = ({ data, colors, totalCount, width, height }) => {
  const chartWidth = width;
  const chartHeight = height;
  const radius = Math.min(width, height) / 2;
  const arc = d3
    .arc()
    .innerRadius(radius * 0.4)
    .outerRadius(radius * 0.65);

  const arcLabel = d3
    .arc()
    .innerRadius(radius)
    .outerRadius(radius * 0.8);

  const drawArea = (svg) => {
    const pie = d3
      .pie()
      .value((d) => d["value"].length)
      .sort(null);
    const arcs = pie(data);

    const path = svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("stroke", "white")
      .selectAll("path")
      .data(arcs);

    path
      .join("path")
      .attr("fill", (d) => colors(d["data"].key))
      .attr("d", arc)
      .on("mouseover", function (d) {
        d3.select(this)
          .transition()
          .duration(500)
          .attr("transform", function (d) {
            var x;
            var y;
            if (d.data._translate) {
              x = d.data._translate.x;
              y = d.data._translate.y;
            } else if (!d.data._expanded) {
              d.data._expanded = true;
              var a =
                d.startAngle + (d.endAngle - d.startAngle) / 2 - Math.PI / 2;
              x = Math.cos(a) * 15;
              y = Math.sin(a) * 15;
              d.data._translate = { x: x, y: y };
            }

            return "translate(" + x + "," + y + ")";
          });
        d3.select("#label-" + d.data.key).style("font-weight", "bold");
      })
      .on("mouseout", function (d) {
        d3.select(this)
          .transition()
          .duration(500)
          .attr("transform", function (d) {
            d.data._expanded = false;

            return (
              "translate(-" + d.data._translate + ",-" + d.data._translate + ")"
            );
          });
        d3.select("#label-" + d.data.key).style("font-weight", "normal");
      });

    svg
      .append("g")
      .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
      .attr("font-family", "sans-serif")
      .attr("font-size", 12)
      .attr("text-anchor", "middle")
      .selectAll("text")
      .data(arcs)
      .join("text")
      .attr("id", (d) => "label-" + d.data.key)
      .attr("transform", (d, i) => {
        return `translate(${arcLabel.centroid(d, i)})`;
      })
      .call((text) =>
        text
          .append("tspan")
          .attr("y", "-0.4em")
          .attr("font-weight", "bold")
          .text((d) => d.data.name)
      )
      .call((text) =>
        text
          .filter((d) => d.endAngle - d.startAngle > 0.25)
          .append("tspan")
          //  .attr("x", 0)
          //  .attr("y", "-0.7em")
          .attr("fill-opacity", 0.7)
          .text(
            (d) =>
              d.data.key.toLocaleString() +
              " (" +
              d3.format(".0%")(d.data.value.length / totalCount) +
              ")"
          )
      );
  };

  const ref = useD3(
    (svg) => {
      drawArea(svg);
    },
    width,
    height,
    [data]
  );

  /*const setHighlighted = (event, value) => {
    if (event === "mouseenter") {
      setHighlightedSubset(value);
    } else if (event === "mousedown") {
      // setHighlightedSubset(value);
    } else if (event === "mouseout") {
      setHighlightedSubset(null);
    }
  };*/

  return (
    <Grid container direction="row" style={{ padding: 0 }}>
      <Grid item>
        <svg ref={ref} />
      </Grid>
    </Grid>
  );
};

export default Doughnut;
