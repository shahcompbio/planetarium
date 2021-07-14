import React from "react";

import * as d3 from "d3";
import { useD3, sortAlphanumeric } from "@shahlab/planetarium";
import _ from "lodash";

const Doughnut = ({ data, colors, width, height, subsetParam }) => {
  const totalCount = data.length;
  const subsetValues = _.uniq(data.map((datum) => datum[subsetParam])).sort(
    sortAlphanumeric
  );

  const subsetCounts = _.groupBy(data, (datum) => datum[subsetParam]);
  const subsets = subsetValues.map((subset) => ({
    key: subset,
    value: subsetCounts[subset],
  }));

  const radius = Math.min(width, height) / 2;
  const arc = d3
    .arc()
    .innerRadius(radius * 0.3)
    .outerRadius(radius * 0.65);

  const arcLabel = d3
    .arc()
    .innerRadius(radius * 0.9)
    .outerRadius(radius * 0.9);

  const drawArea = (svg) => {
    const pie = d3
      .pie()
      .value((d) => d["value"].length)
      .sort(null);
    const arcs = pie(subsets);

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
              x = Math.cos(a) * 10;
              y = Math.sin(a) * 10;
              d.data._translate = { x: x, y: y };
            }

            return "translate(" + x + "," + y + ")";
          });
        d3.select("#label-" + d.data.key).style("font-weight", "bold");
      })
      .on("mouseout", function (d) {
        d3.select(this)
          .transition()
          .delay(200)
          .duration(500)
          .attr("transform", function (d) {
            d.data._expanded = false;
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

  return <svg ref={ref} />;
};

export default Doughnut;
