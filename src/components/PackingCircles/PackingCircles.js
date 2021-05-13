import React, { useState, useEffect } from "react";
import * as d3 from "d3";
//import Grid from "@material-ui/core/Grid";
//import Paper from "@material-ui/core/Paper";
import d3Tip from "d3-tip";
import { useSvg } from "../utils/useSvg";

const LEGEND_HEIGHT = 50;
const radius = 12;
const CAT_LABEL_FONT = "normal 12px Helvetica";
const PADDING = 10;
const LEGEND_SQUARE_LENGTH = 12;
const LEGEND_SQUARE_PADDING = 10;
const scale = d3.scaleOrdinal(d3.schemeCategory10);

const PackingCircles = ({ modifiedData, chartDim }) => {
  const [originalDataLength] = useState(modifiedData.length);
  const canvasWidth = chartDim["width"];
  const canvasHeight = chartDim["height"];

  const [simulation, saveSimulation] = useState();

  const chartWidth = canvasWidth - PADDING - PADDING;
  const chartHeight = canvasHeight - PADDING - PADDING;

  const [xOffSet, setXOffset] = useState(chartWidth / 2);
  const [yOffSet, setYOffset] = useState(chartHeight / 2);
  var tooltip = d3Tip()
    .attr("class", "d3-tip n")
    .attr("id", "circleTip")
    .html(
      data =>
        "Analysis Ticket: " +
        data.jira_ticket +
        "<br/>Cell Count: " +
        data.num_sublibraries +
        "<br/> Description:" +
        data.description
    )
    .offset([-10, 0]);
  const ref = useSvg(
    svgRef => {
      const svg = d3
        .select("#canvasSelection")
        .attr("viewBox", [0, 0, chartWidth, chartHeight]);

      requestAnimationFrame(() => {
        drawCircles(modifiedData, svg, chartWidth, chartHeight);
      });
    },
    canvasWidth,
    canvasHeight,
    [modifiedData]
  );
  const getSelection = modifiedData =>
    modifiedData.map(node => "#" + node["jira_ticket"]).join(",");

  useEffect(() => {
    if (modifiedData) {
      d3.select("#canvasSelection").call(tooltip);
      if (modifiedData.length === originalDataLength) {
        d3.selectAll(".node")
          .on("mouseover", (data, index, element) => {
            tooltip.show(data, element[index]).attr("opacity", 0.2);

            const selection = d3.selectAll(".node").filter(function(node) {
              return data.jira_ticket !== node.jira_ticket;
            });
            selection
              .transition()
              .style("fill", "#aecece")
              .attr("stroke-width", 1);
          })
          .on("mouseleave", function(d) {
            tooltip.hide(d, this);
            d3.selectAll(getSelection(modifiedData))
              .transition()
              .style("fill", "#307ca0");
          })
          .transition()
          .style("fill", "#307ca0");
      } else {
        const ticketText = modifiedData
          .map(node => node["jira_ticket"])
          .join(" ");

        d3.selectAll(".node")
          .on("mouseover", (data, index, element) => {
            tooltip.show(data, element[index]).attr("opacity", 0.2);

            const selection = d3.selectAll(".node").filter(function(node) {
              return data.jira_ticket !== node.jira_ticket;
            });
            selection.transition().style("fill", "#aecece");
          })
          .on("mouseleave", function(d) {
            tooltip.hide(d, this);
            d3.selectAll(getSelection(modifiedData))
              .transition()
              .style("fill", "#307ca0");
          });

        d3.selectAll(".node")
          .filter(node => ticketText.indexOf(node["jira_ticket"]) === -1)
          .on("mouseover", null)
          .on("mouseout", null)
          .transition()
          .style("fill", "#aecece")
          .attr("stroke-width", 0.5);
      }
    }
  }, [modifiedData]);

  const drawCircles = (data, svg, width, height) => {
    svg.call(tooltip);

    const radius = d3
      .scaleLinear()
      .range([12, 200])
      .domain(d3.extent(data, d => d.num_sublibraries));

    var node = svg
      .append("g")
      .selectAll("circle")
      .data(data)
      .enter()
      .append("circle")
      .attr("class", "circles")
      .attr("id", d => d.jira_ticket)
      .attr("class", "node")
      .attr("r", function(d) {
        return radius(d.num_sublibraries);
      })
      .attr("cx", d => d.x)
      .attr("cy", d => d.y)
      .style("fill", function(d) {
        return "#307ca0";
      })
      .style("fill-opacity", 0.8)
      .attr("stroke", "black")
      .style("stroke-width", 1)
      .on("mouseover", (data, index, element) => {
        tooltip.show(data, element[index]).attr("opacity", 0.2);

        const selection = d3.selectAll(".node").filter(function(node) {
          return data.jira_ticket !== node.jira_ticket;
        });
        selection.transition().style("fill", "#aecece");
      })
      .on("mouseleave", function(d) {
        tooltip.hide(d, this);
        d3.selectAll(getSelection(modifiedData))
          .transition()
          .style("fill", "#307ca0");
      });
    //.call(drag(simulation));

    const simulation = d3
      .forceSimulation()
      .alphaTarget(0.1)
      .force(
        "center",
        d3
          .forceCenter()
          .x(width / 2)
          .y(height / 2)
      )
      .force("charge", d3.forceManyBody().strength(0.1))
      .force(
        "collide",
        d3
          .forceCollide()
          .strength(0.2)
          .radius(function(d) {
            return radius(d.num_sublibraries) + 3;
          })
          .iterations(1)
      );

    simulation.nodes(data).on("tick", function(d) {
      node
        .attr("cx", function(d) {
          return d.x;
        })
        .attr("cy", function(d) {
          return d.y;
        });
    });
    saveSimulation(simulation);
  };

  const drag = simulation => {
    function dragstarted(d) {
      if (!d3.event.active) simulation.alphaTarget(0.2).restart();
      d.fx = d.x;
      d.fy = d.y;
    }

    function dragged(d) {
      d.fx = d3.event.x;
      d.fy = d3.event.y;
    }

    function dragended(d) {
      d3.event.sourceEvent.stopPropagation();
      if (!d3.event.active) simulation.alphaTarget(0);
      d.fx = null;
      d.fy = null;
    }

    return d3
      .drag()
      .on("start", dragstarted)
      .on("drag", dragged)
      .on("end", dragended);
  };
  return (
    <div
      style={{
        width: chartWidth,
        height: chartHeight,
        position: "relative"
      }}
      item
    >
      <svg
        id="canvasSelection"
        ref={ref}
        style={{
          background: "rgb(190 214 214)",
          width: chartWidth,
          height: chartHeight,
          position: "relative"
        }}
      />
    </div>
  );
};

export default PackingCircles;
