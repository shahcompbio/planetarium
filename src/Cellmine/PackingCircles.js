import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import _ from "lodash";

import Tooltip from "../Tooltip/Tooltip";
import { useD3 } from "../utils/useD3";

const PackingCircles = ({
  data,
  width,
  height,
  radiusParam,
  idParam,
  highlightedIDs = [],
  onClick,
  tooltipFields,
}) => {
  const [nodes, setNodes] = useState(null);
  const [highlightedNode, setHighlightedNode] = useState(null);

  const radius = d3
    .scaleLinear()
    .range([10, width * 0.15])
    .domain(d3.extent(data, (d) => d[radiusParam]));

  const drawHighlightedCircles = (nodes) => {
    nodes
      .transition()
      .style("fill", "#52a2a2")
      .style("fill-opacity", 0.8)
      .attr("stroke", "black")
      .style("stroke-width", 1);
  };

  const drawUnhighlightedCircles = (nodes) => {
    nodes
      .transition()
      .style("fill", "#aecece")
      .style("fill-opacity", 0.5)
      .style("stroke-width", 0.5);
  };

  const addMouse = (nodes, allNodes) => {
    nodes
      .on("mouseover", (data) => {
        const selection = allNodes.filter(
          (node) => data[idParam] !== node[idParam]
        );

        drawUnhighlightedCircles(selection);
        setHighlightedNode(data);
      })
      .on("mouseleave", (data) => {
        drawHighlightedCircles(nodes);
        setHighlightedNode(null);
      })
      .on("mousedown", (data) => {
        onClick(data);
      });
  };

  const removeMouse = (nodes) => {
    nodes.on("mouseover", null).on("mouseleave", null).on("mousedown", null);
  };

  useEffect(() => {
    // only want to do this after initial draw
    if (nodes !== null) {
      if (highlightedIDs.length === 0) {
        drawHighlightedCircles(nodes);
        addMouse(nodes, nodes);
      } else {
        const high = nodes.filter((node) =>
          highlightedIDs.includes(node[idParam])
        );
        const dark = nodes.filter(
          (node) => !highlightedIDs.includes(node[idParam])
        );

        drawHighlightedCircles(high);
        addMouse(high, nodes);

        drawUnhighlightedCircles(dark);
        removeMouse(dark);
      }
    }
  }, [highlightedIDs, nodes === null]);
  const ref = useD3(
    (svg) => {
      requestAnimationFrame(() => {
        drawCircles(svg, data);
      });
    },
    width,
    height,
    []
  );
  const drawCircles = (svg, data) => {
    const svgNodes = svg
      .append("g")
      .selectAll("circle")
      .data(data)
      .enter()
      .append("circle")
      .attr("id", (d) => d[idParam])
      .attr("r", function (d) {
        return radius(d[radiusParam]);
      })
      .attr("cx", (d) => d.x)
      .attr("cy", (d) => d.y);

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
          .radius(function (d) {
            return radius(d[radiusParam]) + 3;
          })
          .iterations(1)
      );

    simulation.nodes(data).on("tick", () => {
      svgNodes
        .attr("cx", function (d) {
          d.x = Math.max(
            radius(d[radiusParam]) + 10,
            Math.min(width - radius(d[radiusParam]), d.x)
          );
          return d.x;
        })
        .attr("cy", function (d) {
          d.y = Math.max(
            radius(d[radiusParam]) + 10,
            Math.min(height - radius(d[radiusParam]), d.y)
          );
          return d.y;
        });
    });

    setNodes(svgNodes);
  };
  return (
    <div style={{ position: "relative" }}>
      <svg ref={ref} />
      <Tooltip
        getText={(node) => getTooltipText(node, tooltipFields)}
        getX={(node) => node["x"]}
        getY={(node) => node["y"] - radius(node[radiusParam]) + 5}
        data={highlightedNode}
      />
    </div>
  );
};

const getTooltipText = (data, tooltipFields) => {
  return (
    <div>
      {tooltipFields.map((field) => (
        <p>
          <b>{field["label"]}</b>: {data[field["param"]]}
        </p>
      ))}
    </div>
  );
};

export default PackingCircles;
