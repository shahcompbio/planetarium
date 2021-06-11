import React, { useState, useEffect } from "react";
import * as d3 from "d3";
import _ from "lodash";

import d3Tip from "d3-tip";
import { useSvg } from "../utils/useSvg";
import { useD3 } from "../utils/useD3";
import { isValueHighlighted as isHighlighted } from "../utils/isHighlighted";

const PADDING = 10;

const PackingCircles = ({
  data,
  width,
  height,
  radiusParam,
  idParam,
  highlightedIDs = [],
  selectAnalysis,
  tooltipFields,
}) => {
  const [nodes, setNodes] = useState(null);
  const tooltip = d3Tip()
    .attr("class", "d3-tip n")
    .attr("id", "circleTip")
    .html(
      (data) => "hi"
      // tooltipFields
      //   .map(
      //     (field) =>
      //       `<span style='font-weight:bold'>${field["label"]}</span> <span>${
      //         data[field["param"]]
      //       }</span>`
      //   )
      //   .reduce((str, line) => `${str}</br>${line}`)
    )
    .offset([-10, 0]);
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
      .on("mouseover", (data, index, element) => {
        // tooltip.show(data, element[index]).attr("opacity", 0.2);
        const selection = allNodes.filter(
          (node) => data[idParam] !== node[idParam]
        );

        drawUnhighlightedCircles(selection);
      })
      .on("mouseleave", (data) => {
        // tooltip.hide(data, this);
        drawHighlightedCircles(nodes);
      });
  };

  const removeMouse = (nodes) => {
    nodes.on("mouseover", null).on("mouseleave", null);
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
    const radius = d3
      .scaleLinear()
      .range([10, width * 0.15])
      .domain(d3.extent(data, (d) => d[radiusParam]));

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

    svg.call(tooltip);
    setNodes(svgNodes);
  };
  return (
    <div>
      <svg ref={ref} />
    </div>
  );
};

// const PackingCircles = ({ data, width, height, selectAnalysis }) => {
//   const [originalDataLength] = useState(data.length);

//   const [simulation, saveSimulation] = useState();

//   const chartWidth = width - PADDING - PADDING;
//   const chartHeight = height - PADDING - PADDING;

//   var tooltip = d3Tip()
//     .attr("class", "d3-tip n")
//     .attr("id", "circleTip")
//     .html(
//       (data) =>
//         "<span style='font-weight:bold'>Analysis Ticket:</span> <span>" +
//         data.jira_ticket +
//         "</span><br/><span style='font-weight:bold'>Cell Count:</span> <span>" +
//         data.num_sublibraries +
//         "</span><br/><span style='font-weight:bold'>Description:</span> <span>" +
//         data.description +
//         "</span>"
//     )
//     .offset([-10, 0]);

//   const getSelection = (data) =>
//     data.map((node) => "#" + node["jira_ticket"]).join(",");

//   useEffect(() => {
//     if (data) {
//       d3.select("#canvasSelection").call(tooltip);
//       if (data.length === originalDataLength) {
//         d3.selectAll(".node")
//           .style("fill", "#52a2a2")
//           .style("fill-opacity", 0.8)
//           .attr("stroke", "black")
//           .style("stroke-width", 1)
//           .on("mouseover", (data, index, element) => {
//             tooltip.show(data, element[index]).attr("opacity", 0.2);
//             const selection = d3.selectAll(".node").filter(function (node) {
//               return data.jira_ticket !== node.jira_ticket;
//             });
//             selection
//               .transition()
//               .style("fill", "#aecece")
//               .style("fill-opacity", 0.5)
//               //  .attr("stroke-width", 1)
//               //    .attr("stroke", "black")
//               .style("stroke-width", 0.5);

//             d3.select(element[index])
//               //.attr("stroke", "#8dbbb9")
//               .style("stroke-width", 2);
//           })
//           .on("mouseleave", function (d) {
//             tooltip.hide(d, this);
//             d3.selectAll(getSelection(data))
//               .transition()
//               .style("fill", "#52a2a2")
//               .style("fill-opacity", 0.8)
//               .attr("stroke", "black")
//               .style("stroke-width", 1);
//           })
//           .on("mousedown", function (d) {
//             selectAnalysis(d);
//           });
//       } else {
//         const ticketText = data.map((node) => node["jira_ticket"]).join(" ");

//         const nonSelectionOnModified = d3
//           .selectAll(".node")
//           .filter((node) => ticketText.indexOf(node["jira_ticket"]) === -1);

//         d3.selectAll(".node")
//           .on("mouseover", (data, index, element) => {
//             if (ticketText.indexOf(data["jira_ticket"]) !== -1) {
//               tooltip.show(data, element[index]).attr("opacity", 0.2);

//               /*  const nonSelection = d3.selectAll(".node").filter(function (node) {
//               return data.jira_ticket !== node.jira_ticket;
//             });
//           nonSelection
//               .transition()
//               .style("fill", "#aecece")
//               .style("fill-opacity", 0.5)
//               //.attr("stroke", "black")
//               .style("stroke-width", 0);
// */
//               d3.select(element[index]).style("stroke-width", 2);
//             }
//           })
//           .on("mouseleave", function (d) {
//             if (ticketText.indexOf(d["jira_ticket"]) !== -1) {
//               tooltip.hide(d, this);
//               d3.selectAll(".node")
//                 .filter(function (node) {
//                   return d.jira_ticket === node.jira_ticket;
//                 })
//                 .style("fill-opacity", 0.8)
//                 .attr("stroke", "black")
//                 .style("stroke-width", 1);
//             }
//             /*  d3.selectAll(getSelection(data))
//               .transition()
//               .style("fill", "#52a2a2")
//               .style("fill-opacity", 0.8)
//               .attr("stroke", "black")
//               .style("stroke-width", 1);*/
//           })
//           .on("mousedown", function (d) {
//             selectAnalysis(d);
//           })
//           .attr("stroke-width", function (d) {
//             if (ticketText.indexOf(d["jira_ticket"]) !== -1) {
//               return 2;
//             } else {
//               return 0;
//             }
//           });
//         nonSelectionOnModified
//           .style("stroke-width", 0)
//           .style("fill", "#aecece")
//           .style("fill-opacity", 0.5);

//         nonSelectionOnModified.on("mouseover", null).on("mouseout", null);
//       }
//     }
//   }, [data]);

//   const drawCircles = (data, svg, width, height) => {
//     svg.call(tooltip);

//     const radius = d3
//       .scaleLinear()
//       .range([12, 175])
//       .domain(d3.extent(data, (d) => d.num_sublibraries));

//     var node = svg
//       .append("g")
//       .selectAll("circle")
//       .data(data)
//       .enter()
//       .append("circle")
//       .attr("class", "circles")
//       .attr("id", (d) => d.jira_ticket)
//       .attr("class", "node")
//       .attr("r", function (d) {
//         return radius(d.num_sublibraries);
//       })
//       .attr("cx", (d) => d.x)
//       .attr("cy", (d) => d.y)
//       .style("fill", "#52a2a2")
//       .style("fill-opacity", 0.8)
//       .attr("stroke", "black")
//       .style("stroke-width", 1)
//       .on("mouseover", (data, index, element) => {
//         tooltip.show(data, element[index]).attr("opacity", 0.2);

//         const nonSelection = d3.selectAll(".node").filter(function (node) {
//           return data.jira_ticket !== node.jira_ticket;
//         });
//         nonSelection
//           .transition()
//           .style("fill", "#aecece")
//           .style("fill-opacity", 0.5)
//           .style("stroke-width", 0);
//         d3.selectAll(element[index])
//           //  .filter(function (node) {
//           //    return data.jira_ticket === node.jira_ticket;
//           //  })

//           //  .attr("stroke", "#8dbbb9")
//           .style("stroke-width", 2);
//       })
//       .on("mouseleave", function (d) {
//         tooltip.hide(d, this);
//         d3.selectAll(getSelection(data))
//           .transition()
//           .style("fill", "#52a2a2")
//           .attr("stroke", "black")
//           .style("fill-opacity", 0.8)
//           .style("stroke-width", 1);
//       })
//       .on("mousedown", function (d) {
//         selectAnalysis(d);
//       });
//     //.call(drag(simulation));

//     const simulation = d3
//       .forceSimulation()
//       .alphaTarget(0.1)
//       .force(
//         "center",
//         d3
//           .forceCenter()
//           .x(width / 2)
//           .y(height / 2)
//       )
//       .force("charge", d3.forceManyBody().strength(0.1))
//       .force(
//         "collide",
//         d3
//           .forceCollide()
//           .strength(0.2)
//           .radius(function (d) {
//             return radius(d.num_sublibraries) + 3;
//           })
//           .iterations(1)
//       );

//     simulation.nodes(data).on("tick", function (d) {
//       node
//         .attr("cx", function (d) {
//           d.x = Math.max(
//             radius(d.num_sublibraries) + 10,
//             Math.min(width - radius(d.num_sublibraries), d.x)
//           );
//           return d.x;
//         })
//         .attr("cy", function (d) {
//           d.y = Math.max(
//             radius(d.num_sublibraries) + 10,
//             Math.min(height - radius(d.num_sublibraries), d.y)
//           );
//           return d.y;
//         });
//       /*.attr("cx", function (d) {
//           return d.x;
//         })
//         .attr("cy", function (d) {
//           return d.y;
//         });*/
//     });
//     saveSimulation(simulation);
//   };

//   const drag = (simulation) => {
//     function dragstarted(d) {
//       if (!d3.event.active) simulation.alphaTarget(0.2).restart();
//       d.fx = d.x;
//       d.fy = d.y;
//     }

//     function dragged(d) {
//       d.fx = d3.event.x;
//       d.fy = d3.event.y;
//     }

//     function dragended(d) {
//       d3.event.sourceEvent.stopPropagation();
//       if (!d3.event.active) simulation.alphaTarget(0);
//       d.fx = null;
//       d.fy = null;
//     }

//     return d3
//       .drag()
//       .on("start", dragstarted)
//       .on("drag", dragged)
//       .on("end", dragended);
//   };
//   return (
//     <div
//       style={{
//         width: chartWidth,
//         height: chartHeight,
//         position: "relative",
//       }}
//       item
//     >
//       <svg
//         id="canvasSelection"
//         ref={ref}
//         style={{
//           background: "#586773",
//           width: chartWidth,
//           height: chartHeight,
//           position: "relative",
//         }}
//       />
//     </div>
//   );
// };

export default PackingCircles;
