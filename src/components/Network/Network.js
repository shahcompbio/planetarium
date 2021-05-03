import React, { useState } from "react";
import * as d3 from "d3";
import { getData, strokeToFill } from "gradient-path";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useSvg } from "../utils/useSvg";
import { isHighlighted } from "../utils/isHighlighted";

import infoText from "../../Info/InfoText";

const LEGEND_HEIGHT = 50;
const radius = 12;
const PROP_AXIS_FONT = "normal 10px Helvetica";
const CAT_LABEL_SPACE = 150;
const CAT_LABEL_FONT = "normal 12px Helvetica";
const PADDING = 10;
const TITLE_HEIGHT = 30;
const LABEL_PADDING = 5;

const LEGEND_SQUARE_LENGTH = 12;
const LEGEND_SQUARE_PADDING = 10;
const scale = d3.scaleOrdinal(d3.schemeCategory10);

const Network = ({ data, chartDim, chartName }) => {
  const canvasWidth = chartDim["width"];
  const canvasHeight = chartDim["height"];

  //const [simulation, saveSimulation] = useState();

  const chartWidth = canvasWidth - PADDING - PADDING;
  const chartHeight = canvasHeight - PADDING - PADDING;

  const [xOffSet, setXOffset] = useState(chartWidth / 2);
  const [yOffSet, setYOffset] = useState(chartHeight / 2);

  const ref = useSvg(
    svgRef => {
      const svg = d3
        .select("#canvasSelection")
        .attr("viewBox", [0, 0, chartWidth, chartHeight]);

      const nodeColouring = d3
        .scaleOrdinal()
        .domain(["ligand", "receptor"])
        .range(["#dcc6e0", "#e9d460"]);

      const linkColour = d3
        .scaleLinear()
        .domain([1, 25])
        .range(["#bfbfbf", "#000000"]);

      requestAnimationFrame(() => {
        drawNetwork(
          svg,
          data,
          nodeColouring,
          xOffSet,
          yOffSet,
          chartWidth,
          chartHeight,
          linkColour
        );
        drawLabels(
          svg,
          chartWidth + PADDING + 2,
          chartHeight + LEGEND_HEIGHT + LABEL_PADDING
        );
      });
    },
    canvasWidth,
    canvasHeight,
    [data]
  );

  return (
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
              position: "relative"
            }}
          />
        </Grid>
      </Grid>
    </Paper>
  );
};
const getInteractionAtIndex = (row, index) =>
  row["interaction"].split("-")[index];

const drawNetwork = (
  svg,
  data,
  nodeColouring,
  xOffSet,
  yOffSet,
  width,
  height,
  linkColour
) => {
  // calculate fixed radial amount
  const radial = Math.min(width, height) / 2 - 2 * radius;

  /*const scales = {
    color: d3.scaleSequential(d3.interpolateYlGnBu),
    theta: d3.scaleOrdinal().range([0, 2 * Math.PI])
  };*/

  const scales = {
    colour: d3.scaleSequential(d3.interpolateYlOrRd),
    radius: d3.scaleSqrt().range([2, 20]),
    stroke: d3.scaleSqrt().range([2, 8])
  };

  const sortedData = data.sort((a, b) =>
    d3.ascending(parseFloat(a.value), parseFloat(b.value))
  );
  const nodes = [
    ...new Set(
      sortedData
        .reduce(
          (final, row) => [
            ...final,
            getInteractionAtIndex(row, 0),
            getInteractionAtIndex(row, 1)
          ],
          []
        )
        .flat(1)
    )
  ].reduce((final, row, i) => [...final, { id: row }], []);
  // setup color scale using closeness

  const closeness = d3.extent(data, v => v.value);

  //  scales.theta.domain([0, nodes.length]);
  //scales.color.domain(closeness);
  console.log(nodes);
  const lengthScale = d3
    .scaleLinear()
    .domain(closeness)
    .range([10, 100]);

  const totalLinkWeight = sortedData.reduce((final, row) => {
    const source = getInteractionAtIndex(row, 0);
    const target = getInteractionAtIndex(row, 1);

    final[source] = final[source]
      ? {
          value: final[source]["value"] + parseFloat(row["value"]),
          count: final[source]["count"] + 1
        }
      : { value: parseFloat(row["value"]), count: 1 };
    final[target] = final[target]
      ? {
          value: final[target]["value"] + parseFloat(row["value"]),
          count: final[target]["count"] + 1
        }
      : { value: parseFloat(row["value"]), count: 1 };
    return final;
  }, {});

  console.log(totalLinkWeight);
  const links = sortedData.reduce((final, row) => {
    return [
      ...final,
      {
        id: row.interaction,
        source: getInteractionAtIndex(row, 0),
        target: getInteractionAtIndex(row, 1),
        value: parseFloat(row["value"])
      }
    ];
  }, []);

  scales.colour.domain(
    d3.extent(Object.keys(totalLinkWeight), v => {
      //  console.log(totalLinkWeight[v]["count"]);
      //  console.log(totalLinkWeight[v]["value"] / totalLinkWeight[v]["count"]);
      return totalLinkWeight[v]["value"] / totalLinkWeight[v]["count"];
    })
  );
  scales.radius.domain(d3.extent(links, v => v.value));
  scales.stroke.domain(d3.extent(links, e => e.value));
  console.log(scales.colour.domain());
  const simulation = d3
    .forceSimulation()
    .force("center", d3.forceCenter(xOffSet, yOffSet))

    .force("x", d3.forceX().strength(0.2))
    .force("y", d3.forceY().strength(0.2))
    //  .force("forceX", d3.forceX())
    //  .force("forceY", d3.forceY())
    .force(
      "collide",
      d3
        .forceCollide()
        .strength(10)
        .radius(v => scales.radius(v.value) + 2)
        .iterations(15)
    )
    .force("charge", d3.forceManyBody().strength(-700))
    .force("link", d3.forceLink());
  //  simulation.stop();
  console.log(links);
  console.log(d3.extent(links, v => v.value));

  /*const simulation = d3
    .forceSimulation()
    //  .force("collide", d3.forceCollide(25))
    .force(
      "charge",
      d3
        .forceManyBody()
        .strength(-700)
        .distanceMin(50)
        .distanceMax(100)
    )
    .velocityDecay(0.4)
    .alphaDecay(0.03)
    .force("center", d3.forceCenter(xOffSet, yOffSet))
    .force("x", d3.forceX().strength(0.002))
    .force("y", d3.forceY().strength(0.002));
  /*  .force("charge", d3.forceManyBody())
    .force("center", d3.forceCenter(xOffSet, yOffSet))
    .velocityDecay(0.2)
    .alphaDecay(0.02)
    //.alphaTarget(0.9)
    .force("x", d3.forceX().strength(0.002))
    .force("y", d3.forceY().strength(0.002));
*/

  console.log(links);
  console.log(nodes);

  simulation.nodes(nodes).force(
    "link",
    d3
      .forceLink(links)
      .id(d => d.id)
      .strength(2)
      .distance(d => {
        return 10 * scales.radius(d.value);
      })
  );

  var zoom = d3
    .zoom()
    .scaleExtent([0.05, 1])
    .on("zoom", zoomed);

  function zoomed() {
    svg.attr("transform", d3.event.transform);
  }

  //  svg.call(zoom);

  const mainG = svg.append("g");
  const defs = svg.append("defs");

  // Sample the SVG path uniformly with the specified precision.
  function samples(path, precision) {
    var n = path.getTotalLength(),
      t = [0],
      i = 0,
      dt = precision;
    while ((i += dt) < n) t.push(i);
    t.push(n);
    return t.map(function(t) {
      var p = path.getPointAtLength(t),
        a = [p.x, p.y];
      a.t = t / n;
      return a;
    });
  }

  const link = mainG
    .append("g")
    .attr("stroke-opacity", 0.6)
    .selectAll("line")
    .data(links)
    .enter()
    .append("line")
    .attr("stroke", (d, index) => {
      return "white";
      //return linkColour(data.links[index].value);
    })
    .attr("stroke-width", d => scales.stroke(d.value))
    .attr("x1", e => e.source.x)
    .attr("y1", e => e.source.y)
    .attr("x2", e => e.target.x)
    .attr("y2", e => e.target.y);
  //.attr("stroke", d => "url(#" + d.id + "-gradient)")
  //  .attr("id",d=>d.id+"-gradient")
  //  .attr("stroke",)
  //  .attr("stroke-width", d => scales.stroke(d.value));
  /*links.forEach(link => {
    defs
      .append("linearGradient")
      .attr("id", d => link.id + "-gradient")
      .append("stop")
      .attr("offset", "25%")
      .attr("stop-color", d => {
        return scales.colour(
          totalLinkWeight[link.source.id]["value"] /
            totalLinkWeight[link.source.id]["count"]
        );
      })
      .append("stop")
      .attr("offset", "75%")
      .attr("stop-color", d => {
        //console.log(d);
        return scales.colour(
          totalLinkWeight[link.target.id]["value"] /
            totalLinkWeight[link.target.id]["count"]
        );
      })
      .attr("gradientUnits", "userSpaceOnUse");
  });*/
  /*const gradient = defs
    .append("linearGradient")
    .attr("id", "linearGradient")
    .selectAll("#linearGradient")
    .data(links)
    .attr("id", d => d.id + "-gradient")
    .append("stop")
    .attr("offset", "0%")
    .attr("stop-color", d => {
      console.log(d);
      console.log("helo");
      return scales.colour(
        totalLinkWeight[d.source.id]["value"] /
          totalLinkWeight[d.source.id]["count"]
      );
    })
    .append("stop")
    .attr("offset", "100%")
    .attr("stop-color", d => {
      //console.log(d);
      return scales.colour(
        totalLinkWeight[d.target.id]["value"] /
          totalLinkWeight[d.target.id]["count"]
      );
    });*/

  const node = mainG
    .append("g")
    .attr("stroke", "#fff")
    .attr("stroke-width", 1.5)
    .selectAll("circle")
    .data(nodes)
    .join("circle")
    //.attr("r", v => scales.radius(v.value))
    .attr("r", radius)
    .attr("id", d => d.id)
    .attr("cx", v => v.x)
    .attr("cy", v => v.y)
    .attr("fill", (d, index) => {
      return scales.colour(
        totalLinkWeight[d.id]["value"] / totalLinkWeight[d.id]["count"]
      );
      //geneColours[data.nodes[index]]
      //  return "aqua";
    })
    .on("click", function clicked(d, i, allData) {
      if (d3.event.defaultPrevented) return; // dragged

      //  setSelectedNodes(d3.select(this).attr("id"));
    })
    .call(drag(simulation));

  /*const linkLabels = mainG
    .append("g")
    .selectAll("text")
    .data(links)
    .enter()
    .append("text")
    .attr("dy", function(d, i) {
      return d.target.x > d.source.x ? 18 : -11;
    })
    /* .append("textPath")
    .attr("xlink:href", function(_, i) {
      return "#path" + i;
    })

    .attr("text-anchor", "middle")
    .attr("startOffset", "50%")
    .text((d, i) => {
      return `${data[i].interaction}`;
    })
    .style("color", "black");*/

  const nodeLables = mainG
    .append("g")
    .selectAll("text")
    .data(nodes)
    .enter()
    .append("text")
    .text(function(d) {
      return d.id;
    })
    //  .attr("dx", d => d.x + 250)
    //    .attr("dy", d => d.y + 200)
    .style("fill", "white");

  node.append("title").text(d => d.id);
  function toCartesian(radial, theta) {
    var x = radial * Math.cos(theta);
    var y = radial * Math.sin(theta);
    return { x: x, y: y };
  }

  simulation.on("tick", () => {
    var xExtent = d3.extent(node.data(), function(d) {
      return d.x + xOffSet;
    });
    var yExtent = d3.extent(node.data(), function(d) {
      return d.y + yOffSet;
    });
    // get scales:
    var xScale = width / 2 / (xExtent[1] - xExtent[0]);
    var yScale = height / 2 / (yExtent[1] - yExtent[0]);

    // get most restrictive scale
    var minScale = Math.min(xScale, yScale);

    if (minScale < 1 && false) {
      var transform = d3.zoomIdentity.scale(minScale);
      svg.call(zoom.transform, transform);
    }
    link
      .attr("x1", function(d) {
        return Math.max(radius * 2, Math.min(width - radius * 2, d.source.x));
        //  return d.source.x;
      })
      .attr("y1", function(d) {
        return Math.max(radius * 2, Math.min(height - radius * 2, d.source.y));
        return d.source.y;
      })
      .attr("x2", function(d) {
        return Math.max(radius * 2, Math.min(width - radius * 2, d.target.x));
        //  return d.target.x;
      })
      .attr("y2", function(d) {
        return Math.max(radius * 2, Math.min(height - radius * 2, d.target.y));
        //    return d.target.y;
      });

    node
      .attr("cx", (d, i) => {
        /*  d.radial = radial;
        d.theta = scales.theta(i);

        const coords = toCartesian(d.radial, d.theta);
        d.x = coords.x;*/

        return Math.max(radius * 2, Math.min(width - radius * 2, d.x));

        //return d.x;
      })
      .attr("cy", (d, i) => {
        /*  d.radial = radial;
        d.theta = scales.theta(i);

        const coords = toCartesian(d.radial, d.theta);

        d.y = coords.y;*/
        return Math.max(radius * 2, Math.min(height - radius * 2, d.y));

        //return d.y;
      });
    nodeLables
      .attr("dx", d =>
        //d.x + 10
        Math.max(radius * 2, Math.min(width - radius * 2, d.x + 15))
      )
      .attr("dy", d =>
        //  d.y + 10
        Math.max(radius * 2, Math.min(height - radius * 2, d.y + 5))
      );
    /*   */

    /*linkLabels
      .attr("dy", function(d, i) {
        return d.target.x > d.source.x ? 18 : -15;
      })
      .attr("x", function(d) {
        return (d.target.x + d.source.x) / 2 + xOffSet;
      })
      .attr("y", function(d) {
        return (d.target.y + d.source.y) / 2 + yOffSet;
      });*/
  });
};

const drawLabels = (context, highlightedRow) => {
  context.font = PROP_AXIS_FONT;
  context.fillStyle = "black";
  context.textAlign = "center";
  context.lineWidth = 1;
  context.textBaseline = "bottom";
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

export default Network;
