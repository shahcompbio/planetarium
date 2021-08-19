import React, { useEffect, useState } from "react";
import PropTypes from "prop-types";

import * as d3 from "d3";
import { useD3 } from "../utils/useD3";

const SQUARE_LENGTH = 10;
const SQUARE_SPACING = 8;
const STEP = SQUARE_LENGTH + SQUARE_SPACING;

const PADDING = 20;

const Vertical = ({
  ticks,
  colorScale,
  width = 150,
  height = 400,
  title = null,
  disable = false,
  highlight = null,
  onClick = (value) => {},
  onHover = (value) => {},
}) => {
  const isNumerical = typeof ticks === "number";
  const labels = isNumerical ? colorScale.ticks(Math.max(2, ticks)) : ticks;
  const stepWidth = isNumerical ? labels[1] - labels[0] : 0;

  const data = labels.map((datum) =>
    typeof datum === "string" || typeof datum === "number"
      ? { label: datum, value: datum }
      : datum
  );
  const [selected, setSelected] = useState(null);

  useEffect(() => {
    if (highlight !== selected) {
      setSelected(highlight);
    }
  }, [highlight, selected]);

  // useEffect with new param to overwrite selected when needed
  const svgRef = useD3(
    (svg) => {
      const legendHeight = PADDING * 2 + data.length * STEP;
      var offset = 0;
      if (title !== undefined) {
        const titleText = svg
          .append("text")
          .attr("alignment-baseline", "center")
          .attr("dominant-baseline", "center")
          .attr("text-anchor", "middle")
          .attr("font-family", "Helvetica")
          .attr("font-weight", "500")
          .attr("font-size", "12px")
          .attr("fill", "#000000")
          .attr("x", width / 2)
          .attr("y", PADDING)
          .text(title);

        const titleHeight = PADDING + titleText.node().getBBox().height;
        svg.attr("height", legendHeight + titleHeight);
        offset += titleHeight;
      }

      const subsets = svg
        .selectAll("g")
        .attr("pointer-events", "none")
        .data(data)
        .enter()
        .append("g")
        .attr("cursor", disable ? null : "pointer");

      subsets
        .append("rect")
        .attr("pointer-events", "none")
        .attr("width", SQUARE_LENGTH)
        .attr("height", SQUARE_LENGTH)
        .attr("x", 5)
        .attr("y", (d, i) => i * STEP + offset)
        .attr("fill", (d) => colorScale(d["value"]));

      subsets
        .append("text")
        .attr("pointer-events", "none")
        .attr("alignment-baseline", "hanging")
        .attr("dominant-baseline", "hanging")
        .attr("text-align", "left")
        .attr("font-family", "Helvetica")
        .attr("font-weight", "500")
        .attr("font-size", "12px")
        .attr("fill", "#000000")
        .attr("x", SQUARE_LENGTH + 10)
        .attr("y", (d, i) => i * STEP + offset)
        .text((d) => d["label"]);

      const findLabelValue = (mouseY) => {
        const index = Math.round((mouseY - offset) / STEP);

        if (0 <= index && index < data.length) {
          return isNumerical
            ? [data[index]["value"], data[index]["value"] + stepWidth]
            : data[index]["value"];
        }

        return null;
      };

      const mousemove = (d, i, e) => {
        const mouseY = d3.mouse(e[0])[1];

        const label = findLabelValue(mouseY);

        onHover(label);
      };

      const mouseout = () => {
        onHover(null);
      };

      const click = (d, i, e) => {
        const mouseY = d3.mouse(e[0])[1];

        const label = findLabelValue(mouseY);
        const value = isNumerical ? label[0] : label;
        const selectedValue = value === selected ? null : label;
        setSelected(value);
        onClick(selectedValue);
      };

      svg
        .on("mousemove", (d, i, e) => {
          if (disable) {
            return;
          }
          mousemove(d, i, e);
        })
        .on("mouseout", () => {
          if (disable) {
            return;
          }
          mouseout();
        })
        .on("click", (d, i, e) => {
          if (disable) {
            return;
          }
          click(d, i, e);
        });
    },
    width,
    height,
    [title, colorScale, ticks]
  );
  return <svg ref={svgRef} />;
};

Vertical.propTypes = {
  /**
   * width
   */
  width: PropTypes.number.isRequired,
  /**
   * height
   */
  height: PropTypes.number.isRequired,
  /**
   * title for legend
   */
  title: PropTypes.string,
  /**
   * Number of ticks to show
   */
  ticks: PropTypes.oneOfType([
    PropTypes.number,
    PropTypes.arrayOf(PropTypes.object),
    PropTypes.arrayOf(PropTypes.string),
  ]),
  /**
   * Color scale for legend items
   */
  colorScale: PropTypes.func.isRequired,
  /**
   * Whether interactions should be disabled
   */
  disable: PropTypes.bool,
  /**
   * Label to highlight
   */
  highlight: PropTypes.oneOfType([
    PropTypes.number,
    PropTypes.object,
    PropTypes.string,
  ]),
  /**
   * callback for on hover on square
   */
  onHover: PropTypes.func,
  /**
   * callback for on click on square
   */
  onClick: PropTypes.func,
};

export default Vertical;
