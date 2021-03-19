import React, { useRef, useCallback } from "react";
import infoText from "./InfoText.js";
import * as d3 from "d3";
import d3Tip from "d3-tip";
const tooltip = d3Tip()
  .style("font-size", "12px")
  .style("width", 150)
  .style("box-shadow", "1px 1px 4px rgba(0,0,0,0.5)")
  .style("border-radius", "none")
  .attr("class", "d3-tip s")
  .attr("data-placement", "bottom")
  .html(function(d) {
    return "<p>" + d + "</p>";
  })
  .offset([-12, 0]);

const Info = ({ name, direction }) => {
  const [ref] = useHookWithRefCallback();

  function useHookWithRefCallback() {
    const ref = useRef(null);
    const setRef = useCallback(node => {
      if (node) {
        const info = d3.select("#" + name + "-info");
        info.call(tooltip);
        const height =
          (infoText[name]["text"].match(/<br>/g) || []).length * 10.66 - 50;

        d3.select("#" + name + "-info")
          .on("mouseover", function(d) {
            tooltip
              .direction(direction)
              .offset(function() {
                if (direction == "n") {
                  return [-height, 0];
                } else if (direction == "s") {
                  return [height, 0];
                } else if (direction == "e") {
                  return [0, height];
                } else if (direction == "w") {
                  return [0, -height];
                }
              })
              .attr("class", "d3-tip " + direction)
              .attr("id", name + "-tip")
              .show(d, info.node())
              .html(function(d) {
                return "<p>" + infoText[name]["text"] + "</p>";
              });
          })
          .on("mouseout", tooltip.hide);
      }
      ref.current = node;
    }, []);

    return [setRef];
  }

  return (
    <svg
      ref={ref}
      xmlns="http://www.w3.org/2000/svg"
      width="16"
      height="16"
      fill="currentColor"
      class="bi bi-info-circle"
      viewBox="0 0 16 16"
      id={name + "-info"}
    >
      <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z" />
      <path d="M8.93 6.588l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM9 4.5a1 1 0 1 1-2 0 1 1 0 0 1 2 0z" />
    </svg>
  );
};

export default Info;
