/* 
Tooltip is primarily for svg/canvas related components. This requires some setup in the parent component to use properly:

1. svg and Tooltip must be wrapped in parent div with relative positioning:
<div style={{position: "relative"}}>
  <svg/>
  <Tooltip/>
</div>

2. Parent component should use state to keep track of where the tooltip is shown (like highlightedNode or highlightedBin) = default state is null

*/

import React from "react";
import { Tooltip } from "@material-ui/core";

const TooltipP = ({ getText, getX, getY, data }) => {
  return (
    <Tooltip
      title={data ? getText(data) : ""}
      open={data !== null}
      arrow
      placement="top"
    >
      <div
        style={{
          position: "absolute",
          pointerEvents: "none",
          left: data ? getX(data) : null,
          top: data ? getY(data) : null,
        }}
      />
    </Tooltip>
  );
};

export default TooltipP;
