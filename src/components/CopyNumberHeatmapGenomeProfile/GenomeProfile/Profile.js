import React from "react";
import {useD3} from "../../utils/useD3";
import { getGenomeYScale, colorScale } from "./utils";

const Profile = ({ bins, segs, bpTotal, width, height, maxState }) => {
  const ref = useD3((canvas) => {
    const genomeYScale = getGenomeYScale(maxState, height);
    const bpRatio = width / bpTotal;
    let context = canvas.node().getContext("2d");

    bins.forEach((bin) => {
      const x = bin.genomeStart * bpRatio;
      const y = genomeYScale(bin.copy);

      context.fillStyle = colorScale(bin.state);
      context.beginPath();
      context.arc(x, y, 1.5, 0, 2 * Math.PI);
      context.fill();
    });

    segs.forEach((segment) => {
      const x = segment.genomeStart * bpRatio;
      const y = genomeYScale(segment.state);
      const width = segment.length * bpRatio;

      context.fillStyle = "#000000";
      context.fillRect(x, y, width, 2);
    });
  },width,height,[]);

  return (
    <canvas
      id="profile-canvas"
      ref={ref}
      width={width}
      height={height}
      style={{ position: "absolute" }}
    />
  );
};
export default Profile;