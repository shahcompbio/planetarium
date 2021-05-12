import React from "react";
import { useD3 } from "../utils/useD3";
import { useCanvas } from "../utils/useCanvas";
import * as d3 from "d3";
import Grid from "@material-ui/core/Grid";
import { colorScale, BACKGROUND_COLORS } from "./utils";

const Y_AXIS_WIDTH = 25;
const TOP_PADDING = 5;

const GenomeProfile = ({
  bins,
  segs,
  chromosomes,
  width,
  height,
  maxState,
}) => {
  const yScale = d3.scaleLinear().domain([-0.5, maxState]).range([height, 0]);
  const profileWidth = width - Y_AXIS_WIDTH;
  return (
    <Grid container direction="row">
      <Grid item>
        <YAxis width={Y_AXIS_WIDTH} height={height} yScale={yScale} />
      </Grid>
      <Grid item style={{ paddingTop: TOP_PADDING }}>
        <Profile
          bins={bins}
          segs={segs}
          chromosomes={chromosomes}
          width={profileWidth}
          height={height}
          yScale={yScale}
          maxState={maxState}
        />
      </Grid>
    </Grid>
  );
};

const YAxis = ({ width, height, yScale }) => {
  const ref = useD3(
    (svg) => {
      const yAxisTicks = yScale
        .ticks()
        .filter((tick) => Number.isInteger(tick));

      var yAxis = d3
        .axisLeft(yScale)
        .tickValues(yAxisTicks)
        .tickFormat(d3.format("d"));

      svg
        .append("g")
        .attr("class", "genome-y-axis")
        .style("transform", `translate(${width - 2}px, ${TOP_PADDING}px)`)
        .call(yAxis);
    },
    width,
    height,
    []
  );
  return <svg ref={ref} width={width} height={height} />;
};

const Profile = ({
  bins,
  segs,
  chromosomes,
  width,
  height,
  yScale,
  maxState,
}) => {
  const bpTotal = chromosomes.reduce((currSum, chr) => currSum + chr.length, 0);

  const ref = useCanvas(
    (canvas) => {
      const bpRatio = width / bpTotal;
      let context = canvas.getContext("2d");
      context.globalAlpha = 1;
      drawBackground(
        context,
        chromosomes,
        maxState,
        width,
        height,
        yScale,
        bpRatio
      );
      drawProfile(context, bins, segs, yScale, bpRatio);
    },
    width,
    height,
    [bins, segs]
  );

  return <canvas id="profile-canvas" ref={ref} width={width} height={height} />;
};

const drawBackground = (
  context,
  chromosomes,
  maxState,
  width,
  height,
  yScale,
  bpRatio
) => {
  const copyNumbers = Array.from(Array(maxState + 1).keys());

  chromosomes.forEach((chr, index) => {
    context.fillStyle = BACKGROUND_COLORS[index % BACKGROUND_COLORS.length];

    context.fillRect(
      chr.genomeStart * bpRatio,
      0,
      chr.length * bpRatio,
      height
    );
  });

  context.strokeStyle = "#aaaaaa";
  copyNumbers.forEach((cn) => {
    context.beginPath();
    context.setLineDash([1, 2]);
    context.moveTo(0, yScale(cn));
    context.lineTo(width, yScale(cn));
    context.stroke();
  });
};

const drawProfile = (context, bins, segs, yScale, bpRatio) => {
  bins.forEach((bin) => {
    const x = bin.genomeStart * bpRatio;
    const y = yScale(bin.copy);

    context.fillStyle = colorScale(bin.state);
    context.beginPath();
    context.arc(x, y, 1.5, 0, 2 * Math.PI);
    context.fill();
  });

  segs.forEach((segment) => {
    const x = segment.genomeStart * bpRatio;
    const y = yScale(segment.state);
    const width = segment.length * bpRatio;

    context.fillStyle = "#000000";
    context.fillRect(x, y, width, 2);
  });
};

export default GenomeProfile;
