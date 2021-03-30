import React from "react";
import Grid from "@material-ui/core/Grid";
import Profile from "./Profile";
import ProfileBackground from "./ProfileBackground";
import GenomeYAxis from "./GenomeYAxis";
import ChromosomeAxis from "./ChromsomeAxis";
import { Y_AXIS_WIDTH, X_AXIS_HEIGHT, TOP_PADDING } from "./utils";

const GenomeProfile = ({
  bins,
  segs,
  bpTotal,
  chromosomes,
  width,
  height,
  maxState,
}) => {
  const profileWidth = width - Y_AXIS_WIDTH;
  const profileHeight = height - X_AXIS_HEIGHT - TOP_PADDING;

  return (
    <Grid container direction="column">
      <Grid
        item
        container
        direction="row"
        id="profile-wrapper"
        style={{ position: "relative" }}
        width={width}
        height={height}
      >
        <Grid item>
          <GenomeYAxis
            height={profileHeight}
            width={Y_AXIS_WIDTH}
            maxState={maxState}
          />
        </Grid>
        <Grid item style={{ "padding-top": TOP_PADDING }}>
          <ProfileBackground
            chromosomes={chromosomes}
            bpTotal={bpTotal}
            width={profileWidth}
            height={profileHeight}
            maxState={maxState}
          />
          <Profile
            bins={bins}
            segs={segs}
            bpTotal={bpTotal}
            width={profileWidth}
            height={profileHeight}
            maxState={maxState}
          />
        </Grid>
      </Grid>
      <Grid item style = {{marginTop:-22}}>
        <ChromosomeAxis
          chromosomes={chromosomes}
          bpTotal={bpTotal}
          width={profileWidth}
          height={X_AXIS_HEIGHT}
        />
      </Grid>
    </Grid>
  );
};

export default GenomeProfile;