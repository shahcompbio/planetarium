import React, { useState } from "react";
import Grid from "@mui/material/Grid";

import Heatmap, { HEATMAP_SPACING } from "./Heatmap";
import ChromosomeAxis, { HEIGHT as AXIS_HEIGHT } from "./ChromosomeAxis";
import Profile, { Y_AXIS_WIDTH } from "./Profile";

const PROFILE_HEIGHT = 250;
const MIN_HEATMAP_HEIGHT = 300;

const CellScape = ({ data, width, height, chromosomes }) => {
  const [highlighted, setHighlighted] = useState(null);
  const heatmapWidth = width - HEATMAP_SPACING;

  // will probably need a useEffect hook and API prop to query for bin data, as having it all in here will be too much
  return (
    <Grid container direction="column">
      <Grid item style={{ paddingLeft: Y_AXIS_WIDTH }}>
        <Heatmap
          data={data}
          width={width - Y_AXIS_WIDTH}
          height={Math.max(
            height - AXIS_HEIGHT - PROFILE_HEIGHT,
            MIN_HEATMAP_HEIGHT
          )}
          chromosomes={chromosomes}
          onChange={setHighlighted}
        />
      </Grid>
      <Grid item style={{ paddingLeft: Y_AXIS_WIDTH }}>
        <ChromosomeAxis
          chromosomes={chromosomes}
          width={heatmapWidth - Y_AXIS_WIDTH}
        />
      </Grid>
      <Grid item>
        <Profile
          bins={[]}
          segs={highlighted === null ? [] : highlighted["segs"]}
          chromosomes={chromosomes}
          width={heatmapWidth}
          height={PROFILE_HEIGHT}
          maxState={12}
        />
      </Grid>
    </Grid>
  );
};

export default CellScape;
