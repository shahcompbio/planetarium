import React from "react";
import InfoBar from "./InfoBar";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

const MARGIN = 10;
const PADDING = 10;

const Layout = ({ title, infoText, children }) => {
  return (
    <Paper
      style={{
        margin: MARGIN,
      }}
    >
      <Grid
        container
        direction="column"
        justify="flex-start"
        alignItems="stretch"
      >
        <InfoBar title={title} infoText={infoText} />
        <Grid item style={{ padding: PADDING }}>
          {children}
        </Grid>
      </Grid>
    </Paper>
  );
};

export default Layout;
