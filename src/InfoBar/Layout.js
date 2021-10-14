import React from "react";
import PropTypes from "prop-types";

import InfoBar from "./InfoBar";
import Grid from "@mui/material/Grid";
import Paper from "@mui/material/Paper";

const MARGIN = 10;
const PADDING = 10;

const Layout = ({ title, infoText, children, addIcon = null }) => {
  return (
    <Paper
      sx={{
        m: 1,
      }}
    >
      <Grid
        container
        direction="column"
        justifyContent="flex-start"
        alignItems="stretch"
      >
        <InfoBar title={title} infoText={infoText} addIcon={addIcon} />
        <Grid item sx={{ p: 1 }}>
          {children}
        </Grid>
      </Grid>
    </Paper>
  );
};

Layout.propTypes = {
  /**
   * text on information bar
   */
  title: PropTypes.string,
  /**
   * text on tooltip hover
   */
  infoText: PropTypes.string,
  /**
   * additional components to be added to info bar
   */
  addIcon: PropTypes.oneOf([
    PropTypes.elementType,
    PropTypes.arrayOf(PropTypes.elementType),
  ]),
  /**
   * items to put in layout grid
   */
  children: PropTypes.oneOf([
    PropTypes.elementType,
    PropTypes.arrayOf(PropTypes.elementType),
  ]),
};

export default Layout;
