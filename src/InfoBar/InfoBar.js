import React from "react";
import PropTypes from "prop-types";

import Grid from "@material-ui/core/Grid";
import Tooltip from "@material-ui/core/Tooltip";
import SvgIcon from "@material-ui/core/SvgIcon";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
  root: {
    textAlign: "right",
    marginBottom: 10,
    padding: 10,
  },
  svgIcon: { fontSize: 16, marginLeft: 10 },
  arrow: {
    padding: 15,
    color: "#8a8484",
    backgroundColor: "#F7F8FB",
    position: "relative",
    borderRadius: 5,
    "&:after": {
      content: "''",
      position: "absolute",
      left: 50,
      marginLeft: -50,
      borderTop: "solid 50px #F7F8FB",
      borderLeft: "solid 50px transparent",
      borderRight: "solid 50px transparent",
    },
  },
});
const InfoBar = ({ title, infoText, addIcon = null }) => {
  const classes = useStyles();
  return (
    <Grid item className={classes.root}>
      <div className={classes.arrow}>
        <Grid
          container
          direction="row"
          justify="space-between"
          alignItems="center"
        >
          <Grid item xs={3} style={{ textAlign: "left" }}>
            {title}
          </Grid>
          <Grid
            xs={9}
            container
            direction="row"
            justify="flex-end"
            alignItems="center"
          >
            {addIcon}
            <InfoIcon infoText={infoText} classes={classes} />
          </Grid>
        </Grid>
      </div>
    </Grid>
  );
};

const InfoIcon = ({ infoText, classes }) => {
  return (
    <Tooltip title={infoText} arrow>
      <SvgIcon viewBox="0 0 25 25" className={classes.svgIcon}>
        <path d="M 10 0 C 4.476562 0 0 4.476562 0 10 C 0 15.523438 4.476562 20 10 20 C 15.523438 20 20 15.523438 20 10 C 20 4.476562 15.523438 0 10 0 Z M 9.949219 15.59375 C 9.949219 15.59375 9.484375 15.902344 9.167969 15.898438 C 9.140625 15.914062 9.125 15.921875 9.125 15.921875 L 9.125 15.898438 C 9.058594 15.898438 8.992188 15.886719 8.921875 15.871094 L 8.824219 15.847656 C 8.269531 15.707031 7.976562 15.160156 8.113281 14.605469 L 8.886719 11.496094 L 9.234375 10.09375 C 9.558594 8.789062 8.207031 10.371094 7.929688 9.769531 C 7.746094 9.371094 8.980469 8.539062 9.882812 7.910156 C 9.882812 7.910156 10.34375 7.605469 10.660156 7.605469 C 10.691406 7.59375 10.707031 7.585938 10.707031 7.585938 L 10.707031 7.605469 C 10.773438 7.609375 10.839844 7.617188 10.90625 7.632812 L 11.007812 7.660156 C 11.5625 7.796875 11.898438 8.355469 11.761719 8.910156 L 10.992188 12.019531 L 10.640625 13.421875 C 10.320312 14.726562 11.648438 13.140625 11.921875 13.742188 C 12.105469 14.136719 10.851562 14.964844 9.949219 15.59375 Z M 11.90625 5.875 C 11.710938 6.652344 10.925781 7.125 10.152344 6.933594 C 9.375 6.738281 8.902344 5.957031 9.097656 5.179688 C 9.289062 4.40625 10.074219 3.933594 10.847656 4.125 C 11.625 4.316406 12.097656 5.101562 11.90625 5.875 Z M 11.90625 5.875 " />
      </SvgIcon>
    </Tooltip>
  );
};

InfoBar.propTypes = {
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
};

export default InfoBar;
