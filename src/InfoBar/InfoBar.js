import React from "react";
import Grid from "@material-ui/core/Grid";
import Tooltip from "@material-ui/core/Tooltip";
import SvgIcon from "@material-ui/core/SvgIcon";

import { makeStyles } from "@material-ui/core/styles";

const useStyles = makeStyles({
  root: {
    textAlign: "right",
    marginBottom: 10,
  },
  arrow: {
    padding: 15,
    color: "#8a8484",
    backgroundColor: "#F7F8FB",
    borderRadius: 5,
    "&:after": {
      content: "''",
      position: "absolute",
      left: "10%",
      marginLeft: -50,
      borderTop: "solid 50px #F7F8FB",
      borderLeft: "solid 50px transparent",
      borderRight: "solid 50px transparent",
    },
  },
});
const InfoBar = ({ title, infoText }) => {
  const classes = useStyles();
  return (
    <Grid item className={classes.root}>
      <div className={classes.arrow}>
        {title + "    "}

        <InfoIcon infoText={infoText} />
      </div>
    </Grid>
  );
};

const InfoIcon = ({ infoText }) => {
  return (
    <Tooltip title={infoText} arrow>
      <SvgIcon viewBox="0 0 16 16" style={{ fontSize: 16 }}>
        <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14zm0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16z" />
        <path d="M8.93 6.588l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM9 4.5a1 1 0 1 1-2 0 1 1 0 0 1 2 0z" />
      </SvgIcon>
    </Tooltip>
  );
};

export default InfoBar;
