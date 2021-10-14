import React, { useState } from "react";
import PropTypes from "prop-types";
import makeStyles from '@mui/styles/makeStyles';

import Tooltip from "@mui/material/Tooltip";
import SvgIcon from "@mui/material/SvgIcon";
import Collapse from "@mui/material/Collapse";
import MUISearchIcon from "@mui/icons-material/Search";

const useStyles = makeStyles({
  svgIcon: { fontSize: 16, marginLeft: 10 },
});

export const SearchIcon = ({ children }) => {
  const classes = useStyles();
  const [isOpen, setIsOpen] = useState(false);

  return (
    <span style={{ display: "flex", float: "left" }}>
      <Collapse direction={"left"} in={isOpen}>
        {children}
      </Collapse>
      <Tooltip title={"Search"} arrow>
        <SvgIcon
          viewBox="0 0 20 20"
          className={classes.svgIcon}
          onClick={() => setIsOpen(!isOpen)}
        >
          <MUISearchIcon fontSize="medium" />
        </SvgIcon>
      </Tooltip>
    </span>
  );
};

SearchIcon.propTypes = {
  /**
   * components to show/hide on click
   */
  children: PropTypes.oneOf([
    PropTypes.elementType,
    PropTypes.arrayOf(PropTypes.elementType),
  ]),
};

export const DownloadIcon = ({ download }) => {
  const classes = useStyles();

  return (
    <Tooltip title={"Download"} arrow>
      <SvgIcon
        viewBox="0 0 20 20"
        className={classes.svgIcon}
        onClick={download}
      >
        <path d="M 6.089844 8.722656 L 11.378906 3.433594 L 12.847656 4.902344 L 7.558594 10.191406 Z M 6.089844 8.722656 " />
        <path d="M 6.375 0.429688 L 8.757812 0.429688 L 8.757812 7.953125 L 6.375 7.953125 Z M 6.375 0.429688 " />
        <path d="M 7.566406 10.195312 L 2.277344 4.910156 L 3.742188 3.441406 L 9.035156 8.730469 Z M 7.566406 10.195312 " />
        <path d="M 1.125 12.445312 L 14 12.445312 L 14 14.570312 L 1.125 14.570312 Z M 1.125 12.445312 " />
        <path d="M 0.6875 8.570312 L 2.8125 8.570312 L 2.8125 14.570312 L 0.6875 14.570312 Z M 0.6875 8.570312 " />
        <path d="M 12.1875 8.570312 L 14.3125 8.570312 L 14.3125 14.570312 L 12.1875 14.570312 Z M 12.1875 8.570312 " />
      </SvgIcon>
    </Tooltip>
  );
};

DownloadIcon.propTypes = {
  /**
   * download function
   */
  download: PropTypes.func.isRequired,
};
