import React, { useState } from "react";
import PropTypes from "prop-types";
import { makeStyles } from "@material-ui/core/styles";

import Tooltip from "@material-ui/core/Tooltip";
import SvgIcon from "@material-ui/core/SvgIcon";
import Collapse from "@material-ui/core/Collapse";
import MUISearchIcon from "@material-ui/icons/Search";
import Avatar from "@material-ui/core/Avatar";

const useStyles = makeStyles({
  svgIcon: { fontSize: 30 },
});

export const SearchIcon = ({ children }) => {
  const classes = useStyles();
  const [isOpen, setIsOpen] = useState(false);

  return (
    <span style={{ display: "flex", float: "left", paddingRight: 10 }}>
      <Tooltip title={"Search"} arrow>
        <Avatar variant="rounded" style={{ width: isOpen ? 220 : 35 }}>
          <Collapse direction={"left"} in={isOpen}>
            {children}
          </Collapse>
          <SvgIcon
            viewBox="2 0 40 40"
            className={classes.svgIcon}
            onClick={() => setIsOpen(!isOpen)}
          >
            <MUISearchIcon sx={{ fontSize: 40 }} />
          </SvgIcon>
        </Avatar>
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
    <span style={{ paddingRight: 10 }}>
      <Avatar variant="rounded">
        <Tooltip title={"Download"} arrow>
          <SvgIcon
            viewBox="-2 -2 20 20"
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
      </Avatar>
    </span>
  );
};

DownloadIcon.propTypes = {
  /**
   * download function
   */
  download: PropTypes.func.isRequired,
};
