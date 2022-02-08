import React from "react";

import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";
import Paper from "@mui/material/Paper";
import ListItemText from "@mui/material/ListItemText";
import MenuList from "@mui/material/MenuList";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import BiotechIcon from "@mui/icons-material/Biotech";
import BubbleChartIcon from "@mui/icons-material/BubbleChart";
import LayersClearIcon from "@mui/icons-material/LayersClear";
import makeStyles from "@mui/styles/makeStyles";

const useStyles = makeStyles({
  rootWithMarginTop: {
    backgroundColor: "none",
    paddingLeft: 17,
    paddingRight: 6,
    paddingBottom: 0,
    minWidth: 275,
    marginTop: 15,
    marginLeft: 15,
    marginBottom: 0,
  },
});

const MetaData = ({ sample, hasSelection, totalCount, setHighlight }) => {
  const classes = useStyles();
  return (
    <Grid
      container
      direction="column"
      justifyContent="flex-start"
      alignItems="stretch"
    >
      <Grid item>
        <Header classes={classes} sample={sample} totalCount={totalCount} />

        <Paper
          sx={{
            width: "100%",
            maxWidth: "252px",
            margin: 1,
            marginLeft: 8,
            marginTop: 0,
            padding: 2,
            paddingRight: 0,
          }}
        >
          <MenuList style={{ padding: 0 }}>
            <MenuItem
              key={"clear-menuItem"}
              disabled={!hasSelection}
              onClick={() => {
                setHighlight();
              }}
              style={{ height: "100%", width: "100%" }}
            >
              <ListItemIcon key={"clear-listItemIcon"}>
                <LayersClearIcon fontSize="small" />
              </ListItemIcon>
              <ListItemText key={"clear-listItemText"}>Clear</ListItemText>
              <div style={{ paddingRight: 21 }}>
                <Typography variant="body2" color="text.secondary">
                  ⌘X
                </Typography>
              </div>
            </MenuItem>
          </MenuList>
        </Paper>
      </Grid>
    </Grid>
  );
};

const Header = ({ classes, sample, totalCount }) => (
  <div className={classes.rootWithMarginTop}>
    <MenuList key={"header-menuList"}>
      <Paper
        sx={{ width: "100%", maxWidth: "100%", padding: 2, marginBottom: 2 }}
      >
        <MenuItem key={"sample-menuItem"}>
          <ListItemIcon key={"sample-listItemIcon"}>
            <BiotechIcon fontSize="small" key={"sample-bioTechIcon"} />
          </ListItemIcon>
          <ListItemText key={"sample-listItemText"}>{sample}</ListItemText>
        </MenuItem>
      </Paper>
      <Paper sx={{ width: "100%", maxWidth: "100%", padding: 2 }}>
        <MenuItem key={"cellCount-menuItem"}>
          <ListItemIcon key={"cellCount-listItemIcon"}>
            <BubbleChartIcon
              fontSize="small"
              key={"cellCount-bubbleChartIcon"}
            />
          </ListItemIcon>
          <ListItemText key={"cellCount-listItemText"}>
            {totalCount} cells
          </ListItemText>
        </MenuItem>
      </Paper>
    </MenuList>
  </div>
);

export default MetaData;
