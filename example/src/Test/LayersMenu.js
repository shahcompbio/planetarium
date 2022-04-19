import React from "react";
import Card from "@mui/material/Card";
import CardActions from "@mui/material/CardActions";
import CardContent from "@mui/material/CardContent";
import Typography from "@mui/material/Typography";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import ListSubheader from "@mui/material/ListSubheader";
import TreeItem, { treeItemClasses } from "@mui/lab/TreeItem";
import TreeView from "@mui/lab/TreeView";
import SvgIcon from "@mui/material/SvgIcon";
import { alpha, styled } from "@mui/material/styles";
import { removeFilters, setMainStateWithNewFilters } from "./Filters.js";

function CloseSquare(props) {
  return (
    <SvgIcon
      className="close"
      fontSize="inherit"
      style={{ width: 14, height: 14 }}
      {...props}
    >
      {/* tslint:disable-next-line: max-line-length */}
      <path d="M17.485 17.512q-.281.281-.682.281t-.696-.268l-4.12-4.147-4.12 4.147q-.294.268-.696.268t-.682-.281-.281-.682.294-.669l4.12-4.147-4.12-4.147q-.294-.268-.294-.669t.281-.682.682-.281.696 .268l4.12 4.147 4.12-4.147q.294-.268.696-.268t.682.281 .281.669-.294.682l-4.12 4.147 4.12 4.147q.294.268 .294.669t-.281.682zM22.047 22.074v0 0-20.147 0h-20.12v0 20.147 0h20.12zM22.047 24h-20.12q-.803 0-1.365-.562t-.562-1.365v-20.147q0-.776.562-1.351t1.365-.575h20.147q.776 0 1.351.575t.575 1.351v20.147q0 .803-.575 1.365t-1.378.562v0z" />
    </SvgIcon>
  );
}
function MinusSquare(props) {
  return (
    <SvgIcon fontSize="inherit" style={{ width: 14, height: 14 }} {...props}>
      {/* tslint:disable-next-line: max-line-length */}
      <path d="M22.047 22.074v0 0-20.147 0h-20.12v0 20.147 0h20.12zM22.047 24h-20.12q-.803 0-1.365-.562t-.562-1.365v-20.147q0-.776.562-1.351t1.365-.575h20.147q.776 0 1.351.575t.575 1.351v20.147q0 .803-.575 1.365t-1.378.562v0zM17.873 11.023h-11.826q-.375 0-.669.281t-.294.682v0q0 .401.294 .682t.669.281h11.826q.375 0 .669-.281t.294-.682v0q0-.401-.294-.682t-.669-.281z" />
    </SvgIcon>
  );
}

function PlusSquare(props) {
  return (
    <SvgIcon fontSize="inherit" style={{ width: 14, height: 14 }} {...props}>
      {/* tslint:disable-next-line: max-line-length */}
      <path d="M22.047 22.074v0 0-20.147 0h-20.12v0 20.147 0h20.12zM22.047 24h-20.12q-.803 0-1.365-.562t-.562-1.365v-20.147q0-.776.562-1.351t1.365-.575h20.147q.776 0 1.351.575t.575 1.351v20.147q0 .803-.575 1.365t-1.378.562v0zM17.873 12.977h-4.923v4.896q0 .401-.281.682t-.682.281v0q-.375 0-.669-.281t-.294-.682v-4.896h-4.923q-.401 0-.682-.294t-.281-.669v0q0-.401.281-.682t.682-.281h4.923v-4.896q0-.401.294-.682t.669-.281v0q.401 0 .682.281t.281.682v4.896h4.923q.401 0 .682.281t.281.682v0q0 .375-.281.669t-.682.294z" />
    </SvgIcon>
  );
}
const StyledTreeItem = styled((props) => <TreeItem {...props} />)(
  ({ theme }) => ({
    [`& .${treeItemClasses.iconContainer}`]: {
      "& .close": {
        opacity: 0.3,
      },
    },
    [`& .${treeItemClasses.group}`]: {
      marginLeft: 15,
      paddingLeft: 18,
      borderLeft: `1px dashed ${alpha(theme.palette.text.primary, 0.4)}`,
    },
  })
);

const LayersMenu = ({
  layerName,
  thisLayer,
  layerFilters,
  setLayerFilters,
}) => {
  const currentFilters = layerFilters[thisLayer]["filters"];
  const currentLasso = layerFilters[thisLayer]["lasso"];
  return (
    <Card sx={{ minWidth: 275, marginBottom: 5 }} key={"card-" + thisLayer}>
      <CardContent key={"cardContent-" + thisLayer}>
        <Typography variant="h5" component="div" key={"title-" + thisLayer}>
          {layerName}
        </Typography>
      </CardContent>
      <CardActions key={"cardActions-" + thisLayer}>
        {currentFilters && (
          <List
            sx={{
              width: "100%",
              maxWidth: 360,
              bgcolor: "background.paper",
              position: "relative",
              overflow: "auto",
              maxHeight: 300,

              "& ul": { padding: 0 },
            }}
            subheader={<li />}
          >
            <TreeView
              aria-label={"layers-" + thisLayer}
              sx={{
                height: 240,
                flexGrow: 1,
                maxWidth: 400,
                overflowY: "auto",
              }}
              defaultExpanded={["filters-" + thisLayer, "lasso-" + thisLayer]}
              defaultCollapseIcon={<MinusSquare />}
              defaultExpandIcon={<PlusSquare />}
              defaultEndIcon={
                <CloseSquare
                  onClick={(event, d) => {
                    var newFilters = layerFilters;
                    //newFilters[thisLayer]["filters"] =
                  }}
                />
              }
            >
              {currentFilters && (
                <StyledTreeItem
                  key={"filters-" + thisLayer}
                  nodeId={"filters-" + thisLayer}
                  label={"Applied Filters"}
                >
                  {Object.keys(currentFilters).map((filter, index) => (
                    <StyledTreeItem
                      key={filter + index + thisLayer}
                      nodeId={filter + index + thisLayer}
                      label={filter}
                    >
                      {currentFilters[filter].map((item, index) => (
                        <StyledTreeItem
                          key={item + index + thisLayer}
                          nodeId={item + index + thisLayer}
                          label={item}
                          onClick={(event) => {
                            const newFilters = removeFilters(
                              currentFilters,
                              filter,
                              item
                            );
                            setMainStateWithNewFilters(
                              layerFilters,
                              thisLayer,
                              newFilters,
                              setLayerFilters
                            );
                          }}
                        />
                      ))}
                    </StyledTreeItem>
                  ))}
                </StyledTreeItem>
              )}
              {currentLasso && (
                <StyledTreeItem
                  key={"lasso-" + thisLayer}
                  nodeId={"lasso-" + thisLayer}
                  label={"Lasso Selection"}
                >
                  <StyledTreeItem
                    key={"lasso-item-" + thisLayer}
                    nodeId={"lasso-item-" + thisLayer}
                    label={currentLasso["selectionLength"] + " selected"}
                    onClick={(event) => {}}
                  />
                </StyledTreeItem>
              )}
            </TreeView>
          </List>
        )}
      </CardActions>
    </Card>
  );
};

export default LayersMenu;
