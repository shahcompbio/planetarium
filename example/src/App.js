import React, { useState } from "react";
import Button from "@material-ui/core/Button";
import Drawer from "@material-ui/core/Drawer";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";

import StaticFigures from "./StaticFigures/StaticFigures";
import staticFigureFetchData from "./StaticFigures/data/api";

import StaticPeds from "./StaticFigures/StaticPeds";
import pedsFetchData from "./StaticFigures/data/apiPeds";

import StaticFishTail from "./StaticFishTail/StaticFishTail";
import staticFishData from "./StaticFishTail/data/api";

import PaperFormatter from "./PaperFormatter/PaperFormatter";

import UMAP from "./UMAP/UMAP";
import fetchFileData from "./UMAP/data/api";

import Newick from "./NewickTree/Newick";
import fetchTree from "./NewickTree/data/api";

import Test from "./Test/Test";
import * as d3 from "d3";

const componentList = [
  "UMAP",
  "Static Heatmap",
  "Static FishTail",
  "Test",
  "Newick",
  "Peds",
];
const getAppComponent = (selection, data) => {
  console.log(data);
  switch (selection) {
    case "UMAP":
      return (
        <UMAP dashboardID={"Test"} api={"http://localhost:9200"} data={data} />
      );
      break;
    case "Static FishTail":
      return (
        <StaticFishTail
          dashboardID={"Test"}
          api={"http://localhost:9200"}
          data={data["data"]}
          cloneColor={data["clones"]}
        />
      );
      break;
    case "Static Heatmap":
      return (
        <StaticFigures
          dashboardID={"Test"}
          api={"http://localhost:9200"}
          data={data}
        />
      );
      break;
    case "Test":
      return (
        <Test
          dashboardID={"Layers"}
          api={"http://localhost:9200"}
          data={data["metadata"]}
          filters={data["filters"]}
        />
      );
    case "Newick":
      return (
        <Newick
          dashboardID={"Layers"}
          api={"http://localhost:9200"}
          data={data["data"]}
        />
      );
      break;
    case "Peds":
      return (
        <StaticPeds
          dashboardID={"Layers"}
          api={"http://localhost:9200"}
          data={data}
        />
      );
      break;
  }
};

function getAppData(selection) {
  switch (selection) {
    case "UMAP":
      return fetchFileData();
      break;
    case "Static FishTail":
      return staticFishData();
      break;
    case "Static Heatmap":
      return staticFigureFetchData();
      break;
    case "Test":
      return fetchFileData();
    case "Newick":
      return fetchTree();
      break;
    case "Peds":
      return pedsFetchData();
      break;
  }
}
const ComponentList = ({ setSelectedDashboard }) => {
  return (
    <List>
      {componentList.map((component) => (
        <ListItem disablePadding key={"listitem" + component}>
          <ListItemButton
            key={"listbuttom" + component}
            onClick={() => {
              localStorage.setItem("planetariumDashboard", component);
              setSelectedDashboard(component);
            }}
          >
            <ListItemText key={"listtext" + component} primary={component} />
          </ListItemButton>
        </ListItem>
      ))}
    </List>
  );
};
const SideDrawer = ({ setSelectedDashboard }) => {
  const [drawerOpen, toggleDrawer] = useState(false);
  return (
    <React.Fragment key={"bottom"}>
      <Button onClick={() => toggleDrawer(true)}>Select Dashboard</Button>
      <Drawer
        anchor={"bottom"}
        open={drawerOpen}
        onClose={() => toggleDrawer(false)}
      >
        {
          <ComponentList
            setSelectedDashboard={(selection) =>
              setSelectedDashboard(selection)
            }
          />
        }
      </Drawer>
    </React.Fragment>
  );
};
const DevAppWrapper = () => {
  //const [selectedDashboard, setSelectedDashboard] = useState("Static FishTail");
  const [selectedDashboard, setSelectedDashboard] = useState(
    localStorage.getItem("planetariumDashboard") === null
      ? "Static FishTail"
      : localStorage.getItem("planetariumDashboard")
  );
  return (
    <DevApp
      selectedDashboard={selectedDashboard}
      setSelectedDashboard={setSelectedDashboard}
    />
  );
};
const DevApp = ({ selectedDashboard, setSelectedDashboard }) => {
  const data = getAppData(selectedDashboard);
  d3.select("iframe").style("z-index", -10);

  return Object.keys(data).length === 0 ? null : (
    <div>
      <SideDrawer
        setSelectedDashboard={(selection) => setSelectedDashboard(selection)}
      />
      {getAppComponent(selectedDashboard, data)}
    </div>
  );
};

const ProdApp = () => (
  <UMAP
    data={window.isablData}
    dashboardID={window.dashboardID}
    api={window.apiURL}
  />
);

export default process.env.NODE_ENV === "development" ? DevAppWrapper : ProdApp;
