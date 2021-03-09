import React, { useState, useEffect } from "react";
import logo from "./logo.svg";
import LayoutWrapper from "./Layout/LayoutWrapper.js";
import "./App.css";
import _ from "lodash";
import * as d3 from "d3";
import dataSource from "./VDJ_ANALYSIS_8425_metadata.tsv";
import dashboardReducer, {
  initialState
} from "./PlotState/dashboardReducer.js";
import { DashboardProvider } from "./PlotState/dashboardState.js";

const App = ({}) => {
  const [selectedSubtype, setSelectedSubtype] = useState(null);
  const [selectedClonotype, setSelectedClonotype] = useState(null);
  const [data, setData] = useState([]);

  useEffect(() => {
    d3.tsv(dataSource).then(data => {
      setData(data);
    });
  }, []);

  return (
    <div className="App">
      <DashboardProvider
        initialState={{
          ...initialState
        }}
        reducer={dashboardReducer}
      >
        <div className="App">
          <LayoutWrapper data={data} />
        </div>
      </DashboardProvider>
    </div>
  );
};

export default App;
