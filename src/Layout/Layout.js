import React from "react";

import dashboardReducer, {
  initialState
} from "../PlotState/dashboardReducer.js";
import { DashboardProvider } from "../PlotState/dashboardState.js";

import chartDim from "./LayoutConfig.js";
import componentList from "./ComponentList.js";

const App = ({ chartName, data }) => {
  const Component = componentList[chartName];
  return (
    <DashboardProvider initialState={initialState} reducer={dashboardReducer}>
      <div className="App">
        <Component data={data} chartDim={chartDim} />
      </div>
    </DashboardProvider>
  );
};

export default App;
