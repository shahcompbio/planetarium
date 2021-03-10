import React, { useState } from "react";
import LayoutWrapper from "../Layout/LayoutWrapper.js";
import "./App.css";
import dashboardReducer, {
  initialState,
} from "../PlotState/dashboardReducer.js";
import { DashboardProvider } from "../PlotState/dashboardState.js";

const App = ({ data }) => {
  const [selectedSubtype, setSelectedSubtype] = useState(null);
  const [selectedClonotype, setSelectedClonotype] = useState(null);

  return (
    <div className="App">
      <DashboardProvider
        initialState={{
          ...initialState,
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
