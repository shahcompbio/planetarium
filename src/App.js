import React, { useState, useEffect } from "react";
import logo from "./logo.svg";
import LayoutWrapper from "./Layout/LayoutWrapper.js";
import "./App.css";
import _ from "lodash";
import * as d3 from "d3";
import metadataSource from "./metadata.tsv";
import probabilitiesSource from "./probabilities.tsv";
import dashboardReducer, {
  initialState
} from "./PlotState/dashboardReducer.js";
import { DashboardProvider } from "./PlotState/dashboardState.js";

const App = ({}) => {
  const [selectedSubtype, setSelectedSubtype] = useState(null);
  const [selectedClonotype, setSelectedClonotype] = useState(null);
  const [metadata, setMetadata] = useState([]);
  const [probabilities, setProbabilities] = useState([]);

  useEffect(() => {
    Promise.all([d3.tsv(metadataSource), d3.tsv(probabilitiesSource)]).then(
      data => {
        console.log(data[1]);
        setMetadata(data[0]);
        setProbabilities(data[1]);
      }
    );
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
          <LayoutWrapper
            data={{ metadata: metadata, probabilities: probabilities }}
          />
        </div>
      </DashboardProvider>
    </div>
  );
};

export default App;
