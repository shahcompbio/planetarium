import React from "react";

import chartDim from "./LayoutConfig.js";
import componentList from "./ComponentList.js";

const App = ({
  chartName,
  data,
  setSelectedSubtype,
  setSelectedClonotype,
  selectedSubtype,
  selectedClonotype
}) => {
  const Component = componentList[chartName];
  return <Component data={data} chartDim={chartDim} />;
};

export default App;
