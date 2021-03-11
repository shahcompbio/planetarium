import React from "react";

import chartDim from "./LayoutConfig.js";
import componentList from "./ComponentList.js";

const App = props => {
  const Component = componentList[props.chartName];
  return <Component {...props} chartDim={{ ...chartDim, ...props.dim }} />;
};

export default App;
