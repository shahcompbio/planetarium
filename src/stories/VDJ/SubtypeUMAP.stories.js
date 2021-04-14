import React from "react";
import _ from "lodash";

import SubtypeUMAP from "../../VDJ/components/subTypeUmap";

import metadata from "./data/metadata.json";

const Template = (args) => (
  <SubtypeUMAP
    chartName={"SUBTYPEUMAP"}
    data={metadata}
    setSelectedSubtype={() => {}}
    {...args}
  />
);

export default {
  title: "Dashboard/VDJ/UMAP",
  component: SubtypeUMAP,
};

export const Subtype = Template.bind({});
Subtype.args = {
  chartDim: {
    width: 800,
    height: 600,
  },
  selectedSubtype: null,
  hoveredSubtype: null,
};

export const SelectedSubtype = Template.bind({});
SelectedSubtype.args = {
  chartDim: {
    width: 800,
    height: 600,
  },
  selectedSubtype: "Activated CD8",
  hoveredSubtype: null,
};
