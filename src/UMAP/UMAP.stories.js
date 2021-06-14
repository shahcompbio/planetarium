import React from "react";
import _ from "lodash";

import UMAPComponent from "./UMAP";
import metadata from "./data/metadata.json";

const Template = (args) => <UMAPComponent {...args} />;

export default {
  title: "Components/UMAP/UMAP",
  component: UMAPComponent,
};

export const Subtype = Template.bind({});
Subtype.args = {
  width: 800,
  height: 600,
  data: metadata,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
};

export const SelectedSubtype = Template.bind({});
SelectedSubtype.args = {
  width: 800,
  height: 600,
  data: metadata,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
  highlighted: "Activated CD8",
};
