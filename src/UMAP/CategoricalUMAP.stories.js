import React from "react";
import _ from "lodash";

import UMAPComponent from "./CategoricalUMAP";
import metadata from "./data/metadata.json";

const Template = (args) => <UMAPComponent {...args} />;

export default {
  title: "Components/UMAP/CategoricalUMAP",
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
  idParam: "cell_id",
};

export const SelectedSubtype = Template.bind({});
SelectedSubtype.args = {
  width: 800,
  height: 600,
  data: metadata,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  idParam: "cell_id",
  subsetParam: "subtype",
  subset: "Activated CD8",
};

export const Static = Template.bind({});
Static.args = {
  width: 800,
  height: 600,
  data: metadata,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  idParam: "cell_id",
  subsetParam: "subtype",
  disable: true,
};
