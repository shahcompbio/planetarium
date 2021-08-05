import React from "react";
import _ from "lodash";

import UMAPComponent from "./NumericalUMAP";
import data from "./data/cells.json";

const Template = (args) => <UMAPComponent data={data} {...args} />;

export default {
  title: "Components/UMAP/NumericalUMAP",
  component: UMAPComponent,
};

export const VIM = Template.bind({});
VIM.args = {
  width: 800,
  height: 600,
  xParam: "rna_UMAP_1",
  yParam: "rna_UMAP_2",
  subsetParam: "VIM",
  idParam: "cell_id",
};
