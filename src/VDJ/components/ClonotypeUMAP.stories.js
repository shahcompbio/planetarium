import React from "react";
import _ from "lodash";

import ClonotypeUMAP from "./Umap";

import metadata from "../data/metadata.json";

const CLONOTYPE_COLORS = [
  "#674172",
  "#098dde",
  "#fa832f",
  "#0e5702",
  "#c20c1e",
  "#911eb4",
  "#fc97bc",
  "#469990",
  "#b5762a",
  "#5aebed",
  "#8f8f3f",
  "#ed1a1a",
];

const clonotypeCounts = _.countBy(
  metadata.filter((datum) => datum["cdr3s_aa"] !== "None"),
  "cdr3s_aa"
);

const clonotypeLabels = Object.keys(clonotypeCounts)
  .sort((a, b) => clonotypeCounts[b] - clonotypeCounts[a])
  .slice(0, 10)
  .map((value, index) => ({
    value,
    label: `SEQ${index + 1} - ${value}`,
    color: CLONOTYPE_COLORS[index],
  }));

const Template = (args) => (
  <ClonotypeUMAP
    chartName={"UMAP"}
    data={metadata}
    setSelectedClonotype={() => {}}
    {...args}
  />
);

export default {
  title: "Dashboard/VDJ/UMAP",
  component: ClonotypeUMAP,
};

export const Clonotype = Template.bind({});
Clonotype.args = {
  chartDim: {
    width: 800,
    height: 600,
  },
  selectedClonotype: null,
  hoveredClonotype: null,
  clonotypeLabels,
};

export const SelectedClonotype = Template.bind({});
SelectedClonotype.args = {
  chartDim: {
    width: 800,
    height: 600,
  },
  selectedClonotype: "CTCSAEGGGTEVFF",
  hoveredClonotype: null,
  clonotypeLabels,
};
