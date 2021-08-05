import React from "react";
import _ from "lodash";

import VirtualizedSelectComponent from "./VirtualizedSelect";
import genes from "./data/genes.json";

const Template = (args) => (
  <VirtualizedSelectComponent options={genes} {...args} />
);

export default {
  title: "Components/Select/Virtualized",
  component: VirtualizedSelectComponent,
};

export const Default = Template.bind({});
Default.args = {
  value: genes[0],
  title: "Genes",
};

export const SmallWidth = Template.bind({});
SmallWidth.args = {
  value: genes[0],
  width: 200,
  title: "Genes",
};

export const LargeWidth = Template.bind({});
LargeWidth.args = {
  value: genes[0],
  width: 500,
  title: "Genes",
};
