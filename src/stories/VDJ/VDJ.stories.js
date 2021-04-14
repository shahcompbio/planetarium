import React from "react";

import { VDJ as VDJApp } from "../../VDJ/VDJ";

import probabilities from "./data/probabilities.json";
import metadata from "./data/metadata.json";
import degs from "./data/deg.json";

const Template = (args) => <VDJApp {...args} />;

export default {
  title: "Dashboard/VDJ/VDJ",
  component: VDJApp,
};

export const VDJ = Template.bind({});
VDJ.args = {
  probabilities,
  metadata,
  degs,
};
