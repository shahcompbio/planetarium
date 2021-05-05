import React from "react";

import GenomeProfileComponent from "./GenomeProfile";

import { bins, segs, chromosomes, maxState } from "../data/cell_1.json";

const Template = (args) => <GenomeProfileComponent {...args} />;

export default {
  title: "Components/CopyNumber/GenomeProfile",
  component: GenomeProfileComponent,
};

export const GenomeProfile = Template.bind({});
GenomeProfile.args = {
  width: 800,
  height: 500,
  maxState,
  bins,
  segs,
  chromosomes,
};
