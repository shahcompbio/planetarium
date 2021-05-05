import React from "react";

import ProfileBackgroundComponent from "./ProfileBackground";

import { chromosomes, maxState } from "../data/cell_1.json";

const Template = (args) => <ProfileBackgroundComponent {...args} />;

const CHROMOSOMES = [
  { genomeStart: 1, length: 500, chr: "01" },
  { genomeStart: 501, length: 600, chr: "02" },
  { genomeStart: 1101, length: 300, chr: "03" },
  { genomeStart: 1401, length: 500, chr: "04" },
  { genomeStart: 1901, length: 900, chr: "X" },
];

export default {
  title: "Components/CopyNumber/Background",
  component: ProfileBackgroundComponent,
};

export const Background = Template.bind({});
Background.args = {
  width: 800,
  height: 500,
  chromosomes: CHROMOSOMES,
  maxState: 8,
};

export const Human = Template.bind({});
Human.args = {
  width: 800,
  height: 500,
  chromosomes: chromosomes,
  maxState: maxState,
};
