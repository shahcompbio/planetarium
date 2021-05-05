import React from "react";

import ChromosomeAxisComponent from "./ChromosomeAxis";

import { chromosomes } from "./data/cell_1.json";

const Template = (args) => <ChromosomeAxisComponent {...args} />;

const CHROMOSOMES = [
  { genomeStart: 1, length: 500, chr: "01" },
  { genomeStart: 501, length: 600, chr: "02" },
  { genomeStart: 1101, length: 300, chr: "03" },
  { genomeStart: 1401, length: 500, chr: "04" },
  { genomeStart: 1901, length: 900, chr: "X" },
];

export default {
  title: "Components/CopyNumber/ChromosomeAxis",
  component: ChromosomeAxisComponent,
};

export const Sample = Template.bind({});
Sample.args = {
  width: 400,
  height: 20,
  chromosomes: CHROMOSOMES,
};

export const HumanChromosomes = Template.bind({});
HumanChromosomes.args = {
  width: 800,
  height: 20,
  chromosomes: chromosomes,
};
