import React from "react";

import CellScapeComponent from "./CellScape";

import { segs, chromosomes } from "./data/all_seg.json";

const chromMap = chromosomes.reduce(
  (curr, chr) => ({ ...curr, [chr["chr"]]: chr["genomeStart"] }),
  {}
);

const data = segs.map((row) => ({
  ...row,
  segs: row["segs"].map((seg) => ({
    ...seg,
    genomeStart: chromMap[seg["chromosome"]] + seg["start"],
    length: seg["end"] - seg["start"] + 1,
  })),
}));

const Template = (args) => <CellScapeComponent {...args} />;

export default {
  title: "Components/CopyNumber/CellScape",
  component: CellScapeComponent,
};

export const CellScape = Template.bind({});
CellScape.args = { data, width: 800, height: 500, chromosomes };

export const Small = Template.bind({});
Small.args = {
  data: data.slice(0, 100),
  width: 800,
  height: 1000,
  chromosomes,
};
