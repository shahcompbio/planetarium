import React from "react";

import ProfileComponent from "./Profile";

import { bins, segs, bptotal, maxState } from "../data/cell_1.json";

const Template = (args) => <ProfileComponent {...args} />;

const BINS = [
  {
    chr: "09",
    genomeStart: 1,
    length: 100,
    state: 2,
    copy: 2.4171,
  },
  {
    chr: "09",
    genomeStart: 101,
    length: 100,
    state: 2,
    copy: 1.8706,
  },
  {
    chr: "09",
    genomeStart: 201,
    length: 100,
    state: 3,
    copy: 2.7984,
  },
  {
    chr: "09",
    genomeStart: 301,
    length: 100,
    state: 3,
    copy: 3.0765,
  },
  {
    chr: "09",
    genomeStart: 401,
    length: 100,
    state: 3,
    copy: 2.6607,
  },
  {
    chr: "09",
    genomeStart: 501,
    length: 100,
    state: 3,
    copy: 2.5834,
  },
  {
    chr: "09",
    genomeStart: 601,
    length: 100,
    state: 3,
    copy: 3.2014,
  },
  {
    chr: "09",
    genomeStart: 701,
    length: 100,
    state: 2,
    copy: 2.2741,
  },
  {
    chr: "09",
    genomeStart: 801,
    length: 100,
    state: 2,
    copy: 2.2232,
  },
  {
    chr: "09",
    genomeStart: 901,
    length: 100,
    state: 2,
    copy: 2.6474,
  },
  {
    chr: "09",
    genomeStart: 1001,
    length: 100,
    state: 2,
    copy: 2.5307,
  },
];

const SEGS = [
  {
    chr: "09",
    genomeStart: 1,
    state: 2,
    length: 200,
  },
  {
    chr: "09",
    genomeStart: 201,
    state: 3,
    length: 500,
  },
  {
    chr: "09",
    genomeStart: 701,
    state: 2,
    length: 300,
  },
];

export default {
  title: "Components/CopyNumber/Profile",
  component: ProfileComponent,
};

export const Profile = Template.bind({});
Profile.args = {
  width: 100,
  height: 100,
  maxState: 8,
  bpTotal: 1001,
  bins: BINS,
  segs: SEGS,
};

export const Human = Template.bind({});
Human.args = {
  width: 800,
  height: 500,
  maxState: maxState,
  bpTotal: bptotal,
  bins: bins,
  segs: segs,
};
