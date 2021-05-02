import React from "react";
import _ from "lodash";

import Histogram from "./ProbabilityHistogram";

const Template = (args) => <Histogram {...args} />;

const DATA = [
  { prob: 0.05, obs: "a", group: "dog" },
  { prob: 0.05, obs: "a", group: "cat" },
  { prob: 0.05, obs: "a", group: "lion" },
  { prob: 0.15, obs: "b", group: "dog" },
  { prob: 0.05, obs: "c", group: "dog" },
  { prob: 0.05, obs: "c", group: "lion" },
  { prob: 0.2, obs: "d", group: "cat" },
  { prob: 0.2, obs: "d", group: "lion" },
  { prob: 0.01, obs: "e", group: "dog" },
  { prob: 0.7, obs: "f", group: "dog" },
];

export default {
  title: "Components/Bar/ProbabilityHistogram",
  component: Histogram,
};

export const Standard = Template.bind({});
Standard.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  subgroupParam: "group",
  observationParam: "obs",
};

export const SelectedObservation = Template.bind({});
SelectedObservation.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  subgroupParam: "group",
  observationParam: "obs",
  highlightedObservation: "a",
};

export const SelectedSubgroup = Template.bind({});
SelectedSubgroup.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  subgroupParam: "group",
  observationParam: "obs",
  highlightedSubgroup: "dog",
};
