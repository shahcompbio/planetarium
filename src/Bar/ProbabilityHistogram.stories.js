import React from "react";
import * as d3 from "d3";
import _ from "lodash";

import Histogram from "./ProbabilityHistogram";

const Template = (args) => <Histogram {...args} />;

const DATA = [
  { id: 0, prob: 0.05, obs: "a", group: "dog" },
  { id: 1, prob: 0.05, obs: "a", group: "cat" },
  { id: 2, prob: 0.05, obs: "a", group: "lion" },
  { id: 3, prob: 0.15, obs: "b", group: "dog" },
  { id: 4, prob: 0.05, obs: "c", group: "dog" },
  { id: 5, prob: 0.05, obs: "c", group: "lion" },
  { id: 7, prob: 0.2, obs: "d", group: "cat" },
  { id: 6, prob: 0.2, obs: "d", group: "lion" },
  { id: 8, prob: 0.01, obs: "e", group: "dog" },
  { id: 9, prob: 0.7, obs: "f", group: "dog" },
];
const format = d3.format(".3f");
const getTooltipText = (bin) => {
  const subtypeGroups = _.groupBy(bin, "group");

  return (
    <div>
      {Object.keys(subtypeGroups)
        .sort((a, b) => a.length - b.length)
        .map((group) => (
          <p>
            {group}: {format(subtypeGroups[group].length / bin.length) * 100}%
          </p>
        ))}
    </div>
  );
};

export default {
  title: "Components/Bar/ProbabilityHistogram",
  component: Histogram,
};

// Test

export const Standard = Template.bind({});
Standard.args = {
  data: DATA,
  width: 400,
  height: 400,
  probParam: "prob",
  observationParam: "obs",
  idParam: "id",
  highlightedObservation: null,
  highlightedIDs: null,
  getTooltipText,
};

export const HighlightObservation = Template.bind({});
HighlightObservation.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  observationParam: "obs",
  highlightedObservation: "a",
};

export const HighlightIds = Template.bind({});
HighlightIds.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  highlightedIDs: DATA.filter((datum) => datum["group"] === "dog").map(
    (datum) => datum["id"]
  ),
  observationParam: "obs",

  getTooltipText,
};
export const FontChange = Template.bind({});
FontChange.args = {
  width: 400,
  height: 400,
  data: DATA,
  probParam: "prob",
  highlightedIDs: DATA.filter((datum) => datum["group"] === "dog").map(
    (datum) => datum["id"]
  ),
  observationParam: "obs",
  font: "MyFontRegular",
  getTooltipText,
};
