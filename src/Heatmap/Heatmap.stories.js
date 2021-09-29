import React from "react";
import _ from "lodash";

import Heatmap from "./Heatmap";

const Template = (args) => <Heatmap {...args} />;

const DATA = [
  { prob: 0.05, obs: "a", group: "dog" },
  { prob: 0.05, obs: "a", group: "cat" },
  { prob: 0.05, obs: "a", group: "cat" },
  { prob: 0.05, obs: "a", group: "cat" },
  { prob: 0.05, obs: "a", group: "cat" },
  { prob: 0.05, obs: "a", group: "lion" },
  { prob: 0.15, obs: "b", group: "dog" },
  { prob: 0.05, obs: "c", group: "dog" },
  { prob: 0.05, obs: "c", group: "lion" },
  { prob: 0.2, obs: "d", group: "cat" },
  { prob: 0.2, obs: "d", group: "lion" },
  { prob: 0.2, obs: "d", group: "lion" },
  { prob: 0.01, obs: "e", group: "dog" },
  { prob: 0.7, obs: "f", group: "dog" },
  { prob: 0.7, obs: "f", group: "dog" },
  { prob: 0.7, obs: "f", group: "dog" },

  { prob: 0.05, obs: "a", group: "dog1" },
  { prob: 0.05, obs: "a", group: "cat1" },
  { prob: 0.05, obs: "a", group: "cat1" },
  { prob: 0.05, obs: "a", group: "cat1" },
  { prob: 0.05, obs: "a", group: "cat1" },
  { prob: 0.05, obs: "a", group: "lion1" },
  { prob: 0.15, obs: "b", group: "dog1" },
  { prob: 0.05, obs: "c", group: "dog1" },
  { prob: 0.05, obs: "c", group: "lion1" },
  { prob: 0.2, obs: "d", group: "cat1" },
  { prob: 0.2, obs: "d", group: "lion1" },
  { prob: 0.2, obs: "d", group: "lion1" },
  { prob: 0.01, obs: "e", group: "dog1" },
  { prob: 0.7, obs: "f", group: "dog1" },
  { prob: 0.7, obs: "f", group: "dog1" },
  { prob: 0.7, obs: "f", group: "dog1" },
];

export default {
  title: "Components/Heatmap/Heatmap",
  component: Heatmap,
};

export const Standard = Template.bind({});
Standard.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
};

export const SelectedRow = Template.bind({});
SelectedRow.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  highlightedRow: "dog",
};

export const SelectedColumn = Template.bind({});
SelectedColumn.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  highlightedColumn: "c",
};

export const SubsetColumns = Template.bind({});
SubsetColumns.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  columnLabels: ["a", "b", "c"].map((value) => ({
    value,
    label: `Label ${value}`,
  })),
};

export const RowTotal = Template.bind({});
RowTotal.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  rowTotal: { dog: 10, cat: 10, lion: 5, dog1: 20, cat1: 15, lion1: 10 },
};
export const FontChange = Template.bind({});
FontChange.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  highlightedRow: null,
  font: "Noto Sans",
  highlightedColumn: null,
  rowTotal: { dog: 10, cat: 10, lion: 5 },
};
export const Square = Template.bind({});
Square.args = {
  width: 600,
  height: 600,
  data: DATA,
  column: "obs",
  row: "group",
  highlightedRow: null,
  font: "Noto Sans",
  highlightedColumn: null,
  rowTotal: {
    dog: 10,
    cat: 10,
    lion: 5,
    dog1: 12,
    cat1: 15,
    lion1: 10,
    dog2: 14,
    cat2: 13,
    lion2: 9,
  },
};
