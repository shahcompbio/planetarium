import React from "react";
import _ from "lodash";

import Bar from "./StackedHorizontalBar";

const Template = (args) => <Bar {...args} />;

const LABELS = ["small", "medium", "large"];
const DATA = {
  dog: { small: 50, medium: 25, large: 30 },
  cat: { small: 80, medium: 20, large: 0 },
  mouse: { small: 60, medium: 0, large: 10 },
};
export default {
  title: "Components/Bar/StackedHorizontal",
  component: Bar,
};

export const Standard = Template.bind({});
Standard.args = {
  width: 400,
  height: 400,
  data: DATA,
  barLabels: LABELS,
  highlightedRow: null,
};

export const SelectedRow = Template.bind({});
SelectedRow.args = {
  width: 400,
  height: 400,
  data: DATA,
  barLabels: LABELS,
  highlightedRow: "dog",
};

export const FontChange = Template.bind({});
FontChange.args = {
  width: 400,
  height: 400,
  data: DATA,
  barLabels: LABELS,
  font: "MyFontRegular",
  highlightedRow: "dog",
};
