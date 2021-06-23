import FishtailComponent from "./Fishtail";
import data from "./data/timeseries.json";

const Template = (args) => <FishtailComponent {...args} />;

export default {
  title: "Components/TimeSeries/Fishtail",
  component: FishtailComponent,
};

export const Default = Template.bind({});
Default.args = { width: 700, height: 400, data, subsetParam: "clone" };

export const Subset = Template.bind({});
Subset.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  subset: "A",
};

export const Timepoint = Template.bind({});
Timepoint.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  timepoint: "X8",
};

export const Static = Template.bind({});
Static.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  disable: true,
};
