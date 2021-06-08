import FishtailComponent from "./Fishtail";
import data from "./data/timeseries.json";

// const DATA = [
//   { timepoint: 1, subset: "A", count: 100 },
//   { timepoint: 2, subset: "A", count: 50 },
//   { timepoint: 3, subset: "A", count: 5 },
//   { timepoint: 4, subset: "A", count: 20 },
//   { timepoint: 2, subset: "B", count: 50 },
//   { timepoint: 3, subset: "B", count: 200 },
//   { timepoint: 4, subset: "B", count: 5 },
//   { timepoint: 2, subset: "C", count: 50 },
//   { timepoint: 2, subset: "D", count: 50 },
//   { timepoint: 3, subset: "D", count: 50 },
//   { timepoint: 4, subset: "D", count: 50 },
// ];

// const DATA2 = [
//   { timepoint: "X1", A: 100 },
//   { timepoint: "X2", A: 50, B: 50, C: 50, D: 50 },
//   { timepoint: "X3", A: 5, B: 200, D: 50 },
//   { timepoint: "X4", A: 20, B: 5, D: 50 },
// ];

const Template = (args) => <FishtailComponent {...args} />;

export default {
  title: "Components/TimeSeries/Fishtail",
  component: FishtailComponent,
};

export const Default = Template.bind({});
Default.args = { width: 700, height: 400, data, subsetParam: "clone" };
