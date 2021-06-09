import FishtailComponent from "./Fishtail";
import data from "./data/timeseries.json";

const Template = (args) => <FishtailComponent {...args} />;

export default {
  title: "Components/TimeSeries/Fishtail",
  component: FishtailComponent,
};

export const Default = Template.bind({});
Default.args = { width: 700, height: 400, data, subsetParam: "clone" };
