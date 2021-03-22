import React from "react";

import NDVApp from "../NDV/NDV";
import api from "../NDV/data/api";

const DevApp = () => {
  const data = api();

  return Object.keys(data).length === 0 ? null : <NDVApp data={data} />;
};

const Template = (args) => <DevApp {...args} />;

export default {
  title: "Dashboard",
  component: DevApp,
};

export const NDV = Template.bind({});
