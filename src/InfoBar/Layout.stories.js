import React from "react";
import _ from "lodash";

import { TextField } from "@material-ui/core";

import InfoLayoutComponent from "./Layout";

const Template = (args) => (
  <InfoLayoutComponent {...args}>
    <p>Some test content here... nothing to see!</p>
  </InfoLayoutComponent>
);

export default {
  title: "Components/Info/Layout",
  component: InfoLayoutComponent,
};

export const Layout = Template.bind({});
Layout.args = {
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};

export const LayoutWithSearch = Template.bind({});
LayoutWithSearch.args = {
  addIcon: (
    <TextField
      label="Dense"
      id="margin-dense"
      defaultValue="Default Value"
      margin="dense"
    />
  ),
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};
