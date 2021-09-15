import React from "react";
import _ from "lodash";

import InfoBarComponent from "./InfoBar";

import SelectComponent from "../Select/VirtualizedSelect";
import { TextField } from "@material-ui/core";
import SvgIcon from "@material-ui/core/SvgIcon";
import SearchIcon from "@material-ui/icons/Search";

const Template = (args) => <InfoBarComponent {...args} />;

export default {
  title: "Components/Info/Header",
  component: InfoBarComponent,
};

export const Header = Template.bind({});
Header.args = {
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};

export const Select = Template.bind({});
Select.args = {
  addIcon: (
    <SelectComponent
      options={["A", "B", "C", "D"]}
      value={"A"}
      title="Number"
    />
  ),
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};

export const Search = Template.bind({});
Search.args = {
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

export const Multiple = Template.bind({});
Multiple.args = {
  addIcon: [
    <SvgIcon
      viewBox="0 0 20 20"
      fontSize="small"
      style={{ fontSize: 16, marginLeft: 10 }}
    >
      <SearchIcon />
    </SvgIcon>,
    <SvgIcon
      viewBox="0 0 20 20"
      fontSize="small"
      style={{ fontSize: 16, marginLeft: 10 }}
    >
      <path d="M 6.089844 8.722656 L 11.378906 3.433594 L 12.847656 4.902344 L 7.558594 10.191406 Z M 6.089844 8.722656 " />
      <path d="M 6.375 0.429688 L 8.757812 0.429688 L 8.757812 7.953125 L 6.375 7.953125 Z M 6.375 0.429688 " />
      <path d="M 7.566406 10.195312 L 2.277344 4.910156 L 3.742188 3.441406 L 9.035156 8.730469 Z M 7.566406 10.195312 " />
      <path d="M 1.125 12.445312 L 14 12.445312 L 14 14.570312 L 1.125 14.570312 Z M 1.125 12.445312 " />
      <path d="M 0.6875 8.570312 L 2.8125 8.570312 L 2.8125 14.570312 L 0.6875 14.570312 Z M 0.6875 8.570312 " />
      <path d="M 12.1875 8.570312 L 14.3125 8.570312 L 14.3125 14.570312 L 12.1875 14.570312 Z M 12.1875 8.570312 " />
    </SvgIcon>,
  ],
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};
