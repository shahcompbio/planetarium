import React from "react";
import _ from "lodash";

import InfoBarComponent from "./InfoBar";

import { SearchIcon, DownloadIcon } from "./Icons";

import Avatar from "@mui/material/Avatar";
import SelectComponent from "../Select/VirtualizedSelect";
import { TextField } from "@mui/material";

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
    <SearchIcon>
      <TextField
        label="Dense"
        id="margin-dense"
        defaultValue="Default Value"
        margin="dense"
      />
    </SearchIcon>
  ),
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};

export const Multiple = Template.bind({});
Multiple.args = {
  addIcon: [
    <SearchIcon>
      <TextField
        label="Dense"
        id="margin-dense"
        defaultValue="Default Value"
        margin="dense"
      />
    </SearchIcon>,
    <DownloadIcon download={() => {}} />,
  ],
  title: "Info Bar",
  infoText:
    "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla nisl sem, convallis vitae risus sed, facilisis tempor nulla. Aliquam rutrum rhoncus ex, fringilla convallis metus gravida a. Donec bibendum varius convallis. Ut justo lectus, viverra vel dictum tristique, euismod eu odio. Phasellus ultricies lectus id ultricies facilisis. Nam suscipit arcu nibh, ac elementum lorem dictum in. Etiam vel tortor fringilla risus mattis scelerisque. Pellentesque habitant morbi tristique senectus et netus et malesuada fames ac turpis egestas. Praesent ultrices eleifend neque, at sollicitudin orci luctus ut. Cras libero leo, molestie ac erat non, condimentum aliquam leo. Interdum et malesuada fames ac ante ipsum primis in faucibus. Etiam consectetur risus sit amet egestas tincidunt. Vivamus laoreet pharetra nisl, vel rutrum nulla posuere at. Duis sed neque velit. Mauris dignissim, ex eleifend pulvinar eleifend, neque ligula eleifend velit, sit amet condimentum elit turpis vel felis. Donec sit amet gravida mi.",
};
