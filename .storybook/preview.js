import React from "react";

import { theme } from "../src/theme/theme";
import { MuiThemeProvider } from "@material-ui/core/styles";

export const decorators = [
  (Story) => (
    <MuiThemeProvider theme={theme}>
      <Story />
    </MuiThemeProvider>
  ),
];

export const parameters = {
  actions: { argTypesRegex: "^on[A-Z].*" },
};
