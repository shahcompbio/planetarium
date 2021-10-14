import React from "react";

import { theme } from "../src/theme/theme";
import { ThemeProvider, StyledEngineProvider } from "@mui/material/styles";
import { ThemeProvider as EmotionThemeProvider } from "emotion-theming";

export const decorators = [
  (Story) => (
    <EmotionThemeProvider theme={theme}>
      <ThemeProvider theme={theme}>
        <Story />
      </ThemeProvider>
    </EmotionThemeProvider>
  ),
];

export const parameters = {
  actions: { argTypesRegex: "^on[A-Z].*" },
};
