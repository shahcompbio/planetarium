import { createMuiTheme } from "@material-ui/core/styles";
import "../App.css";
/*
  "MyFontBold",
  "MyFontRegular",
  "MyFontLight",
  "Helvetica",*/
export const theme = createMuiTheme({
  typography: {
    fontFamily: ["Noto Sans"].join(","),
  },
  palette: {
    primary: {
      main: "#95d2dc",
      dark: "#618ba0",
    },
    secondary: {
      main: "#f1c023",
    },
    error: {
      main: "#BC4746",
    },
    background: {
      default: "#F5F5F5",
    },
    overrides: {
      MuiFab: {
        root: {
          boxShadow: "none",
        },
      },
    },
  },

  spacing: 4,
});
