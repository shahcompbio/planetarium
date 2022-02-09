import React from "react";
import Card from "@mui/material/Card";
import CardActions from "@mui/material/CardActions";
import CardContent from "@mui/material/CardContent";
import Typography from "@mui/material/Typography";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import ListSubheader from "@mui/material/ListSubheader";

const LayersMenu = ({ layerName }) => {
  return (
    <Card sx={{ minWidth: 275, marginBottom: 5 }} key={layerName + "card"}>
      <CardContent>
        <Typography variant="h5" component="div">
          {layerName}
        </Typography>
      </CardContent>
      <CardActions>
        <List
          sx={{
            width: "100%",
            maxWidth: 360,
            bgcolor: "background.paper",
            position: "relative",
            overflow: "auto",
            maxHeight: 300,

            "& ul": { padding: 0 },
          }}
          subheader={<li />}
        >
          {["Patients", "Subtype"].map((sectionId) => (
            <li key={`section-${sectionId}`}>
              <ul>
                <ListSubheader>{sectionId}</ListSubheader>
                {[0, 1].map((item) => (
                  <ListItem key={`item-${sectionId}-${item}`}>
                    <ListItemText primary={`Item ${item}`} />
                  </ListItem>
                ))}
              </ul>
            </li>
          ))}
        </List>
      </CardActions>
    </Card>
  );
};
export default LayersMenu;
