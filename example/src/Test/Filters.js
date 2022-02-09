import React, { useState, useEffect } from "react";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemButton from "@mui/material/ListItemButton";
import ListItemText from "@mui/material/ListItemText";
import ListSubheader from "@mui/material/ListSubheader";
import Checkbox from "@mui/material/Checkbox";

const addToFilters = (activeFilters, name, value) => {
  if (activeFilters.hasOwnProperty(name)) {
    activeFilters[name] = [...activeFilters[name], value];
  } else {
    activeFilters[name] = [value];
  }
  return activeFilters;
};
const removeFilters = (activeFilters, name, value) => {
  if (activeFilters[name].length === 1) {
    delete activeFilters[name];
  } else {
    const index = activeFilters[name].indexOf(value);
    activeFilters[name].splice(index, 1);
  }
  return activeFilters;
};
const Filters = ({ filters, setLayerFilters, layerFilters, activeLayer }) => {
  const [checked, setChecked] = useState(null);

  const [activeFilters, setActiveFilters] = useState({});
  const handleToggle = (value, name) => () => {
    var newFilters = activeFilters;
    if (checked !== null) {
      const currentIndex = checked.indexOf(value);
      const newChecked = [...checked];

      if (currentIndex === -1) {
        newChecked.push(value);
        newFilters = addToFilters(newFilters, name, value);
      } else {
        newChecked.splice(currentIndex, 1);
        newFilters = removeFilters(newFilters, name, value);
      }
      setActiveFilters(newFilters);
      setChecked(newChecked);
    } else {
      newFilters = addToFilters(newFilters, name, value);
      setActiveFilters(newFilters);

      setChecked([value]);
    }
  };

  useEffect(() => {
    if (activeFilters) {
      var newLayerFilters = layerFilters;
      newLayerFilters[activeLayer]["filters"] = activeFilters;
      setLayerFilters(newLayerFilters);
    }
  }, [activeFilters, setLayerFilters]);

  return (
    <List
      dense
      sx={{ width: "100%", maxWidth: 360, bgcolor: "background.paper" }}
    >
      {filters.map((filter) => {
        const name = filter["name"];

        return (
          <div key={"wrapper-" + name}>
            <ListSubheader key={"subheader-" + name}>{name}</ListSubheader>
            {filter["values"].map((item) => {
              const labelId = `checkbox-list-secondary-label-${item}`;
              return (
                <ListItem
                  key={"listitem-" + item}
                  onClick={handleToggle(item, name)}
                  secondaryAction={
                    <Checkbox
                      edge="end"
                      key={"check-" + item}
                      onChange={handleToggle(item, name)}
                      checked={checked && checked.indexOf(item) !== -1}
                      inputProps={{ "aria-labelledby": item }}
                    />
                  }
                  disablePadding
                >
                  <ListItemButton key={"listitembutton-" + item}>
                    <ListItemText
                      key={"listitemtext-" + item}
                      id={labelId}
                      primary={item}
                    />
                  </ListItemButton>
                </ListItem>
              );
            })}
          </div>
        );
      })}
    </List>
  );
};
export default Filters;
