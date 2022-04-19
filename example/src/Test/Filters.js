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

export const removeFilters = (activeFilters, type, value) => {
  if (activeFilters[type].length === 1) {
    delete activeFilters[type];
  } else {
    const index = activeFilters[type].indexOf(value);
    activeFilters[type].splice(index, 1);
  }
  return activeFilters;
};

export const setMainStateWithNewFilters = (
  layerFilters,
  activeLayer,
  newFilters,
  setLayerFilters
) => {
  var newLayerFilters = layerFilters;
  newLayerFilters[activeLayer]["filters"] = newFilters;
  console.log(newLayerFilters);
  setLayerFilters(newLayerFilters);
};

const Filters = ({ filters, setLayerFilters, layerFilters, activeLayer }) => {
  const [checked, setChecked] = useState([-1]);

  const [activeFilters, setActiveFilters] = useState({});

  useEffect(() => {
    const newActiveFilter = layerFilters[activeLayer]["filters"];
    const filter = newActiveFilter ? newActiveFilter : {};
    const oldChecked = newActiveFilter
      ? Object.keys(newActiveFilter).reduce(
          (final, heading) => {
            final = [...final, ...filter[heading]];
            return final;
          },
          [-1]
        )
      : [-1];
    setChecked([-1, ...oldChecked]);
    setActiveFilters({ ...filter });
  }, [activeLayer]);

  const handleToggle = (value, name) => () => {
    var newFilters = activeFilters;
    if (value) {
      const currentIndex = checked.indexOf(value);
      const newChecked = [...checked];

      if (currentIndex === -1) {
        newChecked.push(value);
        newFilters = addToFilters(newFilters, name, value);
        console.log(newFilters);
      } else {
        newChecked.splice(currentIndex, 1);
        newFilters = removeFilters(newFilters, name, value);
        console.log(newFilters);
      }
      setMainStateWithNewFilters(
        layerFilters,
        activeLayer,
        newFilters,
        setLayerFilters
      );
      setChecked(newChecked);
    }
  };

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
