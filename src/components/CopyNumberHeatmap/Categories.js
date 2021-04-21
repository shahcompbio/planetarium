import React, { useEffect } from "react";

import { scaleOrdinal } from "d3";
import * as d3 from "d3";
import { heatmapConfig } from "./config.js";
import { cleanUpPreviousContent } from "./utils.js";

const getColourScale = (heatmapConfig,types, index) =>
  scaleOrdinal()
    .domain(types)
    .range(heatmapConfig.categories.colours[index]);

const getCategoryWidth = (heatmapConfig,categoriesLength) =>
  categoriesLength * heatmapConfig.categories.squareSize;

const Categories = ({ heatmapConfig, cellStats, yScale, categories }) => {
  const squareSize = heatmapConfig.categories.squareSize;
  useEffect(() => {
    if (cellStats) {
      const categoriesWrapper = d3.select("#categories");
      cleanUpPreviousContent(categoriesWrapper);
      const categoryWidth = getCategoryWidth(heatmapConfig,categories.length);
      categories.forEach((category, index) => {
        const categoryName = category.category;
        const colourScale = getColourScale(heatmapConfig,category.types, index);

        categoriesWrapper
          .selectAll(".category-" + categoryName)
          .data(cellStats)
          .enter()
          .append("rect")
          .attr("id", function(d) {
            return d.id + categoryName;
          })
          .attr("class", "cat-square category-" + categoryName)
          .attr("x", xCordinate)
          .attr("y", function(d, i) {
            return yScale(i);
          })
          .attr("width", squareSize)
          .attr("height", squareSize)
          .attr("fill", function(d) {
            return colourScale(d[categoryName]);
          });

        function xCordinate(d, i) {
          const defaultSpacing = index * heatmapConfig.categories.squareSize;
          const xOffset = categories.length === 1 ? 4 : 0;

          return index === 0
            ? categoryWidth - defaultSpacing - xOffset
            : categoryWidth -
                (defaultSpacing +
                  index * heatmapConfig.categories.squareSpacing);
        }
      });
    }
  }, [categories, cellStats, squareSize, yScale]);

  return <g id="categories"></g>;
};
export default Categories;
