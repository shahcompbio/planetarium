import React, { useState, useEffect, useCallback, useRef } from "react";
import * as d3 from "d3";
import { withStyles } from "@material-ui/core/styles";

import { scalePoint } from "d3";
import { heatmapConfig } from "./config.js";
import { initContext } from "./utils.js";
import {
  getYScale,
  getChromPixelMapping,
  getBPRatio,
  getSegWidth,
  getSegX,
  colorScale,
  drawChromosomeAxis
} from "./utils.js";

import CategoriesLegend from "./CategoriesLegend.js";
import Categories from "./Categories.js";
import Indicator from "./Indicator.js";
import Legend from "./Legend.js";
import Minimap from "./Minimap.js";
import {Grid} from "@material-ui/core";


//const data = require('../test_data/CopyNumberHeatmap/copyNumberHeatmapData.json')

const width = heatmapConfig["width"] + heatmapConfig["paddingLeft"];
const height = heatmapConfig["height"] + heatmapConfig["rowHeight"];

const margin = {
  left: 75,
  top: 37,
  bottom: 90,
  right: 10,
  histogram: 20,
};
const styles = (theme) => ({
  content: {
    flexGrow: 1,
    backgroundColor: "#FFFFFFF",
    padding: theme.spacing.unit * 3,
  },
  container: {
    minHeight: "100vh",
  },
});
const getIndicesFromAllHeatmapOrder = allHeatmapOrder =>
  allHeatmapOrder.filter(
    (order, index) => index < heatmapConfig.height / heatmapConfig.rowHeight - 2
  );


const getSelectedSegsFromIndices = (segs, indices) => {
  return segs.filter(
    seg => seg.index >= indices[0] && seg.index <= indices[indices.length-1]
  )
}

const getSelectedAnalysisCellStatsFromIndices = (cellStats, indices)=>{
  return cellStats.filter(
    stat => stat.heatmap_order >= indices[0] && stat.heatmap_order <= indices[indices.length-1]
    )
}

const CopyNumberHeatmap = ({ analysis, allHeatmapOrder, categoryStats, chromosomes, segs,analysisStats }) => {

  //  const [{ quality, selectedCells, subsetSelection }] = useStatisticsState();
  const [indices, setIndices] = useState([...getIndicesFromAllHeatmapOrder(allHeatmapOrder)]);
  const [hoverCell, setHoverCell] = useState({ cell: {} });
  const [selectedCell, setSelectedCell] = useState({ cell: {} });
  const heatmapOrderToHeatmapIndex = scalePoint()
    .domain([...allHeatmapOrder])
    .range([0, allHeatmapOrder.length - 1]);

  const categoryWidth =
      categoryStats.length * heatmapConfig.categories.squareSize +
      categoryStats.length * heatmapConfig.categories.squareSpacing;

  const yScale = getYScale(
      heatmapConfig.height / heatmapConfig.rowHeight
    );

  const chromMap = getChromPixelMapping(chromosomes);
  const selectedSegs = getSelectedSegsFromIndices(segs, indices)
  const selectedAnalysisStats= {"maxState":analysisStats.maxState,"cellStats":getSelectedAnalysisCellStatsFromIndices(analysisStats.cellStats, indices)}
  return (
          <Grid
            item
            container
            direction="row"
            style={{ position: "relative" }}
            width={heatmapConfig.wrapperWidth}
            height={
              heatmapConfig["height"] - heatmapConfig.chromosome["height"]
            }
          >
            <Grid container direction = "row" >
              <Grid item>
                  <CategoriesLegend choosenStats={categoryStats} />
              </Grid>
              <Grid item>
                    <Legend maxState={analysisStats.maxState} />
              </Grid>
            </Grid>  
            <Grid container direction = "row" >
              <Grid item >
                  <svg
                    width={categoryWidth + heatmapConfig.paddingLeft}
                    height={heatmapConfig["height"]-5}
                  >
                    <Categories
                      categories={categoryStats}
                      cellStats={selectedAnalysisStats.cellStats}
                      yScale={yScale}
                    />

                    {hoverCell.hasOwnProperty("y") && (
                      <Indicator y={hoverCell["y"]} />
                    )}
                  </svg>
              </Grid> 
              <Grid item >
                  <Plot
                    setHoverCellCoordinate={(y, heatmapRow) => {
                      const cell = segs[heatmapRow];
                      if (cell !== undefined) {
                        setHoverCell({
                          y: yScale(heatmapRow),
                          cell: segs[heatmapRow],
                        });
                      }
                    }}
                    chromosomes={chromosomes}
                    analysisStats={selectedAnalysisStats}
                    segs={selectedSegs}
                    heatmapOrderToHeatmapIndex={heatmapOrderToHeatmapIndex}
                    categoryWidth={categoryWidth + heatmapConfig.paddingLeft}
                  />
              </Grid>
              <Grid item>
                  <Minimap
                        triggerHeatmapRequery={(index) => setIndices([...index])}
                        heatmapOrder={allHeatmapOrder}
                        rangeExtent={[
                          heatmapOrderToHeatmapIndex(indices[0]),
                          heatmapOrderToHeatmapIndex(indices[indices.length - 1]),
                        ]}
                        analysis={analysis}
                        chromosomes={chromosomes}
                        chromMap={chromMap}
                        segs = {segs}
                      />
              </Grid>

            </Grid>  
          </Grid>

  );
};

const Plot = ({
  chromosomes,
  analysisStats,
  setHoverCellCoordinate,
  segs,
  categoryWidth,
}) => {
  const [context, saveContext] = useState();

  const [ref] = useHookWithRefCallback();
  const [rowHoverCordinates, setRowHoverCordinates] = useState(null);

  const yScale = getYScale(heatmapConfig.height / heatmapConfig.rowHeight);

  const invertYScale = d3.range(
    yScale.range()[0],
    yScale.range()[1],
    yScale.step()
  );

  const chromMap = getChromPixelMapping(chromosomes);

  useEffect(() => {
    if (rowHoverCordinates !== null) {
      const roundY = Math.max(rowHoverCordinates[1], 0);
      var heatmapRow = yScale.domain()[d3.bisect(invertYScale, roundY) - 1];
      setHoverCellCoordinate(rowHoverCordinates[1], heatmapRow);
    }
  }, [rowHoverCordinates]);

  useEffect(() => {
    if (context) {
      drawHeatmap(segs,chromosomes, context);
    }
  }, [segs]);


  function useHookWithRefCallback() {
    const ref = useRef(null);
    const setRef = useCallback((node) => {
      if (node) {
        const heatmap = d3.select("#heatmap");
        const canvas = heatmap
          .select("canvas")
          .attr("width", width)
          .attr("height", height - heatmapConfig.chromosome["height"])
          .attr("transform", "translate(" + 0 + "," + margin.top + ")");
        const context = initContext(canvas, width, height);

        saveContext(context);
        context.save();
        drawHeatmap(segs, chromosomes,context);
        d3.select("#heatSelection").on("mousemove", function() {
          var coordinates = d3.mouse(this);
          var alreadySelected = d3.select(this).attr("class");
          if (
            alreadySelected === null ||
            parseInt(alreadySelected) !== coordinates[1]
          ) {
            setRowHoverCordinates(coordinates);
            d3.select(this).attr("class", coordinates[1]);
          }
        });
      }
    }, []);

    return [setRef];
  }
  
 
  //chromosomes, chromMap, categoryWidth 
  const drawHeatmap = (segs, chromosomes, context) => {
    context.clearRect(0, 0, width, height);
    let maxY = -1
    segs.forEach((segRow, index) => {
      const y = yScale(index);
      if (y >  maxY){
        maxY = y
      }
      const bpRatio = getBPRatio(chromosomes);
      segRow.segs.forEach((seg) => {
        const x = getSegX(seg, chromMap, bpRatio, false, 0);
        context.fillStyle = colorScale(seg["state"]);
        context.fillRect(
          x,
          y,
          getSegWidth(seg, bpRatio),
          heatmapConfig["rowHeight"]
        );
        context.stroke();
      });
    
    });

    drawChromosomeAxis(context, chromosomes,chromMap, categoryWidth,maxY)
  };

  return (
    <div
      style={{
        width: width - categoryWidth,
        height: height - heatmapConfig.chromosome["height"],
        position: "relative",
      }}
      ref={ref}
    >
      <Grid container direction = "row" >
        <Grid item>
              <div
              id="heatmap"
              style={{
                width: width - categoryWidth,
                height: height - heatmapConfig.chromosome["height"],
                position: "absolute",
                pointerEvents: "all",
              }}
            >
              <canvas />
            </div>
            <svg
              id="heatSelection"
              style={{
                width: width - categoryWidth,
                height: height - heatmapConfig.chromosome["height"],
                position: "relative",
              }}
            />
        </Grid>
      </Grid>
    </div>
  );
};

export default withStyles(styles)(CopyNumberHeatmap);

