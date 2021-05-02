import React, { useState, useEffect, useCallback, useRef } from "react";
import * as d3 from "d3";
import { withStyles } from "@material-ui/core/styles";

import { scalePoint } from "d3";
import {useCanvas} from "../../utils/useCanvas";



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
const getIndicesFromAllHeatmapOrder = (heatmapConfig, allHeatmapOrder) =>
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

const generateConfig = (width,height) =>{


const COLORS = [
  "#2e7aab",
  "#9ECAE1",
  "#CCCCCC",
  "#FDCC8A",
  "#FC8D59",
  "#E34A33",
  "#B30000",
  "#980043",
  "#DD1C77",
  "#DF65B0",
  "#C994C7",
  "#D4B9DA"
];

const CONSTANTS = {
  copyNumberLabels: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
  copyNumberColors: COLORS
};
const config = {
  ...CONSTANTS
};
  const COMPONENTS = {
    width: width,
    height: height,
  
    legendWidth: 100,
    legendX: 0,
  
    miniMapWidth: 100,
    minimapHeight: 400,
  
    heatmapChromHeight: 12
  };
  const componentConfig = {
    ...COMPONENTS,
    legendHeight: COMPONENTS.height,
  
    heatmapWidth:
      COMPONENTS.width - COMPONENTS.miniMapWidth - COMPONENTS.legendWidth,
    heatmapHeight: COMPONENTS.height,
  
    heatmapX: COMPONENTS.legendX + COMPONENTS.legendWidth
  };

  const LEGEND_CONSTANTS = {
    width: componentConfig.legendWidth,
    x: componentConfig.legendX,
    height: componentConfig.legendHeight
  };
  const legendConfig = {
    ...LEGEND_CONSTANTS,
    copyNumberLabels: CONSTANTS.copyNumberLabels,
    copyNumberColors: CONSTANTS.copyNumberColors
  };

const HEATMAP_CONSTANTS = {
  ...CONSTANTS,
  wrapperWidth: componentConfig.miniMapWidth + componentConfig.heatmapWidth,
  width: componentConfig.heatmapWidth,
  height: componentConfig.heatmapHeight,
  x: componentConfig.heatmapX,
  paddingLeft: 3,
  rowHeight:componentConfig.height/77,
  indicatorWidth: 10,
  defaultQuality: "0.75",

  minimap: {
    rowHeight: 0.75,
    width: componentConfig.miniMapWidth,
    brushWidth: componentConfig.miniMapWidth - 4,
    miniMapQueryRangeOffset: 2
  },

  profile: {
    y: componentConfig.heatmapHeight + componentConfig.heatmapChromHeight,
    x: 0,
    chromBoxHeight: 10,
    axisWidth: 25,
    axisTextYOffset: 5,
    axisLineWidth: 1,
    scaleDomain: [-0.5, 8],
    axisDomain: [0, 1, 2, 3, 4, 5, 6, 7],
    segmentColor: "#000000",
    barHeight: 2,
    backgroundColors: ["#fefefe", "#eee"],
    height: 300
  },

  chromosome: {
    height: componentConfig.heatmapChromHeight,
    color: ["#faf9f9", "#e6e6e6"]
  },

  indicator: {
    width: 15
  },

  legend: {
    height: 35,
    width: 300,
    squareSize: 10,
    squareSpacing: 5,
    textOffset: 5,
    titleWidth: 100,
    titleHeight: 15
  },
  categories: {
    legendWidth: 350,
    legendHeight: 55,
    lengendLineHeight: 15,
    squareSize: 6,
    lineSize: 2,
    squareSpacing: 2,
    colours: {
      0: ["#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#005a32"],
      1: [
        "#d9d9d9",
        "#bdbdbd",
        "#969696",
        "#737373",
        "#525252",
        "#252525",
        "#000000"
      ],
      2: ["#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#8c2d04"]
    }
  },
  qualitySliderMarks: [
    { value: 0, label: "0" },
    {
      value: 0.25,
      label: "0.25"
    },
    {
      value: 0.5,
      label: "0.5"
    },
    {
      value: 0.75,
      label: "0.75"
    },
    {
      value: 1.0,
      label: "1.0"
    }
  ]

};

const heatmapConfig = {
  ...HEATMAP_CONSTANTS,
  contentWidth: HEATMAP_CONSTANTS.width - HEATMAP_CONSTANTS.indicatorWidth,
  contentHeight: HEATMAP_CONSTANTS.height
};
return heatmapConfig

}

const CopyNumberHeatmap = ({width,height, analysis, allHeatmapOrder, categoryStats, chromosomes, segs,analysisStats }) => {

  const heatmapConfig = generateConfig(width,height)
  //  const [{ quality, selectedCells, subsetSelection }] = useStatisticsState();
  const [indices, setIndices] = useState([...getIndicesFromAllHeatmapOrder(heatmapConfig, allHeatmapOrder)]);
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
    
  
  
  const chromMap = getChromPixelMapping(heatmapConfig,chromosomes);
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
            <Grid container direction = "row" alignItems="center" >
              <Grid item xs = {2}>
                  <CategoriesLegend heatmapConfig= {heatmapConfig} choosenStats={categoryStats} />
              </Grid>
              <Grid item xs ={2}>
                    <Legend heatmapConfig= {heatmapConfig} maxState={analysisStats.maxState} />
              </Grid>
            </Grid>  
            <Grid container direction = "row">
              <Grid item >
                  <svg
                    width={categoryWidth + heatmapConfig.paddingLeft}
                    height={heatmapConfig["height"]-5}
                  >
                    <Categories
                      heatmapConfig= {heatmapConfig} 
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
                    heatmapConfig= {heatmapConfig} 
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
                        heatmapConfig= {heatmapConfig} 
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
  heatmapConfig,
  chromosomes,
  analysisStats,
  setHoverCellCoordinate,
  segs,
  categoryWidth,
}) => {

  const width =heatmapConfig["width"] + heatmapConfig["paddingLeft"]
  const height = heatmapConfig["height"] + heatmapConfig["rowHeight"];

  const [context, saveContext] = useState();

  const canvasRef = useCanvas((canvas) =>{
    const context = canvas.getContext("2d");
    drawHeatmap(segs, chromosomes, context);
  },width,height, [segs])

  const [rowHoverCordinates, setRowHoverCordinates] = useState(null);

  const yScale = getYScale(heatmapConfig.height / heatmapConfig.rowHeight);

  const invertYScale = d3.range(
    yScale.range()[0],
    yScale.range()[1],
    yScale.step()
  );

  const chromMap = getChromPixelMapping(heatmapConfig, chromosomes);
  
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

  //chromosomes, chromMap, categoryWidth 
  const drawHeatmap = (segs, chromosomes, context) => {

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

    context.clearRect(0, 0, width, height);
    let maxY = -1
    segs.forEach((segRow, index) => {
      const y = yScale(index);
      if (y >  maxY){
        maxY = y
      }
      const bpRatio = getBPRatio(heatmapConfig, chromosomes);
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
              <canvas ref = {canvasRef}/>
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
