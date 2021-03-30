import { config, heatmapConfig } from "./config.js";
import * as d3  from "d3";
import { useRef, useEffect } from "react";


export const cleanUpPreviousContent = wrapper =>
  wrapper.selectAll("*").remove();
/**
 * Returns segment starting x position
 * @param {object} seg
 * @param {object} chromMap
 * @param {number} bpRatio
 * @param {number}
 */
export const getSegX = (seg, chromMap, bpRatio, isMinimap, categoryWidth) =>
  isMinimap
    ? Math.floor(seg.start / bpRatio) + chromMap[seg.chromosome].miniX
    : Math.floor(seg.start / bpRatio) +
      chromMap[seg.chromosome].x +
      categoryWidth;

/**
 * Returns segment width in pixels
 * @param {object} seg
 * @param {number} bpRatio
 */
export const getSegWidth = (seg, bpRatio) =>
  Math.floor((seg.end - seg.start + 1) / bpRatio);

export const getChromPixelMapping = chromosomes => {
  const bpRatio = getBPRatio(chromosomes);
  const miniBpRatio = getMinimapBPRatio(chromosomes);

  let xShift = 0;
  let miniXShift = 0;

  return chromosomes.reduce((map, chrom) => {
    const chromWidth = getChromWidth(chrom, bpRatio);
    const miniWidth = getChromWidth(chrom, miniBpRatio);

    const mapEntry = {
      chrom: chrom,
      x: xShift,
      miniX: miniXShift,
      miniWidth: miniWidth,
      width: chromWidth
    };

    miniXShift += miniWidth;
    xShift += chromWidth;

    return {
      ...map,
      [chrom.id]: mapEntry
    };
  }, {});
};
export const getGenomeYScale = maxState =>
  d3.scaleLinear()
    .domain([-0.5, maxState])
    .range([heatmapConfig.profile.height, 0]);

export const getMinimapBPRatio = chromosomes => {
  const totalBP = getTotalBP(chromosomes);
  return Math.ceil(totalBP / heatmapConfig.minimap.width);
};
/**
 * Gets base pair to pixel ratio
 */
export const getBPRatio = chromosomes => {
  const totalBP = getTotalBP(chromosomes);
  return Math.ceil(totalBP / heatmapConfig["contentWidth"]);
};

/**
 * Gets number of indices that can fit per heatmap row
 */
export const getIndicesPerRow = Math.ceil(config["rowHeight"]);

/**
 * Gets the total number of base pairs in chromosome ranges
 */
export const getTotalBP = chromosomes => {
  const totalBp = chromosomes.reduce(
    (sum, chrom) => sum + chrom.end - chrom.start + 1,
    0
  );
  return totalBp;
};
export const colorScale = d3.scaleLinear()
  .domain(config["copyNumberLabels"])
  .range(config["copyNumberColors"]);
/**
 * Returns the width (in pixels) for chromosome
 * @param {object} chrom - data
 * @param {int} bpRatio
 * @return {int}
 */
const getChromWidth = (chrom, bpRatio) =>
  Math.floor((chrom.end - chrom.start + 1) / bpRatio);

export const getYScale = indicesLength =>
  d3.scalePoint()
    .domain([0, ...Array.from(Array(indicesLength - 1).keys())])
    .range([0, heatmapConfig.rowHeight * (indicesLength - 1)]);

export const heatmapOrderToCellIndex = (heatmapOrder, cellCount) =>
  d3.scalePoint()
    .domain([...heatmapOrder])
    .range([0, cellCount - 1]);

export const getMinimapYScale = indicesLength => {
  return d3.scalePoint()
    .domain([...Array.from(Array(indicesLength).keys())])
    .range([0, indicesLength * heatmapConfig.minimap.rowHeight]);
};
export const heatmapToMinimapScale = (heatmapCellLength, minimapCellLength) =>
  d3.scaleLinear()
    .domain([0, heatmapCellLength - 1])
    .range([0, heatmapConfig.minimap.rowHeight * minimapCellLength - 1]);

export const getChromosomeEndX = chromosomes => {
  const chromPixelMapping = getChromPixelMapping(chromosomes);
  return chromPixelMapping["Y"]["width"] + chromPixelMapping["Y"]["x"];
};
export const getIndicatorXPosition = chromosomes => {
  const annotationX = getAnnotationsX(chromosomes);
  return annotationX;
};

export const getAnnotationsX = chromosomes => {
  const chromosomeX = getChromosomeEndX(chromosomes);
  const heatmapX = 0;
  return chromosomeX + heatmapX + config["spacing"];
};

export const initContext = (canvasSelection, width, height) => {
  const canvas = canvasSelection.node();

  let scale = window.devicePixelRatio;
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  canvas.width = width * scale;
  canvas.height = height * scale;
  var context = canvas.getContext("2d");
  context.scale(scale, scale);
  return context;
};

export const useD3 = (renderChartFn, width, height, dependencies) => {
  const ref = useRef();

  useEffect(() => {
    const svg = d3.select(ref.current);
    console.log("BRUSH WIDTH HEIGHT")
    console.log(width)
    console.log(height)
    svg.style("width", width)
    svg.style("height",height)
    svg.style("position","absolute");
    renderChartFn(svg);
    return () => {};
  }, dependencies);
  
  return ref;
};

export const useCanvas = (renderCanvas, width, height, dependencies) => {
  const ref = useRef(null);

  useEffect(() => {
    const canvas = ref.current
    canvas.width = width 
    canvas.height = height + 100
    canvas.style.position = "absolute"
    canvas.style.marginLeft = 15

    renderCanvas(canvas);
  }, dependencies);

  return ref;
};

export const drawChromosomeAxis= (context, chromosomes, chromMap, categoryWidth, y) =>{
      //create ChromAxis for HeatMap
      chromosomes.forEach((chromosome,index) =>{
        context.beginPath();
        const key = chromosome.id;
        const chromosomeID = chromosome.id;
        const currData = chromMap[chromosomeID];
        //categoryWidth, index, 
        const x = currData["x"] ;
        const width = currData["width"]
        const height = heatmapConfig.chromosome["height"]
        const fill = heatmapConfig.chromosome["color"][index%2]
        //text position
        const text_width = currData["x"] + currData["width"] / 2 + categoryWidth
        context.fillStyle = fill;
        context.fillRect(x,y + 7,width,height)
        context.stroke()
        context.closePath();
        context.fill()
        context.fillStyle = "#000000"
        
        context.fillText(chromosomeID, text_width -32, y + 17);
  
      })


}
