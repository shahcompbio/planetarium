import React, { useRef, useEffect, useCallback, useState } from "react";


import * as d3 from "d3";
import { heatmapConfig } from "./config.js";
import {event} from "d3";
import {useD3} from "../../utils/useD3";
import {useCanvas} from "../../utils/useCanvas";

import {
  getMinimapBPRatio,
  getMinimapYScale,
  getSegX,
  getSegWidth,
  colorScale,
  heatmapToMinimapScale
} from "./utils.js";


const Minimap = ({
  heatmapConfig,
  heatmapOrder,
  rangeExtent,
  analysis,
  chromosomes,
  chromMap,
  triggerHeatmapRequery,
  segs
}) => {
  const cellCount = heatmapOrder.length;

  const [range, setRange] = useState([...rangeExtent]);

  const numRows = Math.min(
    cellCount,
    Math.floor(heatmapConfig.height / heatmapConfig.minimap.rowHeight)
  );
  const ratio = Math.ceil(cellCount / numRows);

  const [data, setData] = useState(segs);

  const bpRatio = getMinimapBPRatio(heatmapConfig,chromosomes);

  useEffect(() => {
    setRange([...rangeExtent]);
  }, [rangeExtent]);

  const canvasRef = useCanvas(
    (canvas)=> {
      canvas.style.position = "absolute"
      const context = canvas.getContext("2d");
      const yScale = getMinimapYScale(data.length);
      data.forEach((record, index) => {
        const y = yScale(index);
        record.segs.forEach(seg => {
          context.fillStyle = colorScale(seg.state);
          context.fillRect(
            getSegX(seg, chromMap, bpRatio, true),
            y,
            getSegWidth(seg, bpRatio),
            heatmapConfig.minimap.rowHeight
          );
        });
      });
  },
  heatmapConfig.minimap.width,
  heatmapConfig.height,
  []
  )

  const brushsvgRef = useD3(
    
    (svg)=>{
    const heatmapToMinimap = heatmapToMinimapScale(cellCount, data.length);
    const brushSvg = svg
    brushSvg.style("position","absolute");

    var brush = d3
      .brushY()
      .extent([
        [0, 0],
        [
          heatmapConfig.minimap.width - 10,
          data.length * heatmapConfig.minimap.rowHeight
        ]
      ])
      .on("end", brushEnd);
    function brushEnd() {
      const selection = d3.event.selection;
      if (!d3.event.sourceEvent || !selection) return;

      if (range[1] !== selection[1]) {
        let heatmapIndex = heatmapOrder.filter(
          (order, index) =>
            index >= Math.round(heatmapToMinimap.invert(selection[0])) &&
            index <= Math.round(heatmapToMinimap.invert(selection[1]))
        );
  
        if (heatmapIndex[heatmapIndex.length-1] === heatmapOrder[heatmapOrder.length-1] && heatmapIndex.length < range[1]-range[0]+1){
          const indexToAdd = heatmapOrder.indexOf(heatmapIndex[0])
          heatmapIndex.unshift(heatmapOrder[indexToAdd-1])
        }
        setRange([...selection]);

        triggerHeatmapRequery(heatmapIndex);
      }
    }
    brushSvg.selectAll(".brush").remove();
    var gBrush = brushSvg.append("g").attr("class", "brush");

    gBrush.call(brush);
    brushSvg.select(".handle").remove();
    brushSvg.select(".brush>.overlay").remove();
    brushSvg.select(".brush>.handle").remove();

    brush.move(gBrush, [
      heatmapToMinimap(range[0]),
      heatmapToMinimap(range[1])
    ]);
  
  },
  heatmapConfig.minimap.width,
  heatmapConfig.height,
  [])


  return (
    <div style={{position: "absolute", marginLeft: 15 }}>
      <canvas
        ref = {canvasRef}
      />
      <svg
        ref ={brushsvgRef}
      />

    </div>

  );
};


export default Minimap;
