import React, { useState, useEffect, useCallback, useRef } from "react";
import Grid from "@material-ui/core/Grid";
import CopyNumberHeatmap from  "./CopyNumberHeatmap/CopyNumberHeatmap";
import GenomeProfile from "./GenomeProfile/GenomeProfile";
import reactDomTestUtilsProductionMin from "react-dom/cjs/react-dom-test-utils.production.min";
const genomeProfileData = require("../../test_data/GenomeProfile/genomeProfileData.json")
const gpchromosomes = genomeProfileData.chromosomes
const gpbins = genomeProfileData.bins
const gpsegs = genomeProfileData.segs


const getBpTotal = (chromosomes)=>{
  return chromosomes.reduce(function(a,b){
    return a + b["end"] 
  }, 0)

}
const filterSegs = (segs,hoverCell) =>{
  if (Object.keys(hoverCell["cell"]).length !== 0){
    return segs.find(seg => seg.id === hoverCell["cell"]["id"])
  } 
  return []

};

const filterBins = (bins,hoverCell) =>{
  if (Object.keys(hoverCell["cell"]).length !== 0){
    return bins.filter(function(bin){
      return bin["cell_id"] === hoverCell["cell"]["id"]
    })
  } 
  return []

}
const CopyNumberHeatmapGenomeProfile  = ({width,height, analysis, allHeatmapOrder, chromosomes, categoryStats, segs, analysisStats, bins, maxState }) =>{

  const [hoverCell, setHoverCell] = useState({ cell: {} });
  const setHoverCellCoordinate =
      (hoverCell) => {
            setHoverCell({
              ...hoverCell
            }); 
      }
    const genomeProfileSegs = filterSegs(segs,hoverCell)
    const genomeProfileBins = filterBins(bins,hoverCell)
    /*set genomeprofile necessary bins and segs based on hovercell here*/
    //<GenomeProfile width = {450} height= {350} bins = {gpbins} segs = {gpsegs} bpTotal = {gpbpTotal} chromosomes = {gpchromosomes} maxState = {gpmaxState}/>
    //console.log(gpbpTotal,gpchromosomes,gpmaxState,gpbins,gpsegs)

    //we change bins and segs
    const bpTotal = getBpTotal(chromosomes)

    //we can just inherit chromosomes, maxstate , not sure what bp total is
    return(      
        <Grid container direction ="column">
          <Grid item>
            <CopyNumberHeatmap width ={width} height = {height} setHoverCellCoordinate ={setHoverCellCoordinate} analysis = {analysis} allHeatmapOrder = {allHeatmapOrder} chromosomes = {chromosomes} categoryStats = {categoryStats} segs = {segs} analysisStats ={analysisStats} />
          </Grid>
          <Grid item style={{marginTop:125}}>
          <GenomeProfile width = {550} height= {350} bins = {gpbins} segs = {gpsegs} bpTotal = {bpTotal} chromosomes = {gpchromosomes} maxState = {analysisStats.maxState}/>
          </Grid>
          <Grid item style = {{marginTop:50}}>
            <div id="heatmapCellID" style={{ height: 15, margin: 5 }}>
                {Object.keys(hoverCell["cell"]).length !== 0 && (
                    <div>Cell ID: {hoverCell["cell"]["id"]}</div>
                )}
              </div>
          </Grid>

    </Grid>
  )
}

export default CopyNumberHeatmapGenomeProfile;