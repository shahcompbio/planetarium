import React from 'react'

import { GenomeProfile,CopyNumberHeatmap,CopyNumberHeatmapGenomeProfile } from 'shahlab-planetarium'

const genomeProfileData = require("./testData/genomeProfileData.json");
const gpbpTotal = genomeProfileData.bptotal;
const gpchromosomes = genomeProfileData.chromosomes;
const gpmaxState = genomeProfileData.maxState;
const gpbins = genomeProfileData.bins;
const gpsegs = genomeProfileData.segs;

const copyNumberHeatmapData = require("./testData/11655_test_data.json");


const analysis =  copyNumberHeatmapData.analysis;
const allHeatmapOrder = copyNumberHeatmapData.heatmapOrder;
const categoryStats= copyNumberHeatmapData.categoriesStats;
const chromosomes = copyNumberHeatmapData.chromosomes;
const segs = copyNumberHeatmapData.segs;
const analysisStats = copyNumberHeatmapData.analysisStats;

const App = () => {
  return   <div>
    <CopyNumberHeatmap width ={750} height = {539} analysis = {analysis} allHeatmapOrder = {allHeatmapOrder} chromosomes = {chromosomes} categoryStats = {categoryStats} segs = {segs} analysisStats ={analysisStats} />
    <GenomeProfile width = {450} height= {350} bins = {gpbins} segs = {gpsegs} bpTotal = {gpbpTotal} chromosomes = {gpchromosomes} maxState = {gpmaxState}/>

  </div>


}

export default App
