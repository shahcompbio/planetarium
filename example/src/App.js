import React from 'react'

import { GenomeProfile } from 'shahlab-planetarium'

const genomeProfileData = require("./genomeProfileData.json")
const gpbpTotal = genomeProfileData.bptotal
const gpchromosomes = genomeProfileData.chromosomes
const gpmaxState = genomeProfileData.maxState
const gpbins = genomeProfileData.bins
const gpsegs = genomeProfileData.segs

const App = () => {
  return         <GenomeProfile width = {450} height= {350} bins = {gpbins} segs = {gpsegs} bpTotal = {gpbpTotal} chromosomes = {gpchromosomes} maxState = {gpmaxState}/>
}

export default App
