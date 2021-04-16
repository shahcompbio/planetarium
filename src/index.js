import React from 'react'
import styles from './styles.module.css'

export const ExampleComponent = ({ text }) => {
  return <div className={styles.test}>Example Component: {text}</div>
}
//bar
export {default as ProbabilityHistogram } from "./components/Bar/ProbabilityHistogram";
export {default as StackedHorizontalBar } from "./components/Bar/StackedHorizontalBar";
//genomeProfile
export {default as GenomeProfile } from "./components/GenomeProfile/GenomeProfile";
//heatmap
export {default as Heatmap } from "./components/Heatmap/Heatmap";
//infoBar
export {default as InfoBar } from "./components/InfoBar/InfoBar";
export {default as Layout } from "./components/InfoBar/Layout";
//VDJ
export {default as ClonotypeExpansion } from "./VDJ/components/ClonotypeExpansion";
export {default as DEGTable } from "./VDJ/components/DEGTable";
export {default as subTypeUmap } from "./VDJ/components/subTypeUmap";
export {default as Umap } from "./VDJ/components/Umap";

