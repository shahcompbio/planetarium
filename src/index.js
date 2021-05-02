import React from 'react'
import styles from './styles.module.css'

export const ExampleComponent = ({ text }) => {
  return <div className={styles.test}>Example Component: {text}</div>
}
//bar
export { default as ProbabilityHistogram } from './Bar/ProbabilityHistogram'
export { default as StackedHorizontalBar } from './Bar/StackedHorizontalBar'
//genomeProfile
export { default as GenomeProfile } from './GenomeProfile/GenomeProfile'
//copynumberheatmap

export { default as CopyNumberHeatmap } from './CopyNumberHeatmap/CopyNumberHeatmap'
export { default as CopyNumberHeatmapGenomeProfile } from './CopyNumberHeatmapGenomeProfile/CopyNumberHeatmapGenomeProfile'
//heatmap
export { default as Heatmap } from './Heatmap/Heatmap'
//infoBar
export { default as InfoBar } from './InfoBar/InfoBar'
export { default as Layout } from './InfoBar/Layout'

export { default as VerticalLegend } from './Legend/VerticalLegend'

export { useCanvas } from './utils/useCanvas'
export { useD3 } from './utils/useD3'

export { isValueHighlighted, isHighlighted } from './utils/isHighlighted'
