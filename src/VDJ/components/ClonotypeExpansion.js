import React from "react";

import StackedHorizontalBar from "../../components/Bar/StackedHorizontalBar";
import Layout from "../../components/InfoBar/Layout";
import infoText from "../InfoText";
import { CONSTANTS } from "../config";
import _ from "lodash";

const ClonotypeExpansion = ({
  data,
  width,
  height,
  chartName,
  highlightedRow,
}) => {
  const { clonotypeParam, subtypeParam } = CONSTANTS;

  const groupedSubtype = _.groupBy(data, subtypeParam);
  const subtypes = Object.keys(groupedSubtype).sort();

  const countedClonotypes = subtypes.reduce((countMap, subtype) => {
    const clonotypeCount = _.countBy(groupedSubtype[subtype], clonotypeParam);

    const countFreq = _.countBy(Object.values(clonotypeCount), (value) =>
      Math.min(value, 10)
    );

    return { ...countMap, [subtype]: countFreq };
  }, {});

  return (
    <Layout
      title={infoText["BARPLOT"]["title"]}
      infoText={infoText["BARPLOT"]["text"]}
    >
      <StackedHorizontalBar
        highlightedRow={highlightedRow}
        data={countedClonotypes}
        width={width}
        height={height}
        barLabels={Array.from(Array(10).keys()).map((value) => ({
          value: value + 1,
          label: value === 9 ? "≥10" : value + 1,
        }))}
        chartName={chartName}
      />
    </Layout>
  );
};

export default ClonotypeExpansion;
