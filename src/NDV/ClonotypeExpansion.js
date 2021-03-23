import React from "react";

import StackedHorizontalBar from "../components/Bar/Barplot";
import { useDashboardState } from "../PlotState/dashboardState";
import _ from "lodash";

const ClonotypeExpansion = ({ data, dim, chartName }) => {
  const [{ clonotypeParam, subtypeParam }] = useDashboardState();

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
    <StackedHorizontalBar
      data={countedClonotypes}
      chartDim={dim}
      barLabels={Array.from(Array(10).keys()).map((value) => ({
        value: value + 1,
        label: value === 9 ? "≥10" : value + 1,
      }))}
      chartName={chartName}
    />
  );
};

export default ClonotypeExpansion;
