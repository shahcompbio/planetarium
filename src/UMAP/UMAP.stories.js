import React, { useState } from "react";
import _ from "lodash";

import UMAP from "./UMAP";
import numData from "./data/cells.json";
import catData from "./data/metadata.json";

const NumTemplate = (args) => <UMAP data={numData} {...args} />;
const CatTemplate = (args) => <UMAP data={catData} {...args} />;
const MoreInfoComp = () => <div>More info</div>;

export default {
  title: "Components/UMAP/UMAP",
  component: UMAP,
};

export const Numerical = NumTemplate.bind({});
Numerical.args = {
  width: 800,
  height: 600,
  xParam: "rna_UMAP_1",
  yParam: "rna_UMAP_2",
  subsetParam: "VIM",
  idParam: "cell_id",
};

export const Categorical = CatTemplate.bind({});
Categorical.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
  idParam: "cell_id",
};

export const CategoricalLabels = CatTemplate.bind({});
CategoricalLabels.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
  idParam: "cell_id",
  labels: (value) => `P ${value}`,
};

export const CategoricalSubset = CatTemplate.bind({});
CategoricalSubset.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  idParam: "cell_id",
  subsetParam: "subtype",
  highlightIDs: catData
    .filter((datum) => datum["subtype"] === "Activated CD8")
    .map((datum) => datum["cell_id"]),
};

export const Static = CatTemplate.bind({});
Static.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  idParam: "cell_id",
  subsetParam: "subtype",
  disable: true,
};

// Test for external reset component
const CatButtonTemplate = (args) => {
  const [highlightData, setHighlightData] = useState(null);
  return (
    <div>
      <UMAP
        data={catData}
        {...args}
        highlightIDs={
          highlightData === null
            ? highlightData
            : highlightData.map((datum) => datum["cell_id"])
        }
        onLasso={setHighlightData}
      />
      <button onClick={() => setHighlightData(null)}>Reset Lasso</button>
    </div>
  );
};
export const CategoricalReset = CatButtonTemplate.bind({});
CategoricalReset.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
  idParam: "cell_id",
};

export const MoreInfo = CatButtonTemplate.bind({});
MoreInfo.args = {
  width: 800,
  height: 600,
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  subsetParam: "subtype",
  idParam: "cell_id",
  MoreInfoComponent: MoreInfoComp,
};
