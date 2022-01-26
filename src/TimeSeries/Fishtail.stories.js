import FishtailComponent from "./Fishtail";
import data from "./data/timeseries.json";
import sankey from "./data/sankey_filtered2.json";

const Template = (args) => <FishtailComponent {...args} />;

export default {
  title: "Components/TimeSeries/Fishtail",
  component: FishtailComponent,
};

export const Default = Template.bind({});
Default.args = { width: 700, height: 400, data, subsetParam: "clone" };

export const HighlightSubset = Template.bind({});
HighlightSubset.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  subset: "A",
};

export const HighlightTimepoint = Template.bind({});
HighlightTimepoint.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  timepoint: "X8",
};

export const Static = Template.bind({});
Static.args = {
  width: 700,
  height: 400,
  data,
  subsetParam: "clone",
  disable: true,
};
export const TwoTimePoints = Template.bind({});
TwoTimePoints.args = {
  width: 700,
  height: 400,
  data: sankey,
  subsetParam: "cell_type",
  timepointParam: "treatment",
  timepointOrder: ["Pre", "Post"],
  disable: true,
};

const TWO_CLONE = [
  /*  {
    cell_id: "SA535X6XB03099_AAACCCAGTTTCACTT-1",
    clone: "A",
    timepoint: "X0.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAAAGAGCATTA-1",
    clone: "A",
    timepoint: "X0.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAAGTGCCTTTC-1",
    clone: "A",
    timepoint: "X0.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X0.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X0.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
*/
  {
    cell_id: "SA535X6XB03099_AAACCCAGTTTCACTT-1",
    clone: "A",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAAAGAGCATTA-1",
    clone: "A",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAAGTGCCTTTC-1",
    clone: "A",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "A",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "A",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Memory",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X1",
    therapy: "Untreated",
    subtype: "Memory",
  },

  /*  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "A",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "A",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Activated",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAATCTTACACT-1",
    clone: "B",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Memory",
  },
  {
    cell_id: "SA535X6XB03099_AAACGCTTCGAGTGGA-1",
    clone: "B",
    timepoint: "X1.5",
    therapy: "Untreated",
    subtype: "Memory",
  },*/
];
export const TwoClone = Template.bind({});
TwoClone.args = {
  width: 700,
  height: 400,
  data: TWO_CLONE,
  subsetParam: "clone",
  timepointParam: "timepoint",
  timepointOrder: ["X0", "X1"],
  addTwoTimepointCurve: true,
  disable: true,
};
