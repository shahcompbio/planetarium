import SankeyComponent from "./Sankey";

const Template = (args) => <SankeyComponent {...args} />;

export default {
  title: "Components/TimeSeries/Sankey",
  component: SankeyComponent,
};

const ONE_CLONE = [
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
];

export const Default = Template.bind({});
Default.args = {
  width: 700,
  height: 400,
  data: ONE_CLONE,
  subsetParam: "subtype",
  cloneParam: "clone",
};

const TWO_CLONE = [
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
];

export const Example2 = Template.bind({});
Example2.args = {
  width: 700,
  height: 400,
  data: TWO_CLONE,
  subsetParam: "subtype",
  cloneParam: "clone",
};

const EXAMPLE_3 = [
  {
    cell_id: "SA535X6XB03099_AAACCCAGTTTCACTT-1",
    clone: "A",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Dysfunctional",
  },
  {
    cell_id: "SA535X6XB03099_AAACGAAAGAGCATTA-1",
    clone: "A",
    timepoint: "X0",
    therapy: "Untreated",
    subtype: "Dysfunctional",
  },
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
];

export const Example3 = Template.bind({});
Example3.args = {
  width: 700,
  height: 400,
  data: EXAMPLE_3,
  subsetParam: "subtype",
  cloneParam: "clone",
};
