import PackingCirclesComponent from "./PackingCircles";
import data from "./data/cellmine.json";

const Template = (args) => <PackingCirclesComponent {...args} />;

export default {
  title: "Components/PackingCircles/PackingCircles",
  component: PackingCirclesComponent,
};

export const Cellmine = Template.bind({});
Cellmine.args = {
  height: 800,
  width: 950,
  data,
  radiusParam: "num_sublibraries",
  idParam: "jira_ticket",
  tooltipFields: [
    { label: "Analysis Ticket", param: "jira_ticket" },
    { label: "Cell Count", param: "num_sublibraries" },
    { label: "Description", param: "description" },
  ],
};

export const HighlightedCellmine = Template.bind({});
HighlightedCellmine.args = {
  height: 800,
  width: 950,
  data,
  radiusParam: "num_sublibraries",
  idParam: "jira_ticket",
  highlightedIDs: ["SC-905", "SC-881"],
  tooltipFields: [
    { label: "Analysis Ticket", param: "jira_ticket" },
    { label: "Cell Count", param: "num_sublibraries" },
    { label: "Description", param: "description" },
  ],
};
