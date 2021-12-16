import all from "./test.json";

const fetchData = () => ({
  patients: Object.keys(all[0]),
  geneData: all[0],
});

export default fetchData;
