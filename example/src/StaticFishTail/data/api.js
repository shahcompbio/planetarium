import all from "./test1.json";

const fetchData = () => ({
  clones: all["clones"],
  data: all["data"],
  //  data: all["data"].filter((d) => d["clone"].indexOf("MK") !== -1),
});

export default fetchData;
