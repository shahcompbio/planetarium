import { useEffect, useState } from "react";
import * as d3 from "d3";
//import all from "./newprotien.tsv";
import data from "./newpro.json";
//import data from "./genes.json";
const pedsFetchData = () => {
  //  const d = await d3.tsv(all);
  //  console.log(d);
  return [...data];
};

export default pedsFetchData;
