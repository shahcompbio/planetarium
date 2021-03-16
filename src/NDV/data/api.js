import { useEffect, useState } from "react";

import * as d3 from "d3";

import metadataSource from "./metadata.tsv";
import probabilitiesSource from "./probabilities.tsv";
import degsSource from "./degs.tsv";

const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    Promise.all([
      d3.tsv(metadataSource),
      d3.tsv(probabilitiesSource),
      d3.tsv(degsSource)
    ]).then(data => {
      setData({ metadata: data[0], probabilities: data[1], degs: data[2] });
    });
  }, []);

  return data;
};

export default useFetchData;