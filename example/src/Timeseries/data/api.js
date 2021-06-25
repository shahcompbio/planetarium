import { useEffect, useState } from "react";

import * as d3 from "d3";

import metadataSource from "./vdjTimeSeries.tsv";

import genesSource from "./genes.tsv";

const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    Promise.all([d3.tsv(metadataSource), d3.tsv(genesSource)]).then((data) => {
      setData({
        metadata: data[0],
        genes: data[1],
      });
    });
  }, []);

  return data;
};

export default useFetchData;
