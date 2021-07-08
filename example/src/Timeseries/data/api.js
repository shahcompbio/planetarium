import { useEffect, useState } from "react";

import * as d3 from "d3";

import metadataSource from "./vdjTimeSeries.tsv";

const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    Promise.all([d3.tsv(metadataSource)]).then((data) => {
      setData({
        metadata: data[0],
      });
    });
  }, []);

  return data;
};

export default useFetchData;
