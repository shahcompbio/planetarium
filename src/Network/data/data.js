import { useEffect, useState } from "react";

import * as d3 from "d3";

import dataSource from "./interactions.csv";

const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    d3.csv(dataSource).then((data) => {
      setData([...data]);
    });
  }, []);

  return data;
};

export default useFetchData;
