import { useEffect, useState } from "react";

import * as d3 from "d3";

import entropy from "./clone_entropy.csv";

const GetCloneApi = () => {
  const [data, setData] = useState({});

  /*useEffect(() => {
    d3.csv(entropy).then((data) => {
      console.log(data);
      setData({
        data: data
          .filter((d) => d["Treatment"] === "Post Treatment")
          .sort(function (a, b) {
            return parseFloat(a["Exhaustion"]) - parseFloat(b["Exhaustion"]);
          }),
      });
    });
  }, []);*/

  return data;
};

export default GetCloneApi;
