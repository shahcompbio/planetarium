import { useEffect, useState } from "react";

//const url = process.env.HOST ? process.env.HOST : "http://127.0.0.1:5000";
//const url = "https://spectrum-staging.shahlab.mskcc.org";
const useFetchData = () => {
  const [data, setData] = useState({});

  useEffect(() => {
    fetch("http://127.0.0.1:5000/getSantoshData/", {
      credentials: "include",
    })
      .then((res) => res.json())
      .then((data) => {
        setData({ metadata: data["metadata"], filters: data["filters"] });
      });
  }, []);

  return data;
};

export default useFetchData;
