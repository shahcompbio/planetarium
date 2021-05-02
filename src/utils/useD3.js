import { useRef, useEffect } from "react";
import * as d3 from "d3";

// https://www.pluralsight.com/guides/using-d3.js-inside-a-react-app
export const useD3 = (renderChartFn, width, height, dependencies) => {
  const ref = useRef();

  useEffect(() => {
    const svg = d3.select(ref.current);
    svg.attr("width", width).attr("height", height);
    svg.selectAll("*").remove();
    renderChartFn(svg);
    return () => {};
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [width, height, ...dependencies]);
  return ref;
};
