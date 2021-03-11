const chartHeight = 550;
const chartWidth = 550;

let margin = { top: 100, right: 10, bottom: 50, left: 20 };
let chartDim = {
  margin: margin,
  width: 650,
  height: 700,
  chart: {
    x1: margin.left,
    y1: margin.top,
    x2: margin.left + chartWidth,
    y2: margin.top + chartHeight
  },
  legend: {
    x1: margin.left,
    y1: margin.top,
    x2: margin.left + chartWidth,
    y2: margin.top + chartHeight
  }
};
export default chartDim;
