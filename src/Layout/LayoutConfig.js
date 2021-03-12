const chartHeight = 500;
const chartWidth = 500;

let margin = { top: 100, right: 10, bottom: 50, left: 20 };
const chartDim = {
  margin: margin,
  width: 550,
  height: 700,
  chart: {
    x1: margin.left,
    y1: margin.top,
    x2: margin.left + chartWidth,
    y2: margin.top + chartHeight
  },
  legend: {
    x1: 10,
    y1: 20,
    x2: 70,
    y2: margin.top + chartHeight
  }
};
export default chartDim;
