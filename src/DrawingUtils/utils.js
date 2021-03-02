import * as d3 from "d3";
export function canvasInit(canvasSelection, width, height) {
  const canvas = canvasSelection.node();

  let scale = window.devicePixelRatio;
  canvas.style.width = width + "px";
  canvas.style.height = height + "px";
  canvas.width = width * scale;
  canvas.height = height * scale;
  var context = canvas.getContext("2d");
  context.scale(scale, scale);
  return context;
}
export function drawAxisLabels(context, x, y, dim, xParam, yParam) {
  context.save();

  context.translate(dim.chart.x1 / 2, dim.chart.x2 / 2);
  context.rotate(-Math.PI / 2);

  context.fillText(yParam, 0, -10);

  context.restore();
  context.fillText(
    xParam,
    dim.chart.x2 / 2,
    dim.chart.y2 + dim.margin.bottom / 2 + 15
  );
  context.stroke();
  context.save();
}
export function drawAxis(context, x, y, dim) {
  context.beginPath();
  context.moveTo(dim.x1, dim.y1);
  context.lineTo(dim.x1, dim.y2);
  context.stroke();

  context.beginPath();
  context.moveTo(dim.x1, dim.y2);
  context.lineTo(dim.x2, dim.y2);
  context.stroke();
}
export function drawAxisTicks(context, x, y, dim) {
  const tickFormat = d3.format(".2s");
  x.ticks(4).forEach(function(d) {
    context.fillStyle = "#000000";
    context.fillText(
      d > 1000 ? tickFormat(d) : d,
      x(d),
      dim.y2 + dim.margin.bottom / 2
    );
  });

  y.ticks(4).forEach(function(d) {
    context.fillStyle = "#000000";
    context.fillText(d > 1000 ? tickFormat(d) : d, dim.x1 - 30, y(d));
  });
}
