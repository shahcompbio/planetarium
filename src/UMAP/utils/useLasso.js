import { useState } from "react";
import * as d3 from "d3";

const LASSO_COLOR = "#c0392b";

const useLasso = (data, xScale, yScale, xParam, yParam) => {
  const [lassoPolys, setLassoPolys] = useState([]);
  const [highlighted, setHighlighted] = useState(null);

  const getHighlighted = (polys) => {
    if (polys.length > 0) {
      return data.filter((datum) =>
        isPointInPoly(polys, xScale(datum[xParam]), yScale(datum[yParam]))
      );
    }
    return null;
  };

  const addLassoHandler = (canvas, disable, onLasso) => {
    const context = canvas.getContext("2d");

    let polys = [];
    d3.select(canvas)
      .on("mousemove", (d, i, e) => {
        if (disable) {
          return;
        }
        if (d3.event.buttons === 1) {
          const poly = d3.mouse(e[0]);
          polys.push(poly);
          context.lineTo(poly[0], poly[1]);
          context.stroke();
        }
      })
      .on("mousedown", function mousedown() {
        if (disable) {
          return;
        }
        polys = [];
        context.lineCap = "round";
        context.strokeStyle = LASSO_COLOR;
        context.lineWidth = 3;

        context.restore();
        context.beginPath();
      })
      .on("mouseup", () => {
        if (disable) {
          return;
        }
        setLassoPolys(polys);
        const lassoedData = getHighlighted(polys);
        setHighlighted(lassoedData);
        onLasso(lassoedData);
      });
  };

  const drawLasso = (context) => {
    context.globalAlpha = 0.2;
    context.fillStyle = "black";
    context.fill();

    const lassoPath = lassoPolys
      .map((poly, idx) => `${idx === 0 ? "M" : " L"} ${poly[0]} ${poly[1]}`)
      .reduce((str, poly) => str + poly, "");

    const lassoSvgPath = new Path2D(lassoPath);
    lassoSvgPath.closePath();

    context.fill(lassoSvgPath, "evenodd");
    context.closePath();
  };

  const resetLasso = () => {
    if (lassoPolys.length > 0) {
      setLassoPolys([]);
      setHighlighted(null);
    }
  };

  // should just returned lasso'd data
  // with canvas, you might need to return the exact handlers, which is annoying

  return [highlighted, drawLasso, addLassoHandler, resetLasso];
};

function isPointInPoly(poly, x, y) {
  for (var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
    ((poly[i][1] <= y && y < poly[j][1]) ||
      (poly[j][1] <= y && y < poly[i][1])) &&
      x <
        ((poly[j][0] - poly[i][0]) * (y - poly[i][1])) /
          (poly[j][1] - poly[i][1]) +
          poly[i][0] &&
      (c = !c);
  return c;
}

export default useLasso;
