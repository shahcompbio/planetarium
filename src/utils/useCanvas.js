import { useRef, useEffect, useState } from "react";

//var FontFaceObserver = require("fontfaceobserver");

/*var bold = new FontFaceObserver("MyFontBold");
var regular = new FontFaceObserver("MyFontRegular");
var light = new FontFaceObserver("MyFontLight");
*/
export const useCanvas = (renderCanvas, width, height, dependencies) => {
  const ref = useRef(null);

  useEffect(() => {
    const canvas = ref.current;
    const context = canvas.getContext("2d");

    let scale = window.devicePixelRatio;
    canvas.style.width = width + "px";
    canvas.style.height = height + "px";
    canvas.width = width * scale;
    canvas.height = height * scale;

    context.scale(scale, scale);
  }, [width, height]);

  useEffect(() => {
    const canvas = ref.current;
    const context = canvas.getContext("2d");

    context.clearRect(0, 0, canvas.width, canvas.height);
    renderCanvas(canvas);

    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [...dependencies]);

  return ref;
};
