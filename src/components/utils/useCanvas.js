import { useRef, useEffect } from "react";

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
  }, []);

  useEffect(() => {
    const canvas = ref.current;
    const context = canvas.getContext("2d");
    context.clearRect(0, 0, canvas.width, canvas.height);
    renderCanvas(canvas);
  }, dependencies);

  return ref;
};
