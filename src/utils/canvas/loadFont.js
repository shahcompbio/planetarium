export const loadFont = (fontname) => {
  const canvas = document.createElement("canvas");

  const ctx = canvas.getContext("2d");

  ctx.fillStyle = "none";
  ctx.font = "1px " + fontname;
  ctx.fillText("text", 0, 8);
};
