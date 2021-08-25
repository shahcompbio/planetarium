export const loadFont = (fontname) => {
  var canvas = document.createElement("canvas");
  //Setting the height and width is not really required
  canvas.width = 16;
  canvas.height = 16;
  var ctx = canvas.getContext("2d");

  //There is no need to attach the canvas anywhere,
  //calling fillText is enough to make the browser load the active font

  //If you have more than one custom font, you can just draw all of them here
  ctx.fillStyle = "none";
  ctx.font = "4px " + fontname;
  ctx.fillText("text", 0, 8);
};
