const AXIS_FONT = "normal 10px Helvetica";
const AXIS_COLOR = "#000000";
const AXIS_LENGTH = 50;

export default (context, xPos, yPos, xParam, yParam) => {
  context.beginPath();
  context.font = AXIS_FONT;
  context.globalAlpha = 1;

  context.fillStyle = AXIS_COLOR;
  context.strokeStyle = AXIS_COLOR;
  context.lineWidth = 1;
  context.lineCap = "butt";
  context.moveTo(xPos, yPos);
  context.lineTo(xPos, yPos - AXIS_LENGTH);
  context.stroke();

  context.beginPath();
  context.moveTo(xPos, yPos);
  context.lineTo(xPos + AXIS_LENGTH, yPos);
  context.stroke();

  context.textAlign = "left";
  context.textBaseline = "middle";
  context.fillText(xParam, xPos + AXIS_LENGTH + 2, yPos);
  context.save();
  context.rotate((270 * Math.PI) / 180);
  context.fillText(yParam, -(yPos - AXIS_LENGTH - 2), xPos);
  context.restore();
};
