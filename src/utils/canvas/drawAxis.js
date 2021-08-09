const LABEL_FONT = "normal 14px Helvetica";
const TICK_FONT = "normal 12px Helvetica";

const drawAxis = ({
  context,
  xScale,
  yScale,
  label = "",
  ticks = 2,
  orientation = "vertical",
  format = (tick) => tick,
  gridlines = true,
}) => {
  context.beginPath();

  context.globalAlpha = 1;
  context.fillStyle = "black";
  context.textAlign = "right";

  context.font = TICK_FONT;
  context.textBaseline = "bottom";

  const scale = orientation === "vertical" ? yScale : xScale;
  const [minValue, maxValue] = scale.domain();
  const minPos = scale(minValue);
  const maxPos = scale(maxValue);

  const lengthScale = orientation === "vertical" ? xScale : yScale;
  const length0 = lengthScale(lengthScale.domain()[0]);
  const length1 = lengthScale(lengthScale.domain()[1]);

  scale.ticks(ticks).forEach((tick) => {
    // Tick text
    context.globalAlpha = 1;
    context.textBaseline = "middle";
    if (orientation === "vertical") {
      context.fillText(format(tick), length0 - 8, scale(tick));
    } else {
      context.textAlign = "center";
      context.fillText(tick, scale(tick), Math.max(length0, length1) + 8);
    }

    if (gridlines) {
      // Tick line
      context.globalAlpha = 0.2;
      context.lineWidth = 0.5;
      context.beginPath();
      if (orientation === "vertical") {
        context.moveTo(length0, scale(tick));
        context.lineTo(length1, scale(tick));
      } else {
        context.moveTo(scale(tick), length0);
        context.lineTo(scale(tick), length1);
      }
      context.stroke();
    }
  });

  // Label
  context.globalAlpha = 1;
  context.font = LABEL_FONT;
  context.textAlign = "center";
  context.textBaseline = "hanging";

  const midPos = (maxPos + minPos) / 2;
  if (orientation === "vertical") {
    context.save();
    context.rotate((270 * Math.PI) / 180);
    context.fillText(label, -midPos, length0 - 50);
    context.restore();
  } else {
    context.fillText(label, midPos, Math.max(length0, length1) + 20);
  }
};

export default drawAxis;
