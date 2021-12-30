import React, { useState } from "react";
import * as d3 from "d3";
import Grid from "@material-ui/core/Grid";
import TextField from "@material-ui/core/TextField";
import Checkbox from "@material-ui/core/Checkbox";

import _ from "lodash";
import { useCanvas } from "@shahlab/planetarium";
import CheckMarkSelect from "./CheckMarkSelect";

const Heatmap = ({ patients, data }) => {
  const [context, saveContext] = useState(null);

  const [settings, setSettings] = useState(
    localStorage.getItem("heatmapSettings") === null
      ? {
          paddingInner: 0.13,
          paddingOuter: 0.25,
          align: 0.5,
          columnPadding: 0.2,
          movePatientsPadding: 0,

          squareWidth: 25,
          squareHeight: 25,
          radius: 0,

          labelPadding: 38,
          labelVertPadding: 9,
          pValPadding: 0,

          topHeadingHeight: 100,

          patientSpacing: 25,
          patientWidth: 200,
          patientHeight: 700,
          blackOutline: false,
          allSectionPadding: 0,
          allLegendPadding: 0,
          allLabelPadding: 0,
          patientNamePadding: 0,
          switchAll: true,

          removePval: false,
          gradWidth: 170,
          gradHeight: 15,

          preLabel: "Pre",
          postLabel: "Post",
          pLabel: "P Val",

          colourScheme: true,
          patientLegendPadding: 0,
        }
      : JSON.parse(localStorage.getItem("heatmapSettings"))
  );
  const {
    radius,
    paddingInner,
    paddingOuter,
    align,
    squareWidth,
    squareHeight,
    labelPadding,
    labelVertPadding,
    columnPadding,

    patientHeight,
    patientWidth,
    patientSpacing,
    topHeadingHeight,
    pValPadding,
    patientLegendPadding,
    allSectionPadding,
    allLegendPadding,
    allLabelPadding,
    patientNamePadding,

    removePval,
    gradWidth,
    gradHeight,

    preLabel,
    postLabel,
    pLabel,
    movePatientsPadding,
  } = settings;

  const canvasHeight = patientHeight;
  const canvasWidth =
    (patients.length + 1) * patientWidth +
    (patients.length + 1) * patientSpacing;
  const [shownGenes, setShownGenes] = useState(Object.keys(data[patients[0]]));

  const patientScale = d3
    .scaleBand()
    .domain(patients)
    .range([40, canvasWidth - 40])
    .paddingInner(paddingInner)
    .paddingOuter(paddingOuter)
    .align(align)
    .round(false);

  const geneYScale = d3
    .scaleBand()
    .domain([...shownGenes])
    .range([topHeadingHeight, patientHeight - 70]);

  const canvasRef = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      saveContext(context);
      drawHeatmap(context, patients, data, patientScale, geneYScale, settings);

      drawLabels(context, shownGenes, geneYScale, settings);
    },
    canvasWidth,
    canvasHeight,
    [settings]
  );
  /*{" "}
  <RangeOption
    label={"Total Width"}
    value={patientWidth}
    setValue={(value) =>
      setSettings({ ...settings, patientWidth: value })
    }
    key={"rangeWidth"}
    min={"100"}
    max={canvasWidth}
    step={"10"}
  />
  */
  return (
    <div>
      <Grid container spacing={1} style={{ width: 700 }}>
        <Grid container item spacing={3}>
          <CheckMarkSelect
            allGenes={Object.keys(data[patients[0]])}
            selectedGenes={shownGenes}
            setSelected={setShownGenes}
          />
        </Grid>
        <Grid container item spacing={3}>
          <RangeOption
            label={"Inner Padding"}
            value={paddingInner}
            setValue={(value) =>
              setSettings({ ...settings, paddingInner: value })
            }
            min={"0"}
            max={"1"}
            step={"0.1"}
          />
          <RangeOption
            label={"Outer Padding"}
            value={paddingOuter}
            setValue={(value) =>
              setSettings({ ...settings, paddingOuter: value })
            }
            key={"range3"}
            min={"0"}
            max={"1"}
            step={"0.1"}
          />
          <RangeOption
            label={"Align"}
            value={align}
            setValue={(value) => setSettings({ ...settings, align: value })}
            key={"range4"}
            min={"0"}
            max={"1"}
            step={"0.1"}
          />
          <RangeOption
            label={"Column Padding"}
            value={columnPadding}
            setValue={(value) =>
              setSettings({ ...settings, columnPadding: value })
            }
            key={"columnPadding"}
            min={"0"}
            max={"1"}
            step={"0.1"}
          />
        </Grid>
        <Grid container item spacing={3}>
          <RangeOption
            label={"Total Height"}
            value={patientHeight}
            setValue={(value) =>
              setSettings({ ...settings, patientHeight: value })
            }
            key={"totalH"}
            min={"100"}
            max={"2000"}
            step={"10"}
          />
        </Grid>
        <Grid container item spacing={3}>
          <RangeOption
            label={"Square Width"}
            value={squareWidth}
            setValue={(value) =>
              setSettings({ ...settings, squareWidth: value })
            }
            key={"range5"}
            min={"2"}
            max={"60"}
            step={"1"}
          />
          <RangeOption
            label={"Square Height"}
            value={squareHeight}
            setValue={(value) =>
              setSettings({ ...settings, squareHeight: value })
            }
            key={"range6"}
            min={"2"}
            max={"60"}
            step={"1"}
          />
          <RangeOption
            label={"Square Radius"}
            value={radius}
            setValue={(value) => setSettings({ ...settings, radius: value })}
            key={"patientRad"}
            min={"0"}
            max={"20"}
            step={"0.1"}
          />
        </Grid>
        <Grid container item spacing={3}>
          <RangeOption
            label={"Patient Width"}
            value={patientWidth}
            setValue={(value) =>
              setSettings({ ...settings, patientWidth: value })
            }
            key={"patientWidth"}
            min={"50"}
            max={"300"}
            step={"10"}
          />
          <RangeOption
            label={"Patient Spacing"}
            value={patientSpacing}
            setValue={(value) =>
              setSettings({ ...settings, patientSpacing: value })
            }
            key={"range8"}
            min={"1"}
            max={"100"}
            step={"1"}
          />
        </Grid>

        <RangeOption
          label={"Header Section"}
          value={topHeadingHeight}
          setValue={(value) =>
            setSettings({ ...settings, topHeadingHeight: value })
          }
          key={"topheading"}
          min={"50"}
          max={"300"}
          step={"1"}
        />
      </Grid>

      <canvas ref={canvasRef} id="heatmapCanvas" />
      <Grid>
        <Grid container item spacing={3}>
          <RangeOption
            label={"Move patient padding"}
            value={movePatientsPadding}
            setValue={(value) =>
              setSettings({ ...settings, movePatientsPadding: value })
            }
            key={"movePatientsPadding"}
            min={"0"}
            max={canvasWidth / 3}
            step={"1"}
          />
          <RangeOption
            label={"All Section Padding"}
            value={allSectionPadding}
            setValue={(value) =>
              setSettings({ ...settings, allSectionPadding: value })
            }
            key={"allSectionPadding"}
            min={"0"}
            max={patientHeight}
            step={"1"}
          />
          <RangeOption
            label={"All Legend Padding"}
            value={allLegendPadding}
            setValue={(value) =>
              setSettings({ ...settings, allLegendPadding: value })
            }
            key={"allLegendPadding"}
            min={"0"}
            max={patientWidth}
            step={"1"}
          />
          <RangeOption
            label={"All Label Padding"}
            value={allLabelPadding}
            setValue={(value) =>
              setSettings({ ...settings, allLabelPadding: value })
            }
            key={"allLabelPadding"}
            min={"0"}
            max={patientWidth}
            step={"1"}
          />
        </Grid>
        <Grid container item spacing={3}>
          <p>Labels</p>
          <RangeOption
            label={"Label Vertical Padding"}
            value={labelVertPadding}
            setValue={(value) =>
              setSettings({ ...settings, labelVertPadding: value })
            }
            key={"vertPad"}
            min={"0"}
            max={squareHeight}
            step={"0.5"}
          />
          <RangeOption
            label={"Label Padding"}
            value={labelPadding}
            setValue={(value) =>
              setSettings({ ...settings, labelPadding: value })
            }
            key={"labelpadding"}
            min={"2"}
            max={"60"}
            step={"1"}
          />
          <RangeOption
            label={"Pval Padding"}
            value={pValPadding}
            setValue={(value) =>
              setSettings({ ...settings, pValPadding: value })
            }
            key={"pvalPadding"}
            min={"0"}
            max={"60"}
            step={"1"}
          />
          <RangeOption
            label={"Legend Padding"}
            value={patientLegendPadding}
            setValue={(value) =>
              setSettings({ ...settings, patientLegendPadding: value })
            }
            key={"patientLegendPadding"}
            min={"0"}
            max={canvasWidth - patientWidth}
            step={"1"}
          />
          <RangeOption
            label={"Patient name Padding"}
            value={patientNamePadding}
            setValue={(value) =>
              setSettings({ ...settings, patientNamePadding: value })
            }
            key={"patientNamePadding"}
            min={"0"}
            max={patientWidth}
            step={"1"}
          />
          <RangeOption
            label={"gradient width"}
            value={gradWidth}
            setValue={(value) => setSettings({ ...settings, gradWidth: value })}
            key={"gradWidth"}
            min={"0"}
            max={"200"}
            step={"1"}
          />
          <RangeOption
            label={"gradient height"}
            value={gradHeight}
            setValue={(value) =>
              setSettings({ ...settings, gradHeight: value })
            }
            key={"gradientHeight"}
            min={"10"}
            max={"150"}
            step={"1"}
          />
        </Grid>
        <Grid container item spacing={3}>
          <TextField
            required
            id="outlined-required"
            label="Pre"
            onChange={(event) =>
              setSettings({ ...settings, preLabel: event.target.value })
            }
            value={preLabel}
          />
          <TextField
            required
            id="outlined-required"
            label="Post"
            onChange={(event) =>
              setSettings({ ...settings, postLabel: event.target.value })
            }
            value={postLabel}
          />
          <TextField
            required
            id="outlined-required"
            label="P Value"
            value={pLabel}
            onChange={(event) =>
              setSettings({ ...settings, pLabel: event.target.value })
            }
          />
        </Grid>
        <Grid container item style={{ marginTop: 20 }}>
          <div style={{ padding: 10 }}>
            <p>Remove P</p>
            <Checkbox
              value={settings.removePval}
              onChange={() =>
                setSettings({
                  ...settings,
                  removePval: !settings.removePval,
                })
              }
              color="default"
            />
          </div>
          <div style={{ padding: 10 }}>
            <p>Switch All</p>
            <Checkbox
              value={settings.switchAll}
              onChange={() =>
                setSettings({
                  ...settings,
                  switchAll: !settings.switchAll,
                })
              }
              color="default"
            />
          </div>
          <div style={{ padding: 10 }}>
            <p>colour scheme</p>
            <Checkbox
              label="Change colour scheme"
              value={settings.colourScheme}
              onChange={() =>
                setSettings({
                  ...settings,
                  colourScheme: !settings.colourScheme,
                })
              }
              color="default"
            />
          </div>
          <div style={{ padding: 10 }}>
            <p>black outline</p>
            <Checkbox
              label="Black Outline"
              value={settings.blackOutline}
              onChange={() =>
                setSettings({
                  ...settings,
                  blackOutline: !settings.blackOutline,
                })
              }
              color="default"
            />
          </div>
          <button
            onClick={() => {
              var image = canvasRef.current
                .toDataURL("image/png")
                .replace("image/png", "image/octet-stream");

              window.location.href = image;
            }}
          >
            Save Image
          </button>
          <button
            onClick={() => {
              localStorage.setItem("heatmapSettings", JSON.stringify(settings));
            }}
          >
            Save Current Settings
          </button>
        </Grid>
      </Grid>
    </div>
  );
};

const drawLabels = (context, genes, geneYScale, settings) => {
  const { patientWidth, labelVertPadding, labelPadding, switchAll } = settings;
  context.fillStyle = "black";
  context.font = "bold 12px Helvetica";

  genes.map((gene) => {
    const x =
      parseFloat(labelPadding) + (switchAll ? parseFloat(patientWidth) : 0);
    const y = geneYScale(gene) + parseFloat(labelVertPadding);
    context.fillText(gene, x, y);
  });
};
const roundedSquares = (context, x, y, width, height, r) => {
  const radius = Math.min(Math.max(width - 1, 1), Math.max(height - 1, 1), r);

  context.save();
  context.lineJoin = "round";
  context.beginPath();
  context.moveTo(x + radius, y);
  context.lineTo(x + width - radius, y);
  context.quadraticCurveTo(x + width, y, x + width, y + radius);
  context.lineTo(x + width, y + height - radius);
  context.quadraticCurveTo(
    x + width,
    y + height,
    x + width - radius,
    y + height
  );
  context.lineTo(x + radius, y + height);
  context.quadraticCurveTo(x, y + height, x, y + height - radius);
  context.lineTo(x, y + radius);
  context.quadraticCurveTo(x, y, x + radius, y);
  context.stroke();
  context.restore();
  context.fill();
};
const addHeaders = (
  context,
  headersScale,
  startingX,
  topHeadingHeight,
  sqaureWidth,
  preLabel,
  postLabel,
  pLabel,
  pvalPadding,
  allXPadding
) => {
  headersScale.domain().map((d) => {
    const padding = d === "p" ? -pvalPadding : sqaureWidth / 2;
    const x = headersScale(d) + startingX + padding - allXPadding;
    //  const x = padding + startingX + headersScale(d) + sqaureWidth / 2;

    context.save();
    context.translate(x, topHeadingHeight - 10);
    context.rotate(-Math.PI / 4);
    context.textAlign = "left";
    context.font = "500 14px Helvetica";
    context.fillStyle = "black";
    const label = d === "Pre" ? preLabel : d === "Post" ? postLabel : pLabel;
    context.fillText(label, 0, 0);
    context.restore();
  });
};
const rangeFunc = (start, end, step) => {
  const len = Math.floor((end - start) / step) + 1;
  return Array(len)
    .fill()
    .map((_, idx) => start + idx * step);
};
const colorFlattenData = (keys, data) =>
  _.flatten(
    keys.reduce(
      (final, d) => [
        ...final,
        Object.keys(data[d]).reduce(
          (final, gene) => [
            ...final,
            data[d][gene]["Pre"],
            data[d][gene]["Post"],
          ],
          []
        ),
      ],
      []
    )
  ).filter((d) => d !== "nan");
const getHeatmapColorScaleByName = (data, colourScheme, name) => {
  const allStats = colorFlattenData(
    Object.keys(data).filter((d) => d === name),
    data
  );
  const allHeatmapExtent = d3.extent(allStats);
  return d3
    .scaleSequential(
      colourScheme ? d3.interpolatePlasma : d3.interpolateViridis
    )
    .domain([...allHeatmapExtent]);
};
const getPatientHeatmapColorScale = (data, colourScheme) => {
  const stats = colorFlattenData(
    Object.keys(data).filter((d) => d !== "All"),
    data
  );

  const heatmapExtent = d3.extent(stats);
  return d3
    .scaleSequential(
      colourScheme ? d3.interpolatePlasma : d3.interpolateViridis
    )
    .domain([...heatmapExtent]);
};
const setLightContextFont = (context) => {
  context.lineWidth = 1;
  context.globalAlpha = 1;
  context.font = "500 14px Helvetica";
  context.textAlign = "center";
  context.textBaseline = "middle";
};
const drawHeatmap = (
  context,
  patients,
  data,
  patientScale,
  geneYScale,
  settings
) => {
  const {
    patientWidth,
    patientHeight,
    radius,
    squareWidth,
    squareHeight,
    columnPadding,
    header,
    topHeadingHeight,
    preLabel,
    postLabel,
    pLabel,
    blackOutline,
    pValPadding,
    colourScheme,
    patientLegendPadding,
    allSectionPadding,
    allLegendPadding,
    switchAll,
    allLabelPadding,
    patientNamePadding,
    movePatientsPadding,
    gradWidth,
    gradHeight,
    removePval,
  } = settings;
  const allColorScale = getHeatmapColorScaleByName(data, colourScheme, "All");
  var colorScale = getPatientHeatmapColorScale(data, colourScheme);

  setLightContextFont(context);
  const columnHeadings = removePval ? ["Pre", "Post"] : ["Pre", "Post", "p"];

  const gradientHeight = parseFloat(gradHeight);
  const gradientWidth = parseFloat(gradWidth);

  patients.map((patient, i) => {
    const startingXGrad = patientScale(patient);
    const x =
      patient === "All"
        ? patientScale(patient) + patientWidth / 2 - 10 - allLabelPadding
        : patientScale(patient) +
          patientWidth / 2 -
          10 +
          parseFloat(patientNamePadding);
    const y = patientHeight - 50;
    setLightContextFont(context);
    //patient name
    context.fillText(patient, x, y);

    if (patient === "All") {
      const padding = 10 - parseFloat(allLegendPadding);

      const allX = padding + parseFloat(startingXGrad);

      var axis = context.createLinearGradient(
        allX,
        90,
        allX + gradientWidth,
        90
      );
      const step = Math.abs(
        (allColorScale.domain()[0] - allColorScale.domain()[1]) / 10
      );
      const range = rangeFunc(
        allColorScale.domain()[0],
        allColorScale.domain()[1],
        step
      );
      const colorRangeStep = 1 / 10;
      const colourRange = rangeFunc(0, 1, colorRangeStep);

      range.map((d, i) => {
        const stop = colourRange[i];
        axis.addColorStop(stop, allColorScale(d));
      });

      context.fillStyle = axis;

      const topCornerX = startingXGrad + padding;
      const topCornerY = y + 18;

      context.fillRect(topCornerX, topCornerY, gradientWidth, gradientHeight);
      context.fill();
      context.fillStyle = "black";
      context.strokeStyle = "black";
      context.font = "bold 10px Helvetica";
      //min
      context.beginPath();
      context.moveTo(topCornerX, topCornerY);
      context.lineTo(topCornerX, topCornerY + 19);
      context.fillText(
        d3.format(",d")(allColorScale.domain()[0]),
        topCornerX,
        topCornerY + 25
      );
      context.stroke();
      //max
      context.beginPath();
      context.moveTo(topCornerX + gradientWidth + 1, topCornerY);
      context.lineTo(topCornerX + gradientWidth + 1, topCornerY + 19);
      context.fillText(
        d3.format(",.1f")(allColorScale.domain()[1]),
        topCornerX + gradientWidth,
        topCornerY + 25
      );
      context.stroke();
    } else if (i === patients.length - 1) {
      const startingXGradPatient =
        parseFloat(startingXGrad) + parseFloat(movePatientsPadding);

      const padding = 10;
      var axis = context.createLinearGradient(
        startingXGradPatient + padding - patientLegendPadding,
        90,
        startingXGradPatient + gradientWidth + padding - patientLegendPadding,
        90
      );
      const step = Math.abs(
        (colorScale.domain()[0] - colorScale.domain()[1]) / 10
      );
      const range = rangeFunc(
        colorScale.domain()[0],
        colorScale.domain()[1],
        step
      );
      const colorRangeStep = 1 / 10;
      const colourRange = rangeFunc(0, 1, colorRangeStep);

      range.map((d, i) => {
        const stop = colourRange[i];
        axis.addColorStop(stop, colorScale(d));
      });

      context.fillStyle = axis;
      const topCornerX = startingXGradPatient + padding - patientLegendPadding;
      const topCornerY = y + 18;
      context.fillRect(topCornerX, topCornerY, gradientWidth, gradHeight);
      context.fill();
      context.fillStyle = "black";
      context.strokeStyle = "black";
      context.font = "bold 10px Helvetica";
      //min
      context.beginPath();
      context.moveTo(topCornerX, topCornerY);
      context.lineTo(topCornerX, topCornerY + 19);
      context.fillText(
        d3.format(",d")(colorScale.domain()[0]),
        topCornerX,
        topCornerY + 25
      );
      context.stroke();
      //max
      context.beginPath();
      context.moveTo(topCornerX + gradientWidth + 1, topCornerY);
      context.lineTo(topCornerX + gradientWidth + 1, topCornerY + 19);
      context.fillText(
        d3.format(",.1f")(colorScale.domain()[1]),
        topCornerX + gradientWidth,
        topCornerY + 25
      );
      context.stroke();
    }
  });
  context.beginPath();
  context.lineWidth = 1;
  context.globalAlpha = 1;
  context.font = "bold 14px Helvetica";
  context.textAlign = "center";
  context.textBaseline = "middle";
  context.fillStyle = "black";
  const perPatient = d3
    .scaleBand()
    .domain([...columnHeadings])
    .range([0, patientWidth])
    .paddingInner(columnPadding);

  const genes = Object.keys(data[patients[0]]);

  patients.map((patient) => {
    const startingX =
      patient === "All"
        ? patientScale(patient)
        : patientScale(patient) + parseFloat(movePatientsPadding);
    const currData = data[patient];
    context.font = "500 12px Helvetica";
    const allXPadding =
      patient === "All" && switchAll ? patientWidth - allSectionPadding : 0;

    addHeaders(
      context,
      perPatient,
      startingX,
      topHeadingHeight,
      squareWidth,
      preLabel,
      postLabel,
      pLabel,
      pValPadding,
      allXPadding,
      switchAll
    );

    genes.map((gene) => {
      const x = startingX;
      const y = geneYScale(gene);
      const currGene = currData[gene];

      perPatient.domain().map((term) => {
        const xBox = perPatient(term) + x - allXPadding;
        if (term !== "p") {
          context.fillStyle =
            patient === "All"
              ? allColorScale(currGene[term])
              : colorScale(currGene[term]);

          if (radius > 0) {
            context.strokeStyle = blackOutline
              ? "black"
              : colorScale(currGene[term]);
            roundedSquares(
              context,
              xBox,
              y,
              parseFloat(squareWidth),
              parseFloat(squareHeight),
              radius
            );
          } else {
            context.fillRect(xBox, y, squareWidth, squareHeight);
          }
        } else {
          context.fillStyle = "black";
          context.font = "bold 25px Helvetica";
          const pval = parseFloat(currGene[term]);
          const pX = xBox - pValPadding;
          const pvalLabel =
            pval <= 0.005
              ? "***"
              : pval <= 0.05
              ? "**"
              : pval <= 0.1
              ? "*"
              : "";
          const pvalY = y + squareHeight / 2 + 3;
          context.fillText(pvalLabel, pX, pvalY);
        }

        return;
      });
    });
  });
};
const RangeOption = ({ label, value, setValue, min, max, step, key }) => (
  <Grid item xs={2} sm={4} md={4} key={key}>
    <label>{label}</label>
    <input
      style={{ position: "relative" }}
      type="range"
      min={min}
      max={max}
      value={value}
      step={step}
      onChange={(event) => {
        setValue(event.target.value);
      }}
    />
  </Grid>
);
export default Heatmap;