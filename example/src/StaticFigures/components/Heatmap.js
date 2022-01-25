import React, { useState } from "react";
import * as d3 from "d3";
import Grid from "@material-ui/core/Grid";
import TextField from "@material-ui/core/TextField";
import Checkbox from "@material-ui/core/Checkbox";
import Collapse from "@material-ui/core/Collapse";

import List from "@material-ui/core/List";
import ListItem from "@material-ui/core/ListItem";
import ListSubheader from "@material-ui/core/ListSubheader";
import Typography from "@material-ui/core/Typography";

import _ from "lodash";
import { useCanvas } from "@shahlab/planetarium";
import CheckMarkSelect from "./CheckMarkSelect";
const sortingArr = [
  "CD74",
  "CIITA",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB5",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQA2",
  "HLA-DQB1",
  "HLA-DQB1-AS1",
  "HLA-DMA",
  "HLA-DMB",
  "HLA-DOA",
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "HLA-E",
  "HLA-G",
  "HLA-F",
  "B2M",
];
const Heatmap = ({ patients, data }) => {
  const [context, saveContext] = useState(null);
  const [inputSettings, setInputSettings] = useState({});
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
          pValVertPadding: 0,

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
          volume: false,
          volumeMaxHeight: 3,

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
    pValVertPadding,
  } = settings;

  const canvasHeight = patientHeight;
  const canvasWidth =
    (patients.length + 1) * patientWidth +
    (patients.length + 1) * patientSpacing;
  const [shownGenes, setShownGenes] = useState(
    Object.keys(data[patients[0]]).sort(
      (a, b) => sortingArr.indexOf(a) - sortingArr.indexOf(b)
    )
  );

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
      context.webkitImageSmoothingEnabled = false;
      context.mozImageSmoothingEnabled = false;
      context.imageSmoothingEnabled = false;
      saveContext(context);
      drawHeatmap(context, patients, data, patientScale, geneYScale, settings);

      drawLabels(context, shownGenes, geneYScale, settings);
    },
    canvasWidth,
    canvasHeight,
    [settings, shownGenes]
  );

  return (
    <div>
      <Grid container spacing={1} style={{ width: 700 }}>
        <FlexList title={"Omit Genes"}>
          {[
            <CheckMarkSelect
              allGenes={Object.keys(data[patients[0]])}
              selectedGenes={shownGenes}
              setSelected={setShownGenes}
            />,
          ]}
        </FlexList>

        <FlexList title={"Padding"}>
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
            label={"Second Column Padding"}
            value={columnPadding}
            setValue={(value) =>
              setSettings({ ...settings, columnPadding: value })
            }
            key={"columnPadding"}
            min={"0"}
            max={"1"}
            step={"0.1"}
          />
        </FlexList>
        <FlexList title={"Overall"}>
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
          <Typography>Width:{canvasWidth}</Typography>
          <Typography>Height:{canvasHeight}</Typography>
        </FlexList>
        <FlexList title={"Square"}>
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
        </FlexList>
        <FlexList title={"Patients"}>
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
        </FlexList>
        <FlexList title={"Header"}>
          {[
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
            />,
          ]}
        </FlexList>
      </Grid>

      <canvas
        ref={canvasRef}
        id="heatmapCanvas"
        width={canvasWidth * 300}
        height={canvasHeight * 300}
      />
      <Grid>
        <FlexList title={"All Section"}>
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
        </FlexList>
        <FlexList title={"Patient Section"}>
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
        </FlexList>
        <FlexList title={"Label Padding"}>
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
            max={"200"}
            step={"1"}
          />
        </FlexList>
        <FlexList title={"P Value"}>
          <RangeOption
            label={"Pval Padding"}
            value={pValPadding}
            setValue={(value) =>
              setSettings({ ...settings, pValPadding: value })
            }
            key={"pvalPadding"}
            min={"0"}
            max={"200"}
            step={"1"}
          />
          <RangeOption
            label={"Pval Vertical Padding"}
            value={pValVertPadding}
            setValue={(value) =>
              setSettings({ ...settings, pValVertPadding: value })
            }
            key={"pvalVertPadding"}
            min={"0"}
            max={"100"}
            step={"1"}
          />
        </FlexList>
        <FlexList title={"Gradient"}>
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
        </FlexList>
        <FlexList title={"Labels Text"}>
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
        </FlexList>
        <FlexList title={"Volume"}>
          <div style={{ padding: 10 }}>
            <p>Add Volume</p>
            <Checkbox
              value={settings.volume}
              onChange={() =>
                setSettings({
                  ...settings,
                  volume: !settings.volume,
                })
              }
              color="default"
            />
          </div>
          <RangeOption
            label={"Volume Max Height"}
            value={settings.volumeMaxHeight}
            setValue={(value) =>
              setSettings({ ...settings, volumeMaxHeight: value })
            }
            key={"volumeMaxHeight"}
            min={"0"}
            max={"100"}
            step={"1"}
          />
        </FlexList>
        <FlexList title={"Options"}>
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
        </FlexList>
        <FlexList title={"Settings"}>
          <button
            onClick={() => {
              localStorage.setItem("heatmapSettings", JSON.stringify(settings));
            }}
          >
            Save Current Settings
          </button>
          <button
            onClick={() => {
              window.prompt(
                "Copy to clipboard: Ctrl+C, Enter",
                JSON.stringify(settings)
              );
            }}
          >
            Export Current Settings
          </button>
          <textarea
            name="settings-input"
            cols="50"
            rows="10"
            id="settings-input"
            onChange={(event) => {
              setInputSettings({ ...JSON.parse(event.target.value) });
            }}
          >
            {}
          </textarea>
          <button
            onClick={() => {
              setSettings({ ...inputSettings });
            }}
          >
            Set Settings
          </button>
        </FlexList>
      </Grid>
    </div>
  );
};

const drawLabels = (context, genes, geneYScale, settings) => {
  const { patientWidth, labelVertPadding, labelPadding, switchAll } = settings;
  context.fillStyle = "black";
  context.textAlign = "right";
  context.font = "bold 16px Helvetica";

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
  )
    .filter((d) => d !== "nan")
    .map((d) => parseFloat(d));
const getGradientColorScaleByName = (data, name) => {
  const allStats = colorFlattenData(
    Object.keys(data).filter((d) => d === name),
    data
  );
  const allHeatmapExtent = d3.extent(allStats);

  //90%
  d3.piecewise(d3.interpolateHsl, [
    d3.rgb("#C58CE6"),
    d3.rgb("#A595E6"),
    d3.rgb("#937EE6"),
    d3.rgb("#E6878F"),
    d3.rgb("#E6494F"),
  ]);

  //darker
  d3.piecewise(d3.interpolateHsl, [
    d3.rgb("#996DB3"),
    d3.rgb("#8174B3"),
    d3.rgb("#7262B3"),
    d3.rgb("#B3696F"),
    d3.rgb("#B3393D"),
  ]);
  /*      d3.piecewise(d3.interpolateHsl, [
        d3.rgb("#A273BD"),
        d3.rgb("#B8A7FF"),
        d3.rgb("#EB97D7"),
        d3.rgb("#E6868F"),
        d3.rgb("#DD464C"),
      ])*/
  /*d3.interpolate("rgb(108,99,255)", "red")*/

  //d3.interpolate("rgb(108,99,255)", "#FFD653")

  return d3
    .scaleSequential(d3.interpolate("white", "#E03B2F"))
    .domain([...allHeatmapExtent]);
};
const getGradientByPatient = (data, name) => {
  const allStats = colorFlattenData(
    Object.keys(data).filter((d) => d !== "All"),
    data
  );
  const allHeatmapExtent = d3.extent(allStats);
  //  .scaleSequential(d3.piecewise(d3.interpolateHsl, ["#99165D", "#FFD653"]))
  //  .scaleSequential(d3.interpolate("rgb(108,99,255)", "#FFD653"))
  //"rgb(108,99,255)"
  return d3
    .scaleSequential(d3.interpolate("white", "#E03B2F"))
    .domain([...allHeatmapExtent]);
};
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
const getVolumeRange = (volumeMaxHeight, colorScale) => {
  const minMaxDiff = Math.abs(colorScale.domain()[0] - colorScale.domain()[1]);
  return d3
    .scaleLinear()
    .domain([0, minMaxDiff])
    .range([0, parseFloat(volumeMaxHeight)]);
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
    volume,
    volumeMaxHeight,
  } = settings;
  //  const allColorScale = getHeatmapColorScaleByName(data, colourScheme, "All");
  //var colorScale = getPatientHeatmapColorScale(data, colourScheme);
  const allColorScale = getGradientColorScaleByName(data, "All");
  const colorScale = getGradientByPatient(data);

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
    context.font = "normal 18px Helvetica";
    const patientLabel = patient === "All" ? "Combined (n = 4)" : "UPN" + i;
    context.fillText(patientLabel, x, y);
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
      context.font = "normal 10px Helvetica";
      //min
      context.lineWidth = 1;
      context.beginPath();
      context.moveTo(topCornerX, topCornerY);
      context.lineTo(topCornerX, topCornerY + 19);
      context.fillText(
        d3.format(",.3f")(allColorScale.domain()[0]),
        topCornerX,
        topCornerY + 25
      );
      context.stroke();
      //max
      context.beginPath();
      context.moveTo(topCornerX + gradientWidth + 1, topCornerY);
      context.lineTo(topCornerX + gradientWidth + 1, topCornerY + 19);
      context.fillText(
        d3.format(",.3f")(allColorScale.domain()[1]),
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
      context.font = "normal 10px Helvetica";
      //min
      context.lineWidth = 1;
      context.beginPath();
      context.moveTo(topCornerX, topCornerY);
      context.lineTo(topCornerX, topCornerY + 19);
      //",.3f"
      //",d"
      context.fillText(
        d3.format(",d")(colorScale.domain()[0]),
        topCornerX,
        topCornerY + 25
      );
      context.stroke();
      //max
      //",.3f"
      //",.1f"
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
  context.font = "bold 18px Helvetica";
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
    context.font = "500 18px Helvetica";
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

    const volumeRange = getVolumeRange(volumeMaxHeight, allColorScale);

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
            if (volume) {
              const diff = Math.abs(currGene["Pre"] - currGene["Post"]);

              const vol = patient === "All" ? volumeRange(diff) : 3;
              //  console.log(vol);
              if (currGene["Pre"] >= currGene["Post"] && term === "Pre") {
                roundedSquares(
                  context,
                  xBox + vol,
                  y + vol,
                  parseFloat(squareWidth),
                  parseFloat(squareHeight),
                  radius
                );
              } else if (
                currGene["Pre"] <= currGene["Post"] &&
                term === "Post"
              ) {
                roundedSquares(
                  context,
                  xBox + vol,
                  y + vol,
                  parseFloat(squareWidth),
                  parseFloat(squareHeight),
                  radius
                );
              }
            }
          } else {
            context.fillRect(xBox, y, squareWidth, squareHeight);
            if (volume) {
              context.fillRect(
                xBox + parseFloat(volumeMaxHeight),
                y + parseFloat(volumeMaxHeight),
                squareWidth,
                squareHeight
              );
            }
          }
        } else {
          context.fillStyle = "black";
          context.font = "500 20px Helvetica";
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
          /*  console.log(pval);
          console.log(currGene);
          if (pvalLabel !== "") {
            console.log(pvalLabel);
            console.log(currGene);
          }*/

          const pvalY = y + squareHeight / 2 + 3;
          context.fillText(pvalLabel, pX, pvalY);
        }

        return;
      });
    });
  });
};
const FlexList = (props) => (
  <Grid container item spacing={3} key={props.title + "-grid"}>
    <List
      key={props.title + "-list"}
      component="div"
      disablePadding
      style={{ display: "flex", flexDirection: "row", padding: 0 }}
      subheader={
        <ListSubheader
          style={{ lineHeight: "20px" }}
          key={props.title + "-label"}
        >
          {props.title}
        </ListSubheader>
      }
    >
      {props.children.map((child) => (
        <ListItem> {child} </ListItem>
      ))}
    </List>
  </Grid>
);
const RangeOption = ({ label, value, setValue, min, max, step, key }) => (
  <Grid item xs={2} sm={4} md={4} key={key}>
    <label style={{ lineHeight: "20px" }}>{label}</label>
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
