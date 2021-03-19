import React from "react";
import * as d3 from "d3";
import _ from "lodash";
import { useDashboardState } from "../PlotState/dashboardState";

import Grid from "@material-ui/core/Grid";
import Paper from "@material-ui/core/Paper";

import { useCanvas } from "../components/utils/useCanvas";
import Info from "../Info/Info.js";
import infoText from "../Info/InfoText.js";
import { changeFontSize } from "../DrawingUtils/utils.js";

/*

Stacked bar chart of proportions (so they all add to 100%)

*/

const BAR_COLORS = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#ABDDA4",
  "#E6F598",
  "#FFFFBF",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#9E0142",
];

const LEGEND_WIDTH = 250;
const Y_AXIS_FONT = "normal 10px Helvetica";

const CAT_LABEL_SPACE = 100;
const CAT_LABEL_FONT = "bold 12px Helvetica";
const TOP_PADDING = 5;

const DataWrapper = ({ data, chartDim }) => {
  const [{ clonotypeParam, subtypeParam, fontSize }] = useDashboardState();

  const groupedSubtype = _.groupBy(data, subtypeParam);
  const subtypes = Object.keys(groupedSubtype).sort();

  const countedClonotypes = subtypes.reduce((countMap, subtype) => {
    const clonotypeCount = _.countBy(groupedSubtype[subtype], clonotypeParam);

    const countFreq = _.countBy(Object.values(clonotypeCount), (value) =>
      Math.min(value, 10)
    );

    return { ...countMap, [subtype]: countFreq };
  }, {});

  return (
    <StackedBarProportion
      data={countedClonotypes}
      chartDim={chartDim}
      barValues={Array.from(Array(10).keys()).map((value) => value + 1)}
    />
  );
};

const StackedBarProportion = ({ data, chartDim, barValues }) => {
  // make this flexible
  const categoryValues = Object.keys(data);

  const chartWidth = chartDim["width"] - LEGEND_WIDTH;
  const chartHeight = chartDim["height"] - CAT_LABEL_SPACE - TOP_PADDING;

  const catScale = d3
    .scaleBand()
    .domain(categoryValues)
    .range([0, chartWidth])
    .paddingInner(0.03)
    .paddingOuter(0.25);

  const barYScale = d3
    .scaleLinear()
    .domain([1, 0])
    .range([TOP_PADDING, chartHeight + TOP_PADDING]);
  const barSizeScale = d3
    .scaleLinear()
    .domain([0, 1])
    .range([0, chartHeight]);

  const colors = d3
    .scaleOrdinal()
    .domain(barValues)
    .range(BAR_COLORS.slice(0, barValues.length));

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");

      drawBars(
        context,
        data,
        categoryValues,
        barValues,
        catScale,
        barYScale,
        barSizeScale,
        colors
      );
      drawLabels(
        context,
        data,
        categoryValues,
        barValues,
        catScale,
        barYScale,
        colors
      );
    },

    chartDim["chart"]["x2"] - chartDim["chart"]["x1"],
    chartDim["chart"]["y2"] - chartDim["chart"]["y1"],
    []
  );

  return (
    <div
      style={{
        width: chartDim["width"],
        height: chartDim["height"],
        position: "relative",
      }}
    >
      <div
        id="barchart"
        style={{
          position: "absolute",
          pointerEvents: "all",
          display: "flex",
        }}
      >
        <canvas ref={ref} />
      </div>
    </div>
  );
};

const drawBars = (
  context,
  data,
  categoryValues,
  barValues,
  catScale,
  barYScale,
  barSizeScale,
  colors
) => {
  categoryValues.map((cValue) => {
    const categoryData = data[cValue];
    const total = Object.values(categoryData).reduce((sum, x) => sum + x, 0);

    var currHeight = barYScale(0);
    const xPos = catScale(cValue);

    barValues.map((bValue) => {
      if (categoryData.hasOwnProperty(bValue)) {
        context.fillStyle = colors(bValue);
        const barHeight = barSizeScale(categoryData[bValue] / total);

        context.fillRect(
          xPos,
          currHeight - barHeight,
          catScale.bandwidth(),
          barHeight
        );

        currHeight += barHeight * -1;
      }
    });
  });
};

const drawLabels = (
  context,
  data,
  categoryValues,
  barValues,
  catScale,
  barScale,
  colors
) => {
  const values = barScale.ticks(10);
  context.font = Y_AXIS_FONT;
  context.fillStyle = "black";
  context.textAlign = "right";
  context.lineWidth = 1;
  context.textBaseline = "middle";

  values.forEach((value) => {
    context.fillText(
      value * 100,
      catScale(categoryValues[0]) - 2,
      barScale(value)
    );
  });

  context.font = CAT_LABEL_FONT;

  categoryValues.forEach((cValue) => {
    context.save();
    context.translate(
      catScale(cValue) + catScale.bandwidth() / 2,
      barScale(0) + 7
    );
    context.rotate((322 * Math.PI) / 180);
    if (cValue.indexOf("/") !== -1) {
      context.fillText(cValue.split("/")[0] + "/", -5, 0);
      context.fillText(cValue.split("/")[1], -5, 10);
    } else {
      context.fillText(cValue, -5, 5);
    }
    context.stroke();
    context.restore();
  });
};

const StackedBar = ({ chartName, data, chartDim }) => {
  const [{ clonotypeParam, subtypeParam, fontSize }] = useDashboardState();

  const barWidth = 67;
  const groupedData = _.groupBy(data, subtypeParam);
  // console.log(groupedData);
  const subtypes = Object.keys(groupedData).sort();
  const stackedBarData = subtypes.reduce((final, subtype) => {
    const groupedClonotypes = _.groupBy(groupedData[subtype], clonotypeParam);
    // console.log(_.countBy(groupedData[subtype], clonotypeParam));
    var counter = {};
    const hitList = Object.entries(groupedClonotypes).map(
      (cellHits) => cellHits[1].length
    );
    // console.log(hitList);
    hitList.forEach((x) => (counter[x] = (counter[x] || 0) + 1));
    final[subtype] = {
      total: hitList.length,
      counts: counter,
    };
    return final;
  }, []);

  // X axis
  const x = d3
    .scaleBand()
    .domain(subtypes)
    .range([chartDim["chart"]["x1"], chartDim["chart"]["x2"] - 30]);
  // Y axis
  const y = d3
    .scaleLinear()
    .domain([100, 0])
    .range([chartDim["chart"]["y1"], chartDim["chart"]["y2"]]);

  var colors = d3
    .scaleOrdinal()
    .domain([...Array.from(Array(10).keys())])
    .range(BAR_COLORS);

  const ref = useCanvas(
    (canvas) => {
      const context = canvas.getContext("2d");
      drawLegend(context);
      drawBars(context);
      drawAxisLabels(context);
      drawYAxisLabels(context);
    },
    chartDim["width"],
    chartDim["height"],
    []
  );

  function drawBars(context) {
    context.beginPath();
    context.lineWidth = 1;
    context.strokeStyle = "black";
    subtypes.forEach((subtype, subIndex) => {
      var currentHeight = 0;
      [...Array.from(Array(10).keys())].map((key, index) => {
        const { counts, total } = stackedBarData[subtype];
        var height;
        if (index === 9) {
          const allOther = Object.entries(counts)
            .filter((row) => row[0] > 9)
            .map((row) => row[1]);
          height =
            allOther.length > 0
              ? (allOther.reduce((a, b) => a + b) / total) * 100
              : 0;
        } else {
          height = counts[key + 1] ? (counts[key + 1] / total) * 100 : 0;
        }
        context.fillStyle = colors(key);
        context.fillRect(
          x(subtype),
          y(height + currentHeight),
          barWidth,
          y(0) - y(height)
        );
        currentHeight += height;
        context.fill();
      });
    });
  }
  function drawLegend(context) {
    changeFontSize(context, fontSize.legendFontSize);
    [...Array.from(Array(10).keys())]
      .sort((a, b) => b - a)
      .map((key, index) => {
        context.fillStyle = colors(key);
        context.fillRect(
          chartDim["chart"]["x2"] - 20,
          chartDim["chart"]["y1"] + index * 14 + index * 2,
          9,
          9
        );
        context.fillStyle = "#000000";
        const legendText = key + 1 >= 10 ? "≥10" : key + 1;
        context.fillText(
          legendText,
          chartDim["chart"]["x2"] - 5,
          chartDim["chart"]["y1"] + index * 14 + index * 2 + 8
        );
        context.fill();
      });

    context.fillRect(
      chartDim["chart"]["x2"] + 20,
      chartDim["chart"]["y1"],
      5,
      5
    );
  }

  function drawYAxisLabels(context) {
    context.fillStyle = "black";
    context.textAlign = "right";
    context.lineWidth = 1;
    context.textBaseline = "middle";
    const ticks = y.ticks(10);

    context.beginPath();
    ticks.forEach(function(d) {
      changeFontSize(context, fontSize["tickLabelFontSize"]);
      //  context.moveTo(chartDim["margin"]["left"] + 17, y(d));
      //  context.lineTo(chartDim["margin"]["left"] + 27, y(d));
      context.fillText(d, chartDim["margin"]["left"] + 2, y(d));
      context.stroke();
    });
  }
  function drawAxisLabels(context) {
    context.beginPath();
    context.globalAlpha = 1;
    context.lineWidth = 1;
    context.fillStyle = "black";
    context.textAlign = "right";

    changeFontSize(context, fontSize["axisLabelFontSize"]);
    subtypes.map((subtype) => {
      context.save();
      context.translate(x(subtype) + barWidth / 2, y(0) + 7);
      context.rotate((322 * Math.PI) / 180);
      if (subtype.indexOf("/") !== -1) {
        context.fillText(subtype.split("/")[0] + "/", -5, 0);
        context.fillText(subtype.split("/")[1], -5, 10);
      } else {
        context.fillText(subtype, -5, 5);
      }
      context.stroke();
      context.restore();
    });
  }
  //style={{ margin: 10, width: chartDim["width"] }}
  return (
    <Paper style={{ margin: 10 }}>
      <Grid
        container
        direction="row"
        justify="flex-start"
        alignItems="flex-start"
        style={{
          width: chartDim["width"],
          height: chartDim["height"],
          position: "relative",
        }}
      >
        <Grid
          item
          xs={17}
          sm={8}
          id="barchart"
          style={{
            pointerEvents: "all",
            paddingRight: 0,
          }}
        >
          <canvas ref={ref} />
        </Grid>
        <Grid
          item
          xs={7}
          sm={4}
          style={{
            textAlign: "right",
            marginTop: 10,
            paddingRight: 15,
            pointerEvents: "all",
            curser: "pointer",
            zIndex: 100,
          }}
        >
          {infoText[chartName]["title"] + "    "}

          <Info name={chartName} direction="s" />
        </Grid>
      </Grid>
    </Paper>
  );
};
export default StackedBar;
