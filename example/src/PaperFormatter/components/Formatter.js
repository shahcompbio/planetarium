import React, { useState, useRef, useEffect } from "react";
import * as d3 from "d3";
import Grid from "@material-ui/core/Grid";

import { Stage, Layer, Image as ImageKonva, Text } from "react-konva";
import Konva from "konva";

import _ from "lodash";
import { useCanvas } from "@shahlab/planetarium";

const Formatter = ({ x, y, type }) => {
  const [context, saveContext] = useState(null);
  const [layer, saveLayer] = useState(null);
  const canvasWidth = x * 96;
  const canvasHeight = y * 96;
  const [settings, setSettings] = useState(
    localStorage.getItem("formatterSettings") === null
      ? {}
      : JSON.parse(localStorage.getItem("formatterSettings"))
  );

  const canvasRef = useRef();

  useEffect(() => {
    var stage = canvasRef.current.getStage();

    var l = new Konva.Layer();
    stage.add(l);
    stage.draw();
    saveLayer(l);
    // main API:
  }, []);

  /*  imageObj.onload = function () {
    var img = new Konva.Image({
      x: 50,
      y: 50,
      image: imageObj,
      width: 106,
      height: 118,
    });
    console.log(img);
    console.log("load");
    // add the shape to the layer
    layer.add(img);
  };*/
  //  imageObj.src = "/Users/vbojilova/Projects/planetarium/example/screenshot.png";

  const imageChange = (event) => {
    if (event.target.files && event.target.files[0]) {
      Array.from(event.target.files).map((file) => {
        var img = new Image();
        var URL = window.webkitURL || window.URL;
        var url = URL.createObjectURL(file);
        img.src = url;
        img.onload = function () {
          var img_width = img.width;
          var img_height = img.height;

          // calculate dimensions to get max 300px
          var max = 300;
          var ratio =
            img_width > img_height ? img_width / max : img_height / max;

          // now load the Konva image
          var theImg = new Konva.Image({
            image: img,
            x: 50,
            y: 30,
            width: img_width / ratio,
            height: img_height / ratio,
            draggable: true,
            dragBoundFunc: function (pos) {
              return {
                x: pos.x,
                y: pos.y,
              };
            },
          });
          var text = new Konva.Text({
            x: theImg.x() + 5,
            y: theImg.y() + 6,
            fontStyle: "bold",
            fontSize: 20,
          });
          var textBackground = new Konva.Rect({
            x: theImg.x(),
            y: theImg.y(),
            width: 120,
            height: 50,
            fill: "white",
            stroke: "black",
            strokeWidth: 1,
          });
          theImg.on("dragmove", function () {
            updateText();
          });
          theImg.on("transform", function () {
            updateText();
          });

          var tr1 = new Konva.Transformer({
            nodes: [theImg],
            keepRatio: false,
            rotateEnabled: false,
            enabledAnchors: [
              "top-left",
              "top-right",
              "bottom-left",
              "bottom-right",
            ],
          });

          layer.add(tr1);
          layer.add(theImg);
          layer.add(textBackground);
          layer.add(text);
          updateText();
          layer.draw();

          function updateText() {
            const width = Math.round(theImg.width() * theImg.scaleX());
            const height = Math.round(theImg.height() * theImg.scaleY());
            var lines = ["width: " + width, "height: " + height];

            text.text(lines.join("\n"));
            text.x(theImg.x() + 5);
            text.y(theImg.y() + 6);
            textBackground.x(theImg.x());
            textBackground.y(theImg.y());
          }
        };
      });
    }
  };

  return (
    <div>
      <Grid row spacing={1} style={{ width: 700 }}>
        <Grid>
          <Stage
            ref={canvasRef}
            width={canvasWidth}
            height={canvasHeight}
          ></Stage>
        </Grid>
        <Grid>
          <h1>Select Image</h1>
          <input
            type="file"
            name="myImage"
            onChange={imageChange}
            multiple="multiple"
          />
        </Grid>
      </Grid>
    </div>
  );
};
//  <canvas ref={canvasRef} id="formatter" />
export default Formatter;
