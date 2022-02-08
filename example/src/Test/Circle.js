import { Component } from "react-regl";

// ReglComponents components support a `drawCommand` Method which return a regl draw command

export default class Circle extends Component {
  drawCommand(regl) {
    return regl({
      frag: `
        precision mediump float;
        //attribute vec4 color;
        //uniform vec4 color;
        varying vec4 fcolor;
        void main () {
          gl_FragColor = fcolor;
        }`,

      vert: `
        precision mediump float;
        attribute vec2 position;
        attribute vec4 color;
        varying vec4 fcolor;

        // @change acquire the pointWidth uniform
        //  this is set by the uniforms section below
        uniform float pointWidth;
        uniform float stageWidth;
        uniform float stageHeight;

        vec2 normalizeCoords(vec2 position) {
        // read in the positions into x and y vars
        float x = position[0];
        float y = position[1];
        return vec2(
          2.0 * ((x / stageWidth) - 0.5),
          // invert y to treat [0,0] as bottom left in pixel space
          -(2.0 * ((y / stageHeight) - 0.5)));
        }

        void main () {
          // @change Set gl_PointSize global to
          //  configure point size
          fcolor = color;
          gl_PointSize = pointWidth;
          gl_Position = vec4(normalizeCoords(position), 0, 1);
        }`,

      attributes: {
        color: function (context, props) {
          return props.data.map((d) => {
            return JSON.parse(d["color"]);
          });
        },
        position: function (context, props, i) {
          return props.data.map((d) => {
            return [d.x, d.y];
          });
        },
      },

      uniforms: {
        stageWidth: 1000,

        stageHeight: 800,
        // @change: Add a pointWidth uniform -
        //  set by a prop
        pointWidth: (context, props) => 2,
        color: function (context, props) {
          return props.data.map((d) => {
            return JSON.parse(d["color"]);
          });
        },
      },

      count: function (context, props) {
        return props.data.length;
      },
      // @change: Set our primitive to points
      primitive: "points",
    });
  }
}
//export default Circle;
