import React, { useRef, useState, useMemo, Suspense } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import * as THREE from "three";
import "./App.css";
import {
  Environment,
  OrbitControls,
  PerspectiveCamera,
} from "@react-three/drei";
import Model from "./Model.js";

import Cellmine from "./Cellmine.js";

import five from "./assets/five.png";

const Box = (props) => {
  // This reference will give us direct access to the mesh
  const mesh = useRef();

  // Set up state for the hovered and active state
  const [active, setActive] = useState(false);

  // Rotate mesh every frame, this is outside of React without overhead
  /*  useFrame(() => {
    mesh.current.rotation.x = mesh.current.rotation.y += 0.01;
  });*/

  const texture = useMemo(() => new THREE.TextureLoader().load(five), []);
  /*    <mesh
      {...props}
      ref={mesh}
      scale={active ? [2, 2, 2] : [1.5, 1.5, 1.5]}
      onClick={(e) => setActive(!active)}
    >
      <boxBufferGeometry args={[1, 1, 1]} />
      <meshBasicMaterial attach="material" transparent side={THREE.DoubleSide}>
        <primitive attach="map" object={texture} />
      </meshBasicMaterial>    </mesh>*/
  return (
    <mesh
      {...props}
      ref={mesh}
      scale={active ? [2, 2, 2] : [1.5, 1.5, 1.5]}
      onClick={(e) => {
        console.log(e);
        setActive(!active);
      }}
    >
      <Cellmine />
    </mesh>
  );
};
/*    <hemisphereLight />
    <ambientLight intensity={0.5} />
    <spotLight position={[10, 10, 10]} angle={0.15} penumbra={1} />
    <pointLight position={[-10, -10, -10]} />
*/
const Blender = () => {
  return (
    <div style={{ width: "100vw", height: "100vh" }}>
      <Canvas flat linear>
        <ambientLight intensity={2} />
        <pointLight position={[40, 40, 40]} />
        <Box position={[2.5, 0, 0]} />
        <PerspectiveCamera makeDefault position={[10, 10, 5]} />
        <OrbitControls makeDefault />
      </Canvas>
    </div>
  );
};

export default Blender;
