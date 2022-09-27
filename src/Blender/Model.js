import React, { useRef, useEffect } from "react";
import { useThree } from "@react-three/fiber";

import { useGLTF, useAnimations } from "@react-three/drei";

const deg2rad = (degrees) => degrees * (Math.PI / 180);
//var THREE = require("three");
//const raycaster = new THREE.Raycaster();
//const pointer = new THREE.Vector2();

//var THREE = require("three");
//var initializeDomEvents = require("threex.domevents");

export default function Model(props) {
  const group = useRef();
  const { nodes, materials, animations } = useGLTF("/untitled.glb");
  const { actions, mixer, names } = useAnimations(animations, group);

  /*  function onPointerMove(event) {
    console.log(event);
    // calculate pointer position in normalized device coordinates
    // (-1 to +1) for both components

    pointer.x = (event.clientX / window.innerWidth) * 2 - 1;
    pointer.y = -(event.clientY / window.innerHeight) * 2 + 1;
  }*/
  useThree(({ camera, raycaster, scene, mouse }) => {
    console.log(mouse);
    console.log(raycaster);
    //  camera.rotation.set(deg2rad(30), 0, 0);
  });
  useEffect(() => {
    /*  for (let i = 0; i < intersects.length; i++) {
      intersects[i].object.material.color.set(0xff0000);
    }

    renderer.render(scene, camera);*/

    actions[names[names.length - 1]].play();
  });
  //window.addEventListener("pointermove", onPointerMove);

  //  window.requestAnimationFrame(render);

  return (
    <>
      <group ref={group} {...props} dispose={null}>
        <group name="Scene">
          <group name="Charge" position={[-3.34, 0, 0]} />
          <mesh
            name="Suzanne"
            castShadow
            receiveShadow
            geometry={nodes.Suzanne.geometry}
            material={nodes.Suzanne.material}
          />
          <mesh
            name="Cube"
            castShadow
            receiveShadow
            geometry={nodes.Cube.geometry}
            material={nodes.Cube.material}
            position={[1.65, 1.39, -2.98]}
          />
        </group>
      </group>
    </>
  );
}

useGLTF.preload("/untitled.glb");
