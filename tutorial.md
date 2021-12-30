# Tutorial Intro
This tutorial should show you how to integrate Planetarium into a new React project.

We will be using [yarn](https://yarnpkg.com/), but you can also use [npm](https://docs.npmjs.com/). We will also use [create-react-app](https://create-react-app.dev/), but you can use any React toolchain (including the [boilerplates](https://github.com/shahcompbio/viz-flask-boilerplate))


## Setup

Create a new create-react-app project

```
npx create-react-app my-app
cd my-app
```

Install planetarium

```
yarn add @shahlab/planetarium
```


## Add Planetarium component

In `App.js`
NOTE: may need to supply default data

```
import data from "./data.json";
import { Sankey } from "@shahlab/planetarium";

const App = () => {
  return (
    <Sankey
      data={data}
      width={800}
      height={700}
      subsetParam="cell_type"
      cloneParam="clone_id"
      timepointOrder={["Pre", "Post"]}
      timepointParam="treatment"
    />
  );
}

export default App;
```

Then you can run the app with `yarn start`

Celebrate!


## Interactive
The idea of this is to add state to keep track of what is interactivable. Here, we'll be adding that will highlight a UMAP plot when a Sankey node is clicked.

First, we will add the UMAP plot.

```
import React from 'react'
import data from "./data.json";
import { Sankey, UMAP } from "@shahlab/planetarium";

const App = () => {
    const [highlightNode, setHighlightNode] = useState(null)

  return (
    <div>
        <Sankey
            data={data}
            width={800}
            height={700}
            subsetParam="cell_type"
            cloneParam="clone_id"
            timepointOrder={["Pre", "Post"]}
            timepointParam="treatment"
        />
        <UMAP
            height={800}
            width={800}
            data={data}
            xParam={"UMAP_1"}
            yParam={"UMAP_2"}
            subsetParam={"cell_type"}
            idParam={"cell_id"}
            disable={true}
        />
    </div>
  );
}

export default App;
```

Second, we'll add state to track and filter the data, passing that filtered data over to the UMAP.

```
import React, { useState } from 'react'
import data from "./data.json";
import { Sankey } from "@shahlab/planetarium";

const App = () => {
    const [highlightNode, setHighlightNode] = useState(null)

    const highlightedIDs = highlightNode ? data.filter(datum => datum["treatment"] === highlightNode["treatment"] && datum["cell_type"] === highlightNode["cell_type"]).map(datum => datum["cell_id"]) : null


  return (
    <div>
        <Sankey
            data={data}
            width={800}
            height={700}
            subsetParam="cell_type"
            cloneParam="clone_id"
            timepointOrder={["Pre", "Post"]}
            timepointParam="treatment"
        />
        <UMAP
            height={800}
            width={800}
            data={data}
            highlightedIDs={highlightedIDs}
            xParam={"UMAP_1"}
            yParam={"UMAP_2"}
            subsetParam={"cell_type"}
            idParam={"cell_id"}
            disable={true}
        />
    </div>
  );
}

export default App;
```


Then add an event handler to set the state



```
import React, { useState } from 'react'
import data from "./data.json";
import { Sankey } from "@shahlab/planetarium";

const App = () => {
    const [highlightNode, setHighlightNode] = useState(null)

    const highlightedIDs = highlightNode ? data.filter(datum => datum["treatment"] === highlightNode["treatment"] && datum["cell_type"] === highlightNode["cell_type"]).map(datum => datum["cell_id"]) : null

  const onNodeClick = (node) => {
    if (!node) {
      setHighlightNode(null);
    } else {
      setHighlightNode(node);
    }
  };
  return (
    <div>
        <Sankey
            data={data}
            width={800}
            height={700}
            subsetParam="cell_type"
            cloneParam="clone_id"
            timepointOrder={["Pre", "Post"]}
            timepointParam="treatment"
            onNodeClick={onNodeClick}
        />
        <UMAP
            height={800}
            width={800}
            data={data}
            highlightedIDs={highlightedIDs}
            xParam={"UMAP_1"}
            yParam={"UMAP_2"}
            subsetParam={"cell_type"}
            idParam={"cell_id"}
            disable={true}
        />
    </div>
  );
}

export default App;
```

