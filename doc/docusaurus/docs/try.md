---
id: try
title: Try it online!
sidebar_label: Try it online!
---

import WasmInterface from '@site/src/components/wasminterface.js'

# WASM Interface

You can use DuneCopasi without installation directly here, in your browser.
This is still work in progress, not every feature might work in the browser and performance may be impacted.
Local installation is still recommended for intensive applications.

Upload at most one config file and optionally any number of data files and or type the configuration in the textarea.
Please be aware that all files will be placed in the same directory, so you might have to change the data file paths in the config file.
Then click on "compute" and as soon as the computation is done you can download the output.

<!--
    set of currently supported examples
    each key into `examples` identifies one example, 
    each corresponding object consists of paths to the config file and the data files
    Expand with more examples: 
      - add the files to the static/examples folder,
      - add a unique key for the example to the `examples` object,
      - for the new config, add paths relative to `static` folder as below
*/}
-->

<WasmInterface examples={{
    "gauss": {
        config: "/test/gauss.ini",
        dataFiles: [],
        dataFilesPath: ""
    },
    "mitchell_schaefer": {
        config: "/test/mitchell_schaefer.ini",
        dataFiles: [],
        dataFilesPath: ""
    },
    "two_disks": {
        config: "/test/two_disks.ini",
        dataFiles: ["/test/data/grids/two_disks.msh"],
        dataFilesPath: "data/grids"
    }
}}/>


You can visualize your results by opening them in [Glance](https://kitware.github.io/glance/app/) by Kitware and clicking the folder icon on the website.