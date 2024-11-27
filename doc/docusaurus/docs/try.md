---
id: try
title: Try it online!
sidebar_label: Try it online!
tags:
  - CLI
---

import WasmInterface from '@site/src/components/WasmInterface.js'

# Browser Interface

You can use DuneCopasi without installation directly here, in your browser, as you would use the executable on your machine.
Local installation is still recommended for intensive applications as performance might me impacted.

There are two "windows": the editor and the terminal.
The terminal is a limited UNIX-like environment, with a full virtual file system.
The editor is a VSCode like editor.
It is simply a buffer you can copy to and from (with ``edit`` and ``save``).

A typical workflow looks as follows:
  - ``upload`` a file from your machine to the virtual file system
  - ``edit FILE`` in the editor, this simply copies the file content into the editor, you have to manually save to a path after editing
  - ``save FILE`` to the virtual file system
  - ``dune-copasi --config=PATH --model.writer.vtk.path=OUTPUT`` run dune-copasi with a config file and overriding the output
  - ``download PATH`` from the virtual file system to your machine

You can visualize your results by opening them in [Glance](https://kitware.github.io/glance/app/) by Kitware and clicking the folder icon on the website.

<WasmInterface version="git"
  examples={{
    "/ini/next/heat.ini": "/dunecopasi/heat.ini",
    "/ini/next/grey_scott.ini": "/dunecopasi/grey_scott.ini",
    "/ini/next/volka_terra.ini": "/dunecopasi/volka_terra.ini"
  }}
/>
