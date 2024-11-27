[![Netlify Status](https://api.netlify.com/api/v1/badges/6fc6d371-87df-49b5-8e72-e1873fa5d54b/deploy-status)](https://app.netlify.com/sites/dune-copasi/deploys)

# DuneCopasi Website

This website is built using [Docusaurus 2](https://v2.docusaurus.io/), a modern static website generator.

### Installation

This command installs all the dependencies required to build the web page

```bash
yarn
```

### Local Development

This command starts a local development server and open up a browser window. Most changes are reflected live without having to restart the server.


```bash
yarn start
```

### Build

This command generates static content into the `build` directory and can be served using any static contents hosting service.

```bash
yarn build
```

### Deployment

This is automatically done by Netlify on every branch of the GitHub mirror. See
https://app.netlify.com/teams/soilros/builds/. In this case, we run a small script
[`build.sh`](build.sh) which checkout every git tag and generates different versions
of the documentation. Finally, it builds the documentation as a normal build.

```bash
./build.sh
```


### Versioning

#### Documentation
The usual command for versioning just works:

```bash
yarn docusaurus docs:version
```

This will copy the current docs into a folder and create a versioned sidebar.

#### Doxygen

Doxygen versions are automatically downloaded from the GitLab registry and simply put into the static assets.
If the doxygen documentation does not exist, links to that page (e.g. `DoxygenButton.js`) will fail to load.
This typically means that the versioned page must be done *after* the release is done.

#### Ini Files

Versioning does not include the asstes folder -_-, so we have to manage them manually.
Since we want the ini files for examples statically so that they can be served to web CLI, they are served
from static folder rather than from the docs folder. Links must also be updated manually.

#### Wasm Binary

This also needs to be done manually, just add the new version to the `package.json` and to the `WasmWorker.js`
