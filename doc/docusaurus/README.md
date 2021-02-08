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
# WARNING: be aware that this modifies your git status and may screw your local changes!
./build.sh
```

