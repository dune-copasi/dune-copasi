const fs = require("fs");
const path = require("path");

async function setUpVersionRoutes(context, options) {
  return {
    name: "docusaurus-plugin-doxygen",
    async contentLoaded({ content, actions }) {
      const { createData, addRoute } = actions;

      // Path to the directory containing doxygen version folders
      const versionsDir = path.join(__dirname, "../../static/doxygen");

      // Function to list directories
      function listDirectories(dirPath) {
        return new Promise((resolve, reject) => {
          fs.readdir(dirPath, { withFileTypes: true }, (err, items) => {
            if (err) {
              return reject(err);
            }

            const directories = items
              .filter((item) => item.isDirectory())
              .map((folder) => folder.name);

            resolve(directories);
          });
        });
      }

      const versions = await listDirectories(versionsDir);

      for (const version of versions) {
        const versionJsonPath = await createData(
          `api/${version}.json`,
          JSON.stringify([`${version}`])
        );
        addRoute({
          path: `/api/${version}`,
          component: "@site/plugins/docusaurus-plugin-doxygen/components/APIVersionPage",
          exact: true,
          modules: {
            version: versionJsonPath,
          },
        });
      }
    },
  };
}

module.exports = setUpVersionRoutes;
