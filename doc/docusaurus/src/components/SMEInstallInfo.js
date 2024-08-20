import React from "react";
import Admonition from "@theme/Admonition";
import CodeBlock from "@theme/CodeBlock";
import Tabs from "@theme/Tabs";
import TabItem from "@theme/TabItem";
import clsx from "clsx";

import styles from "./styles.module.css";
import Link from "@docusaurus/Link";

export default function SMEInstallInfo() {
  const button_cls = "button button--outline button--secondary button--lg";
  return (
    <div>
      <div>
        The{" "}
        <a href="https://spatial-model-editor.github.io/">
          {" "}
          Spatial Model Editor (SME){" "}
        </a>{" "}
        is a <b>user friendly</b> GUI editor to create and edit 2D and 3D
        spatial <i>Systems Biology Markup Language</i>{" "}
        <a href="https://en.wikipedia.org/wiki/SBML">(SBML)</a> models of
        bio-chemical reactions. Additionally, it can run simulations with{" "}
        DuneCopasi. A big adventage of this package is that is
        tailored for biologists and is available with just a click on the major
        plataforms. You can download it here:
      </div>
      <Tabs
        defaultValue="windows"
        values={[
          { label: "Windows", value: "windows" },
          { label: "MacOS Arm64", value: "macos-arm" },
          { label: "MacOS Intel", value: "macos-intel" },
          { label: "Linux", value: "linux" },
        ]}
      >
        <TabItem value="linux">
          <div className={styles.buttons}>
            <Link
              className={clsx(button_cls)}
              to="https://github.com/spatial-model-editor/spatial-model-editor/releases/latest/download/spatial-model-editor"
            >
              <svg
                xmlns="http://www.w3.org/2000/svg"
                version="1.1"
                x="0px"
                y="0px"
                width="100"
                height="100"
                viewBox="0 0 1000 1000"
              >
                <metadata>
                  {" "}
                  Svg Vector Icons : http://www.onlinewebfonts.com/icon{" "}
                </metadata>
                <g>
                  <g>
                    <g>
                      <path d="M896,809.4c-25.9-10.6-47.2-27.4-45.7-59.4c1.5-32-22.9-53.2-22.9-53.2s21.3-70.1,1.5-128c-19.8-58-85.3-150.8-135.5-220.9c-50.2-70.1-7.6-150.9-53.2-254.4C594.4-10,475.7-4,411.7,40.2c-64,44.1-44.2,153.8-41.2,205.6c3.1,51.7,1.4,88.6-4.5,102c-6,13.4-47.2,62.5-74.7,103.6c-27.4,41.1-47.2,126.4-67.1,161.5c-19.8,35-6.1,67-6.1,67s-13.7,4.5-24.4,27.5c-10.7,22.8-32,33.4-70.1,41c-38.1,7.6-38.1,32.1-29,59.5c9.2,27.4,0,42.7-10.6,77.6c-10.7,34.9,42.6,45.7,94.4,51.7c51.8,6.2,109.7,39.7,158.5,45.7c48.7,6.1,63.9-33.5,63.9-33.5s54.8-12.3,112.7-13.7c57.9-1.5,112.7,12.2,112.7,12.2s10.7,24.4,30.5,35c19.8,10.7,62.5,12.2,89.9-16.7c27.5-29,100.6-65.5,141.7-88.4C929.5,855,921.9,820,896,809.4z M539.6,137.7c26.1,0,47.2,25.9,47.2,57.8c0,22.7-10.6,42.2-26,51.7c-3.9-1.7-8.1-3.5-12.5-5.4c9.3-4.6,15.9-16.4,15.9-30.3c0-18-11.1-32.6-24.9-32.6c-13.6,0-24.8,14.6-24.8,32.6c0,6.7,1.6,13.1,4.3,18.3c-8.1-3.2-15.6-6.2-21.5-8.5c-3.2-7.8-5-16.6-5-25.9C492.3,163.6,513.4,137.7,539.6,137.7z M536.2,259.6c13.1,4.5,27.6,13,26.1,21.4c-1.5,8.5-8.4,8.5-26.1,19.3c-17.7,10.7-56,34.5-68.3,36c-12.3,1.5-19.2-5.4-32.3-13.8c-13.1-8.5-37.6-28.5-31.4-39.2c0,0,19.1-14.6,27.5-22.3c8.5-7.7,30-26.1,43-23.7C487.8,239.6,523.1,255,536.2,259.6z M418.4,146.8c20.6,0,37.3,24.5,37.3,54.8c0,5.6-0.5,10.7-1.6,15.7c-5,1.7-10.1,4.5-15.1,8.7c-2.5,2.1-4.7,4-6.9,5.9c3.3-6.1,4.6-14.8,3.1-24c-2.8-16.5-13.8-28.6-24.7-26.9c-10.9,1.9-17.5,16.7-14.7,33.4c2.8,16.6,13.8,28.7,24.6,26.9c0.6-0.1,1.2-0.3,1.8-0.5c-5.3,5.1-10.2,9.5-15.2,13.2c-15.1-7-26.1-27.8-26.1-52.4C381.2,171.3,397.9,146.8,418.4,146.8z M378.2,916.5c-4.9,21.8-30.4,37.7-30.4,37.7c-23.2,7.3-87.6-20.7-116.8-32.9c-29.2-12.1-103.4-15.9-113.2-26.7c-9.7-11,4.9-35.4,8.6-58.4c3.6-23.2-7.3-37.7-3.7-53.6c3.7-15.8,51.1-15.8,69.3-26.7c18.3-11,21.9-42.6,36.5-51.1c14.6-8.6,41.3,21.8,52.3,39c10.9,16.9,52.3,90,69.3,108.3C367.3,870.2,383.1,894.6,378.2,916.5z M647.6,704.2c-4.4,21.5-4.4,99.1-4.4,99.1s-47.2,65.4-120.4,76.1c-73.1,10.7-109.7,3-109.7,3L372,835.3c0,0,31.9-4.6,27.4-36.6c-4.6-32-97.5-76.2-114.2-115.8c-16.7-39.5-3-106.6,18.3-140.2c21.3-33.5,34.9-106.5,56.3-131c21.3-24.3,38-76.1,30.4-99c0,0,45.7,54.9,77.6,45.8c32-9.2,103.7-62.5,114.2-53.3c10.6,9.2,102,210.2,111.1,274.2c9.2,63.9-6.1,112.7-6.1,112.7S652.1,682.9,647.6,704.2z M881.4,847.7c-14.2,13.1-93.4,45-117.1,70c-23.6,24.7-54.4,44.9-73.3,39c-19-6-35.5-32-27.2-69.9c8.2-37.8,15.4-79.2,14.2-102.9c-1.2-23.7-6-55.7,0-60.4c5.9-4.6,15.3-2.3,15.3-2.3s-4.6,44.9,22.5,56.8c27.2,11.7,66.2-4.7,78.1-16.7c11.9-11.8,20.2-29.5,20.2-29.5s11.8,6,10.6,24.9c-1.2,18.9,8.2,46.2,26.1,55.6C868.4,821.7,895.6,834.8,881.4,847.7z" />
                    </g>
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                    <g />
                  </g>
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                  <g />
                </g>
              </svg>
              <br />
              Download
            </Link>
          </div>
          <br />
          <Admonition type="info" title="permissions">
            <CodeBlock>
              chmod +x spatial-model-editor <br />
            </CodeBlock>
            <Admonition type="caution" title="Fedora/RHEL/CentOS">
              On linux some additional system libraries are required that may
              not be installed by default. To install them:
              <CodeBlock>
                sudo yum install xcb-util-image xcb-util-keysyms
                xcb-util-renderutil xcb-util-wm
              </CodeBlock>
            </Admonition>
          </Admonition>
        </TabItem>
        <TabItem value="macos-arm">
          <div className={styles.buttons}>
            <Link
              className={clsx(button_cls)}
              to="https://github.com/ssciwr/sme-osx-arm64/releases/latest/download/spatial-model-editor.dmg"
            >
              <svg
                xmlns="http://www.w3.org/2000/svg"
                x="0px"
                y="0px"
                width="100"
                height="100"
                viewBox="0 0 50 50"
              >
                <path d="M 44.527344 34.75 C 43.449219 37.144531 42.929688 38.214844 41.542969 40.328125 C 39.601563 43.28125 36.863281 46.96875 33.480469 46.992188 C 30.46875 47.019531 29.691406 45.027344 25.601563 45.0625 C 21.515625 45.082031 20.664063 47.03125 17.648438 47 C 14.261719 46.96875 11.671875 43.648438 9.730469 40.699219 C 4.300781 32.429688 3.726563 22.734375 7.082031 17.578125 C 9.457031 13.921875 13.210938 11.773438 16.738281 11.773438 C 20.332031 11.773438 22.589844 13.746094 25.558594 13.746094 C 28.441406 13.746094 30.195313 11.769531 34.351563 11.769531 C 37.492188 11.769531 40.8125 13.480469 43.1875 16.433594 C 35.421875 20.691406 36.683594 31.78125 44.527344 34.75 Z M 31.195313 8.46875 C 32.707031 6.527344 33.855469 3.789063 33.4375 1 C 30.972656 1.167969 28.089844 2.742188 26.40625 4.78125 C 24.878906 6.640625 23.613281 9.398438 24.105469 12.066406 C 26.796875 12.152344 29.582031 10.546875 31.195313 8.46875 Z"></path>
              </svg>
              <br />
              Download
            </Link>
          </div>
          <br />
          <Admonition type="info" title="permissions">
            <b>Right-Click </b> &gt; <b>Open</b> &gt; <b>Ok</b>
          </Admonition>
        </TabItem>
        <TabItem value="macos-intel">
          <div className={styles.buttons}>
            <Link
              className={clsx(button_cls)}
              to="https://github.com/spatial-model-editor/spatial-model-editor/releases/latest/download/spatial-model-editor.dmg"
            >
              <svg
                xmlns="http://www.w3.org/2000/svg"
                x="0px"
                y="0px"
                width="100"
                height="100"
                viewBox="0 0 50 50"
              >
                <path d="M 44.527344 34.75 C 43.449219 37.144531 42.929688 38.214844 41.542969 40.328125 C 39.601563 43.28125 36.863281 46.96875 33.480469 46.992188 C 30.46875 47.019531 29.691406 45.027344 25.601563 45.0625 C 21.515625 45.082031 20.664063 47.03125 17.648438 47 C 14.261719 46.96875 11.671875 43.648438 9.730469 40.699219 C 4.300781 32.429688 3.726563 22.734375 7.082031 17.578125 C 9.457031 13.921875 13.210938 11.773438 16.738281 11.773438 C 20.332031 11.773438 22.589844 13.746094 25.558594 13.746094 C 28.441406 13.746094 30.195313 11.769531 34.351563 11.769531 C 37.492188 11.769531 40.8125 13.480469 43.1875 16.433594 C 35.421875 20.691406 36.683594 31.78125 44.527344 34.75 Z M 31.195313 8.46875 C 32.707031 6.527344 33.855469 3.789063 33.4375 1 C 30.972656 1.167969 28.089844 2.742188 26.40625 4.78125 C 24.878906 6.640625 23.613281 9.398438 24.105469 12.066406 C 26.796875 12.152344 29.582031 10.546875 31.195313 8.46875 Z"></path>
              </svg>
              <br />
              Download
            </Link>
          </div>
          <br />
          <Admonition type="info" title="permissions">
            <b>Right-Click </b> &gt; <b>Open</b> &gt; <b>Ok</b>
          </Admonition>
        </TabItem>
        <TabItem value="windows">
          <div className={styles.buttons}>
            <Link
              className={clsx(
                "button button--outline button--secondary button--lg",
              )}
              to="https://github.com/spatial-model-editor/spatial-model-editor/releases/latest/download/spatial-model-editor.exe"
            >
              <svg
                xmlns="http://www.w3.org/2000/svg"
                x="0px"
                y="0px"
                width="100"
                height="100"
                viewBox="0 0 50 50"
              >
                <path d="M19.852 7.761l-15 2.25C4.362 10.085 4 10.505 4 11v12c0 .553.448 1 1 1h15c.552 0 1-.447 1-1V8.75c0-.291-.127-.567-.348-.758C20.432 7.803 20.139 7.721 19.852 7.761zM45.652 4.242c-.22-.189-.512-.271-.801-.231l-21 3.15C23.362 7.235 23 7.655 23 8.15V23c0 .553.448 1 1 1h21c.552 0 1-.447 1-1V5C46 4.709 45.873 4.433 45.652 4.242zM20 26H5c-.552 0-1 .447-1 1v12c0 .495.362.915.852.989l15 2.25c.05.007.099.011.148.011.238 0 .47-.085.652-.242C20.873 41.817 21 41.541 21 41.25V27C21 26.447 20.552 26 20 26zM45 26H24c-.552 0-1 .447-1 1v14.85c0 .495.362.915.852.989l21 3.15C44.901 45.996 44.951 46 45 46c.238 0 .47-.085.652-.242C45.873 45.567 46 45.291 46 45V27C46 26.447 45.552 26 45 26z"></path>
              </svg>
              <br />
              Download
            </Link>
          </div>
          <br />
          <Admonition type="info" title="permissions">
            <b>More info</b> &gt; <b>Run anyway</b>
          </Admonition>
        </TabItem>
      </Tabs>
      <div>
        Or take a look at the Spatial Model Editor (SME){" "}
        <a href="https://spatial-model-editor.github.io/">
          website
        </a>{" "}
         or its {" "}
        <a href="https://spatial-model-editor.readthedocs.io/en/stable/quickstart/get-started.html">
          documentation
        </a>
        . Note that the version of DuneCopasi and the Spatial
        Model Editor are not synchronized!
      </div>
    </div>
  );
}
