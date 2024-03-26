import { useState, useEffect } from "react";
import { ReactTerminal, TerminalContextProvider } from "react-terminal";
import wasm from "dune-copasi-wasm-git"

export default function WasmTerminal() {
  const theme = "custom";
  const themes = {
    custom: {
      themeBGColor: "rgb(41,45,62)",
      themeToolbarColor: "#ffffff1a",
      themeColor: "rgb(191, 199, 213)",
      themePromptColor: "rgb(130, 170, 255)",
    },
  }

  const home = "/dunecopasi"
  const [cwd, setCwd] = useState(home)

  // setup file system
  const preRun = (instance) => {
    instance.FS.mkdir(home)
    instance.FS.chdir(home)
  }

  // create wasm instance
  const [instance, setInstance] = useState(null)
  useEffect(() => {
    wasm({preRun}).then(inst => setInstance(inst))
  }, [])

  // Define commands here
  const commands = {
    whoami: "dunecopasi",
    where: () => {
      return instance.FS.cwd()
    },
    cd: (path) => {
      // guard against multiple paths
      if (path.split(" ").length > 1) 
        return "Error: usage: cd <path>"

      // no argument means home
      if (path.length === 0)
        path = startdir

      // switch on absolute vs relative path
      var newpath = path.startsWith("/") ? path : cwd + "/" + path

      // normalize path to not end with slash
      if (newpath.endsWith("/") && newpath !== "/")
        newpath = newpath.slice(0, -1)

      // guard against nonexistent paths
      if (!instance.FS.analyzePath(newpath).exists)
        return `Error: ${newpath} is not a valid path`

      setCwd(newpath)
      instance.FS.chdir(newpath)
    },
    ls: (path) => {
      // guard against multiple paths
      if (path.split(" ").length > 1) 
        return "Error: usage: ls <path>"

      // no argument means cwd
      if (path.length === 0)
        path = cwd

      return instance.FS.readdir(path).filter((dir) => dir !== "." && dir !== "..").join(" ")
    },
  }

  return (
    <TerminalContextProvider>
      <div style={{ height: "30vh" }}>
        <ReactTerminal
          prompt={`${cwd} $`}
          commands={commands}
          themes={themes}
          theme={theme}
        />
      </div>
    </TerminalContextProvider>
  )
}