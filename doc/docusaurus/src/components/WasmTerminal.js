import { useState, useEffect } from "react";
import { ReactTerminal, TerminalContextProvider } from "react-terminal";
import wasm from "dune-copasi-wasm-git"

// NOTE: update when command list changes
const helpMessage = `
help - print this message
whoami - print username
cd [PATH] - change current working directory to home or given path
ls [PATH] - list current working directory or given path
mkdir PATH - create directory
rmdir DIR - remove empty directory
rm FILE - remove file
edit FILE - load file into editor
save FILE - save editor contents to file
`

export default function WasmTerminal({setEditorText, getEditorText}) {
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
    help: () => helpMessage,
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
        return "Error: usage: ls PATH"

      // no argument means cwd
      if (path.length === 0)
        path = cwd

      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`

      return instance.FS.readdir(path).filter((dir) => dir !== "." && dir !== "..").join(" ")
    },
    mkdir: (path) => {
      // guard against multiple paths
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: mkdir DIR"

      // guard against invalid paths
      if (!instance.FS.analyzePath(path).parentExists)
        return `Error: parent of ${path} is not a valid path`

      instance.FS.mkdir(path)
    },
    rmdir: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: rmdir DIR"

      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`
      
      try {
        instance.FS.rmdir(path)
      } catch (error) {
        return `Error: ${error}`
      }
    },
    rm: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: rm FILE"

      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`
      
      instance.FS.unlink(path)
    },
    edit: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: edit FILE"

      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`

      try {
        setEditorText(instance.FS.readFile(path, {encoding: 'utf8'}))
      } catch (error) {
        return `Error: ${error}, ${error.msg}`
      }
    },
    save: (path) => {
      return getEditorText()
    }
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