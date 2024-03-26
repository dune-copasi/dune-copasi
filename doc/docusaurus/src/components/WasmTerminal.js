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

  // instanciate a worker
  var wasmWorker = {}
  if (typeof Worker !== 'undefined') {
    wasmWorker = new Worker(new URL("@site/src/components/wasmworker.js", import.meta.url))
  }

  wasmWorker.onmessage = (e) => {
    console.log(e.data)
    if ("printText" in e.data) {
      let elem = document.getElementById(e.data.id)
      if (elem.innerHTML == "")
        elem.innerHTML = e.data.printText
      else
        elem.innerHTML = elem.innerHTML + "<br />" + e.data.printText
    }
    if ("error" in e.data) {
      alert("There was an unexpected error! See console log.")
      console.error(e.data.error)
    }
  }

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

      const id = `dune-copasi-${Date.now()}` 
      wasmWorker.postMessage({id, cmd: "cd", path})
      setCwd(newpath)
    },
    ls: (path) => {
      // guard against multiple paths
      if (path.split(" ").length > 1) 
        return "Error: usage: ls PATH"

      // no argument means cwd
      if (path.length === 0)
        path = cwd

      const id = `dune-copasi-${Date.now()}` 
      wasmWorker.postMessage({id, cmd: "ls", path})
      return <span id={id} />
    },
    mkdir: (path) => {
      // guard against multiple paths
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: mkdir DIR"

      // guard against invalid paths
      if (!instance.FS.analyzePath(path).parentExists)
        return `Error: parent of ${path} is not a valid path`

      // TODO: use wasmworker
      instance.FS.mkdir(path)
    },
    rmdir: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: rmdir DIR"

      // TODO: use wasmworker
      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`
      
      // TODO: use wasmworker
      try {
        instance.FS.rmdir(path)
      } catch (error) {
        return `Error: ${error}`
      }
    },
    rm: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: rm FILE"

      // TODO: use wasmworker
      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`
      
      // TODO: use wasmworker
      instance.FS.unlink(path)
    },
    edit: (path) => {
      if (path.split(" ").length > 1 || path.length === 0) 
        return "Error: usage: edit FILE"

      // TODO: use wasmworker
      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).exists)
        return `Error: ${path} is not a valid path`

      try {
      // TODO: use wasmworker
        setEditorText(instance.FS.readFile(path, {encoding: 'utf8'}))
      } catch (error) {
        return `Error: ${error}, ${error.msg}`
      }
    },
    save: (path) => {
      if (path.split(" ").length > 1 || path.length === 0)
        return "Error: usage: save FILE"

      // TODO: use wasmworker
      // guard against nonexistent paths
      if (!instance.FS.analyzePath(path).parentExists)
        return `Error: parent of ${path} is not a valid path`

      // TODO: use wasmworker
      instance.FS.writeFile(path, getEditorText())
    },
    "dune-copasi": (args) => {
      // generate a unique id for this dune-copasi call
      const id = `dune-copasi-${Date.now()}` 
      // wasmWorker.postMessage({id, args})
      return <span id={id}></span>
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