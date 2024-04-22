import BrowserOnly from '@docusaurus/BrowserOnly'
import React from "react"

function basename(path) {
  if (path === "/")
    return "/"
  else
    return path.split('/').reverse()[0];
}

function offerDownload(filename, content) {
    var blob = new Blob([content], { type: 'application/octet-stream' })

    // Create a temporary link element and trigger a download
    var a = document.createElement('a')
    a.href = URL.createObjectURL(blob)
    a.download = filename
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
}

const offerUpload = async () => {
  var input = document.createElement("input")
  input.type = "file"
  input.style.display = "none"
  document.body.appendChild(input)
  input.click()
  var resolve
  input.addEventListener('change', function() {
    resolve(input.files)
  })
  return new Promise((res, rej) => resolve = res)
}

// NOTE: need this extra class to call functions on mount
//       as we cannot call `useEffect` inside BrowserOnly
class WasmXtermImpl extends React.Component {
  constructor(props) {
    super(props)
    this.xtermRef = React.createRef()
    this.home = "/dunecopasi"
    this.state = {
      input: "",
      cwd: this.home,
      cursorPos: 0,
      running: false,
      history: [],
      historyCursor: 0,
    }
    this.wasmWorker = {}

    this.wasmCommands = {
      print: (data) => {
        if (!this.state.running)
          console.error("Not running but received a print!")
        this.writeln(data.msg.replace(/\n/g, "\r\n"))
      },
      exit: (data) => {
        this.setState({running: false}, () => this.prompt())
      },
      cwd: (data) => {
        this.setState({
          cwd: data.cwd
        })
      },
      edit: (data) => {
        this.props.setEditorText(data.text)
      },
      getEditorText: (data) => {
        this.wasmWorker.postMessage({cmd: "response", id: data.id, resolve: this.props.getEditorText()})
      },
      download: (data) => {
        offerDownload(data.filename, data.blob)
      },
      upload: async (data) => {
        const files = await offerUpload()
        this.wasmWorker.postMessage({cmd: "response", id: data.id, resolve: files})
      },
      open: async (data) => {
        const types = {
          "svg": "image/svg+xml",
          "png": "image/png",
          "jpeg": "image/jpeg",
          "jpg": "image/jpeg",
          "txt": "text/plain",
          "ini": "text/plain",
        }
        const extension = data.name.split(".").pop()
        if (!extension in types) {
          this.writeln(`Error: extension '${extension}' not found in ${types}`)
        } else {
          const file = new File([data.buf], data.name, {type: types[extension]})
          const url = window.URL.createObjectURL(file)
          window.open(url)
          window.URL.revokeObjectURL(url)
        }
      }
    }
  }

  write(t) {
    this.xtermRef.current.terminal.write(t)
  }

  writeln(t) {
    this.xtermRef.current.terminal.writeln(t)
  }

  newline() {
    this.write("\r\n")
  }

  prompt() {
    const dir = this.state.cwd === this.home ? "~" : basename(this.state.cwd)
    this.write(`${dir} $ `)
  }

  componentDidMount() {
    window.addEventListener("beforeunload", (event) => {
        event.preventDefault()
        event.returnValue = true
    })

    // instantiate a worker
    if (typeof Worker !== 'undefined') {
      this.wasmWorker = new Worker(new URL("@site/src/components/WasmWorker.js", import.meta.url))
    }

    this.write("Welcome to DuneCopasi shell v")
    this.setState({running: true}, () => {
      this.wasmWorker.postMessage({
        cmd: "instantiate", 
        args: {
          version: this.props.version, 
          files: this.props.examples,
        }
      })
    })

    this.wasmWorker.onmessage = (e) => {
      const data = e.data
      const cmd = data.cmd

      if (cmd in this.wasmCommands) {
        try {
          this.wasmCommands[cmd](data)
        } catch (error) {
          this.writeln(`Error in processing wasm command: ${error}`)
        }
      }
       else {
        this.writeln(`WasmXterm: Error: received unknown command ${cmd} from WasmWorker`)
      }

    }
  }

  render() {
    // NOTE: cannot import xterm globally as it accesses `document` which fails during SSG
    const XTerm = require("xterm-for-react").XTerm

    const onKey = (e) => {
      const printable = !e.altKey && !e.altGraphKey && !e.ctrlKey && !e.metaKey
      const key = e.domEvent.key

      if (key === "Enter") {
        this.newline()
        // NOTE: we do not handle space
        const inputSplit = this.state.input.split(" ").filter(arg => arg !== "")
        this.setState((state, props) => {return {
          input: "",
          cursorPos: 0,
          running: inputSplit.length > 0,
          history: [...state.history, state.input],
          historyCursor: state.history.length+1, 
        }}, () => {
          if (inputSplit.length > 0) {
            const [cmd, ...args] = inputSplit
            this.wasmWorker.postMessage({cmd, args})
          } else {
            this.prompt()
          }
        })
      }
      else if (key === "Backspace") {
        if (this.state.input.length > 0 && this.state.cursorPos > 0) {
          // move cursor back, override buffer
          // update state
          this.setState((state,props) => { 
            this.write("\b")
            this.write(state.input.slice(state.cursorPos))
            this.write(" ")
            this.write("\b".repeat(state.input.length - state.cursorPos + 1))
            return {
              input: state.input.slice(0,state.cursorPos-1) + state.input.slice(state.cursorPos),
              cursorPos: state.cursorPos - 1,
              historyCursor: state.history.length,
            }
          })
        }
      }
      else if (key === "ArrowLeft") {
        if (this.state.cursorPos > 0) {
          // move cursor back
          this.write("\b")
          this.setState((state, props) => {return {
            cursorPos: state.cursorPos-1,
          }})
        }
      }
      else if (key === "ArrowRight") {
        if (this.state.cursorPos < this.state.input.length) {
          // write what is currently under cursor (moves cursor forward)
          this.write(this.state.input[this.state.cursorPos])
          this.setState((state, props) => {return {
            cursorPos: state.cursorPos+1,
          }})
        }
      }
      else if (key === "ArrowUp") {
        this.setState((state,props) => {
          const historyCursor = state.historyCursor > 0 ? state.historyCursor-1 : 0
          const input = historyCursor < state.history.length ? state.history[historyCursor] : state.input
          // clear line
          this.write("\b \b".repeat(state.input.length))
          this.write(input)
          return {
            historyCursor,
            input,
            cursorPos: input.length,
          }
        })
      }
      else if (key === "ArrowDown") {
        this.setState((state,props) => {
          const historyCursor = state.historyCursor < state.history.length ? state.historyCursor+1 : state.history.length
          const input = historyCursor < state.history.length ? state.history[historyCursor] : state.input
          // clear line
          this.write("\b \b".repeat(state.input.length))
          this.write(input)
          return {
            historyCursor,
            input,
            cursorPos: input.length,
          }
        })
      }
      else if (e.domEvent.ctrlKey && key === "u") {
        this.write("\b \b".repeat(this.state.input.length))
        this.setState((state,props) => {return {
          input: "",
          cursorPos: 0,
          historyCursor: state.history.length,
        }})
      }
      else if (key === "Tab") {
        // unsupported currently
      }
      else if (printable) {
        this.setState((state,props) => {
          this.write(e.key + state.input.slice(state.cursorPos))
          this.write("\b".repeat(state.input.length-state.cursorPos))
          return {
            input: state.input.slice(0,state.cursorPos) + e.key + state.input.slice(state.cursorPos),
            cursorPos: state.cursorPos+1,
          }
        })
      }
    }

    return (
      <XTerm
        ref={this.xtermRef}
        onKey={onKey}
      />
    )
  }
}

export default function WasmXterm(props) {
  return (<BrowserOnly fallback={<div>Loading terminal...</div>}>
    {() => {
      return <WasmXtermImpl {...props} />
    }}
  </BrowserOnly>)
}