import BrowserOnly from '@docusaurus/BrowserOnly'
import React from "react"

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
    }
    this.wasmWorker = {}
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
    this.write(this.state.cwd + " $ ")
  }

  componentDidMount() {
    // Add the starting text to the terminal
    this.writeln(`Welcome to DuneCopasi shell ${this.props.version}`);
    this.prompt()

    // instanciate a worker
    if (typeof Worker !== 'undefined') {
      this.wasmWorker = new Worker(new URL("@site/src/components/wasmworker.js", import.meta.url))
    }

    this.wasmWorker.onmessage = (e) => {
      const data = e.data
      const cmd = data.cmd

      console.log("WasmXterm received:", data)

      if (cmd === "print") {
        if (!this.state.running)
          console.error("Not running but received a print!")
        this.writeln(data.msg)
      } else if (cmd === "exit") {
        this.setState({
          running: false,
        }, () => this.prompt())
      } else if (cmd === "cwd") {
        this.setState({
          cwd: data.cwd
        })
      } else if (cmd === "getEditorText") {
          this.wasmWorker.postMessage({cmd: "response", id: data.id, resolve: this.props.getEditorText()})
      } else {
        console.log(`WasmXterm: ${cmd} not implemented`)
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
        const [cmd, ...args] = this.state.input.split(" ").filter(arg => arg !== "")
        this.setState({
          input: "",
          cursorPos: 0,
          running: true,
        }, () => {
          console.log("WasmXterm sending:", {cmd, args})
          this.wasmWorker.postMessage({cmd, args})
        })
      }
      else if (key === "Backspace") {
        if (this.state.input.length > 0 && this.state.cursorPos > 0) {
          // move cursor back, override buffer
          // update state
          this.setState((state,props) => { 
            this.write("\b \b" + state.input.slice(state.cursorPos) + "\b".repeat(state.input.length - state.cursorPos))
            return {
              input: state.input.slice(0,state.cursorPos-1) + state.input.slice(state.cursorPos),
              cursorPos: state.cursorPos - 1,
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
      else if (e.domEvent.ctrlKey && key === "u") {
        this.write("\b \b".repeat(this.state.input.length))
        this.setState({
          input: "",
          cursorPos: 0,
        })
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