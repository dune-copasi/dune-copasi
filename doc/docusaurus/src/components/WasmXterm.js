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
    }
  }

  componentDidMount() {
    // Add the starting text to the terminal
    this.xtermRef.current.terminal.writeln(
      `Welcome to DuneCopasi shell ${this.props.version}`
    );
    this.xtermRef.current.terminal.write("echo> ");
  }

  render() {
    // NOTE: cannot import xterm globally as it accesses `document` which fails during SSG
    const XTerm = require("xterm-for-react").XTerm
    const write = (t) => this.xtermRef.current.terminal.write(t)
    const writeln = (t) => this.xtermRef.current.terminal.writeln(t)
    const newline = () => write("\r\n")

    const prompt = () => write(this.state.cwd + " $ ")

    const onKey = (e) => {
      const printable = !e.altKey && !e.altGraphKey && !e.ctrlKey && !e.metaKey
      const key = e.domEvent.key
      console.log(e)

      if (key === "Enter") {
        newline()
        writeln(this.state.input)
        prompt()
        this.setState({ input: "", cursorPos: 0 })
      }
      else if (key === "Backspace") {
        if (this.state.input.length > 0 && this.state.cursorPos > 0) {
          // move cursor back, override buffer
          write("\b \b" + this.state.input.slice(this.state.cursorPos-1))
          // update state
          this.state((state,props) => { return {
            input: state.input.slice(0,state.cursorPos-1) + input.slice(state.cursorPos),
            cursorPos: state.cursorPos-1, 
          }})
        }
      }
      else if (key === "ArrowLeft") {
        if (this.state.cursorPos > 0) {
          // move cursor back
          write("\b")
          this.setState((state, props) => {return {
            cursorPos: state.cursorPos-1,
          }})
        }
      }
      else if (key === "ArrowRight") {
        if (this.state.cursorPos < this.state.input.length) {
          // write what is currently under cursor (moves cursor forward)
          write(this.state.input[this.state.cursorPos])
          this.setState((state, props) => {return {
            cursorPos: state.cursorPos+1,
          }})
        }
      }
      else if (printable) {
        this.setState((state,props) => {
          write(e.key + state.input.slice(state.cursorPos))
          write("\b".repeat(state.input.length-state.cursorPos))
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