import BrowserOnly from '@docusaurus/BrowserOnly'
import React from "react"

const ASCII = {
  ENTER: 13,
  DELETE: 127,
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
      cwd: this.home
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

    const onData = (data) => {
      const code = data.charCodeAt(0);

      // handle enter
      if (code === ASCII.ENTER && this.state.input.length > 0) {
        this.xtermRef.current.terminal.write(
          "\r\nYou typed: '" + this.state.input + "'\r\n"
        )
        this.xtermRef.current.terminal.write("echo> ")
        this.setState({ input: "" })
        return
      } 
      // ignore other control characters
      else if (code < 32) { 
        console.error("WasmXterm: unhandled code", code)
        return
      } 
      // handle delete
      else if (code === ASCII.DELETE) {
        if (this.state.input.length > 0) {
          this.xtermRef.current.terminal.write("\b \b")
          this.setState({input: this.state.input.slice(0,-1)})
        }
      }
      // handle regular printable char
      else if (code < 128) {
        this.xtermRef.current.terminal.write(data)
        this.setState({
          input: this.state.input + data,
        })
      }
    }

    return (
      <XTerm
        ref={this.xtermRef}
        onData={onData}
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