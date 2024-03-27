import { useEffect, useRef, useState } from "react"
import BrowserOnly from '@docusaurus/BrowserOnly'

export default function WasmXterm() {
  return (<BrowserOnly fallback={<div>Loading terminal...</div>}>
    {() => {
      const XTerm = require("xterm-for-react").XTerm

      const xtermRef = useRef()
      const [input, setInput] = useState("")

      useEffect(() => {
        xtermRef.current.terminal.writeln("Please press any string then press enter:")
        xtermRef.current.terminal.write("echo> ")
      }, [])

      const onData = (data) => {
        const code = data.charCodeAt(0);
        // If the user hits empty and there is something typed echo it.
        if (code === 13 && input.length > 0) {
          xtermRef.current.terminal.write(
            "\r\nYou typed: '" + input + "'\r\n"
          );
          xtermRef.current.terminal.write("echo> ");
          setInput("")
        } else if (code < 32 || code === 127) { // Disable control Keys such as arrow keys
          return;
        } else { // Add general key press characters to the terminal
          xtermRef.current.terminal.write(data);
          setInput(input => input + data)
        }
      }

      return <XTerm
        ref={xtermRef}
        onData={onData}
      />
    }}
  </BrowserOnly>)
}