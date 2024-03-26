import { useState } from "react";
import { ReactTerminal, TerminalContextProvider } from "react-terminal";

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

  const [cwd, setCwd] = useState("/dunecopasi")

  // Define commands here
  const commands = {
    whoami: "dunecopasi",
    cd: (args) => {
      if (args.split(" ").length > 1) 
        return "Error: usage: cd <path>"
      const path = args
      setCwd(path.startsWith("/") ? path : cwd + "/" + path)
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
  );
}