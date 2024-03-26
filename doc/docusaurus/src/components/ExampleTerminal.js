import React from "react";
import { ReactTerminal, TerminalContextProvider } from "react-terminal";

// Simple window to render of terminal instructions and its output in the terminal
export default function ExampleTerminal({ children, prompt, title }) {
  prompt = prompt ? prompt : "~$";
  const theme = "custom";
  const themes = {
    custom: {
      themeBGColor: "rgb(41, 45, 62)",
      themeToolbarColor: "#ffffff1a",
      themeColor: "rgb(191, 199, 213)",
      themePromptColor: "rgb(130, 170, 255)",
    },
  };

  const makePrompt = (prompt) => {
    return (
      <span
        style={{
          color: themes[theme].themePromptColor,
          WebkitUserSelect: "none",
          KhtmlUserSelect: "none",
          MozUserSelect: "none",
          msUserSelect: "none",
          OUserSelect: "none",
          userSelect: "none",
        }}
      >
        {prompt}
      </span>
    );
  };

  const welcomeMessage = (
    <span style={{ whiteSpace: "pre-wrap" }}>
      {children.map((item, index) => (
        <span key={index}>
          {makePrompt(item.prompt ? item.prompt : prompt)}{" "}
          {/* handle input sequences */}
          <span>
            {Array.isArray(item.input) ? (
              item.input.map((result, innerIndex) => (
                <span key={innerIndex}>
                  {result} <br />
                </span>
              ))
            ) : (
              <span>
                {" "}
                {item.input} <br />
              </span>
            )}
          </span>
          {/* handle output sequences */}
          {item.output &&
            item.output.map((result, innerIndex) => (
              <span key={innerIndex}>
                {result} <br />
              </span>
            ))}
        </span>
      ))}
    </span>
  );

  const buttonStyle = {
    display: "inline-block",
    width: "10px",
    float: "left",
    height: "10px",
    background: "#f9f9f9",
    borderRadius: "50%",
    margin: "8px 4px 0 0",
  };

  return (
    <>
      <div
        style={{
          width: "100%",
          borderRadius: "7px 7px 7px 7px",
          WebkitBoxShadow: "rgba(0, 0, 0, 0.3) 0 1px 3px",
          MozBoxShadow: "rgba(0,0,0,0.3) 0 1px 3px",
          BoxShadow: "rgba(0, 0, 0, 0.3) 0 1px 3px",
          WebkitFontSmoothing: "antialiased",
          textRendering: "optimizelegibility",
        }}
      >
        <div
          style={{
            width: "100%",
            position: "relative",
            WebkitBoxSizing: "border-box",
            borderRadius: "7px 7px 7px 7px",
          }}
        >
          <div
            style={{
              background: "rgb(244,244,239)",
              width: "100%",
              borderRadius: "7px 7px 7px 7px",
            }}
          >
            <div
              style={{
                height: "25px",
                padding: "0px 10px",
                color: "#4d494d",
                textAlign: "center",
                fontFamily:
                  "HelveticaNeue, 'Helvetica Neue', 'Lucida Grande', Arial, sans-serif",
              }}
            >
              <div
                style={{
                  ...buttonStyle,
                  background: "#FF6057",
                  border: "1px solid #E14640",
                }}
              />
              <div
                style={{
                  ...buttonStyle,
                  background: "#FFBD2E",
                  border: "1px solid #DFA123",
                }}
              />
              <div
                style={{
                  ...buttonStyle,
                  background: "rgba(218, 218, 213)",
                  border: "1px solid rgba(203, 203, 196)",
                }}
              />
              <div
                style={{
                  height: "10px",
                  textAlign: "center",
                  fontSize: 14,
                  opacity: 0.75,
                  WebkitFontSmoothing: "antialiased",
                  textRendering: "optimizelegibility",
                }}
              >
                {title && <>{title}</>}
              </div>
            </div>
            <div
              style={{
                width: "100%",
                background: "rgb(41,45,62)",
                margin: "0px auto",
                borderRadius: "0px 0px 7px 7px",
                padding: "7px 7px 7px",
                position: "relative",
                WebkitBoxSizing: "border-box",
                boxSizing: "border-box",
              }}
            >
              <TerminalContextProvider>
                <ReactTerminal
                  prompt={prompt}
                  theme={"custom"}
                  themes={themes}
                  showControlButtons={false}
                  welcomeMessage={welcomeMessage}
                  enableInput={false}
                  showControlBar={false}
                />
              </TerminalContextProvider>
            </div>
          </div>
        </div>
      </div>
      <br />
    </>
  );
}
