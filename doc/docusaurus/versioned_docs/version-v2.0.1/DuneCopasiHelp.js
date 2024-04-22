import ExampleTerminal from '@site/src/components/ExampleTerminal'

export default function DuneCopasiHelp({prefix}) {
  prefix = prefix ? prefix : "";
  return (
    <ExampleTerminal title="Command Line Interface of DuneCopasi">
      {[
        {
          input: `${prefix}dune-copasi --help`,
          output: [
            "USAGE: dune-copasi [options]",
            "",
            "Numerical simulator for diffusion-reaction systems on single or multiple compartments",
            "",
            "Website: <https://dune-copasi.netlify.app>",
            "",
            "OPTIONS:",
            "",
            "  -h / --help          - Display this help.",
            "  --help-full          - Display this help with long descriptions.",
            "  --version            - Display the version of this program.",
            "  --parser-default     - Display the default parser.",
            "  --parser-list        - Display the parsers available for this program.",
            "  --dimension-list     - Display the grid dimensions available for this program.",
            "  --dump-config        - Dumps configuration in the INI format to stdout.",
            "  --config=<string>    - Specifies a config file in INI format. See Configuration Options.",
            "  --{key}={value}      - Overrides key=value sections of the config file. See Configuration Options.",
            "",
            "...",
          ],
        },
      ]}
    </ExampleTerminal>
  );
}
