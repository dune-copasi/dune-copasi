import React from 'react';
import MDXDetails from '@theme/MDXComponents/Details';
import MDXContent from '@theme/MDXContent';
import Admonition from '@theme/Admonition';
import Code from '@theme/MDXComponents/Code';
import ConfigOptsJSON from './assets/config_opts.json';
import Markdown from 'react-markdown';

function parseOptions(props, opts, prefix) {
    var content = [];
    for (const key in opts) {
        content.push(
            (
                <div key={prefix + "-" + key + "-div"}>
                    {opts[key][".brief"] && <Markdown>{"#### "+opts[key][".brief"]}</Markdown>}
                    <MDXDetails {...props} key={prefix + "-" + key + "-details"}>
                        <summary><Code>{key}</Code> </summary>
                        {opts[key][".type"] && <Admonition type="note" title="Type"><Code>{opts[key][".type"]} </Code></Admonition>}
                        {opts[key][".default"] && <Admonition type="tip" title="Default"><Code>{opts[key][".default"]}</Code></Admonition>}
                        {opts[key][".details"] && <Markdown>{ Array.isArray(opts[key][".details"]) ? opts[key][".details"].join('\n') : opts[key][".details"]}</Markdown> }
                        {parseOptions(props, opts[key][".options"], prefix + key)}
                    </MDXDetails>
                </div>
            )
        );
    }
    return (
        <div key={prefix + "-div"}>
            {content}
        </div>
    )
}


export default function ConfigOpts(props) {
    return (
        <MDXContent>
            <MDXDetails {...props} open key="cfg">
                <summary><b>Configuration Options</b></summary>
                {parseOptions(props, ConfigOptsJSON, "cfg")}
            </MDXDetails>
        </MDXContent>
    )

}