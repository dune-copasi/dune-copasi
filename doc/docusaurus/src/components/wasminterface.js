import {useRef} from "react"
import Editor from '@monaco-editor/react';

function basename(path) {
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

function Config({ setConfigText, setDataFiles, setDataFilesPath, examples }) {
    const editorRef = useRef(null)
    const dataFilesPathRef = useRef(null)

    const handleEditorDidMount = (editor, monaco) => {
        editorRef.current = editor
    }

    const fileChange = async e => {
        e.target.files[0].text().then(text => {
            editorRef.current.setValue(text)
            setConfigText(text)
        })
    } 

    const selectExample = async e => {
        // the selected example name
        const example = e.target.value

        // reset if selected empty example
        if (example == "") {
            editorRef.current.setValue("")
            setConfigText("")
            setDataFiles([])
            setDataFilesPath("")
            return
        }

        // fetch and set corresponding config ...
        fetch(examples[example].config)
        .catch(reason => alert(`error while fetch'ing example ${example}: ${reason}`))
        .then(resp => resp.text())
        .then(text => {
            editorRef.current.setValue(text)
            setConfigText(text)
        })

        // ... and data files
        Promise.all(
            examples[example].dataFiles
            .map(path =>
                fetch(path)
                .catch(reason => alert(`error while fetch'ing example ${example}: ${reason}`))
                .then(resp => resp.blob())
                .then(blob => new File([blob], basename(path)))
            )
        ).then(files => setDataFiles(files))

        const path = examples[example].dataFilesPath
        dataFilesPathRef.current.value = path
        setDataFilesPath(path)
    }

    return (<>
        <div style={{ display: "grid", gridTemplateColumns: 'auto auto', width: "fit-content", gap: "0.5em"}}>
            <label htmlFor="examples">Example:</label>
            <select name="examples"
                    onChange={selectExample}
                    style={{width: "fit-content"}}>
                <option value="">None</option>
                {Object.keys(examples).map(
                    example => <option value={example} key={example}>{example}</option>
                )}
            </select>
            <label htmlFor="configfile">Config:</label>
            <input type="file" id="configfile"
                    onChange={fileChange}/>
            <label htmlFor="datafiles">Data:</label>
            <input type="file" multiple id="datafiles"
                    onChange={e => setDataFiles(e.target.files)} />
            <label htmlFor="datafilespath">Data folder:</label>
            <input type="text" onChange={e => setDataFilesPath(e.target.value)}
                   ref={dataFilesPathRef}></input>
        </div>
        <Editor 
            onChange={(value,event) => setConfigText(value)}
            onMount={handleEditorDidMount}
            height="30vh"
            language="ini"
        />
    </>)
}

export default function WasmInterface({examples}) {
    const terminal = useRef(null)
    var configText = ""
    var dataFiles  = []
    var zipBlob = null
    var dataFilesPath = "."

    const updateTerminal = (content) => {
        terminal.current.value += content + "\n"
        terminal.current.scrollTop = terminal.current.scrollHeight;
    }

    var wasmWorker = {}
    if (typeof Worker !== 'undefined') {
        wasmWorker = new Worker(new URL("@site/src/components/wasmworker.js", import.meta.url))
    }

    const run = async () => {
        if (configText.trim() == "") {
            alert('Please input a config before computation!')
            return
        }

        // reset zipBlob --> not yet available
        zipBlob = null;
        updateTerminal("DuneCopasi |>\n")

        wasmWorker.postMessage({
            configText,
            dataFiles,
            dataFilesPath
        })
    }

    wasmWorker.onmessage = (e) => {
        // switch on different messages
        if ("printText" in e.data) {
            updateTerminal(e.data.printText)
        }
        if ("zipBlob" in e.data) {
            zipBlob = e.data.zipBlob
        }
        if ("error" in e.data) {
            alert("There was an unexpected error! See console log.")
            console.error(e.data.error)
        }
    }

    const handleDownload = () => {
        if (zipBlob == null) alert("Error: output files not (yet) available")
        else offerDownload("output.zip", zipBlob)
    }

    return (
        <form encType="multipart/form-data"
              style={{display: "flex", flexDirection: "column", gap: "1em"}}>
            <Config setConfigText={text => configText = text} 
                    setDataFiles={files => dataFiles = files}
                    setDataFilesPath={path => dataFilesPath = path}
                    examples={examples}
            />
            <button onClick={run} type="button"
                    style={{width: "fit-content"}}>
                Run
            </button>
            <div>
            <label>This is the stdout of dune-copasi.</label>
            <textarea ref={terminal} readOnly
                      style={{ height: "30vh", width: "100%" }} />
            </div>
            <button onClick={handleDownload} type="button" id="download"
                    style={{ width: "fit-content"}}>
                Download results
            </button>
        </form>
    )
}