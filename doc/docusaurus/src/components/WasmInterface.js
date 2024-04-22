import {useRef} from "react"
import Editor from '@monaco-editor/react';
import WasmXterm from "./WasmXterm"

export default function WasmInterface({ version, examples }) {
    const editorRef = useRef(null)

    const handleEditorDidMount = (editor, monaco) => {
        editorRef.current = editor
    }

    return (<div style={{display: "flex", flexDirection: "column", gap: "1em"}}>
        <Editor
            onMount={handleEditorDidMount}
            height="30vh"
            language="ini"
        />
        <WasmXterm 
            version={version}
            examples={examples}
            setEditorText={text => editorRef.current.setValue(text)}
            getEditorText={() => editorRef.current.getValue()}
        />
    </div>)
}