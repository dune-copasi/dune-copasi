import {BlobWriter, ZipWriter, BlobReader} from "@zip.js/zip.js"
import wasm_git from "dune-copasi-wasm-git"
// import wasm_local from "/dune-copasi.js"
import wasm_v2 from "dune-copasi-wasm-v2"

const versions = {
    "git": wasm_git,
    // "local": wasm_local,
    "v2": wasm_v2,
}

const home = "/dunecopasi"

// setup file system
const preRun = (instance) => {
    instance.FS.mkdir(home)
    instance.FS.chdir(home)
}

const print = (msg) => {
    postMessage({cmd: "print", msg: msg})
}

const wasm_options = {print, preRun}
var instance

function basename(path) {
    return path.split('/').reverse()[0];
}

const postError = (error) => {
    postMessage({cmd: "print", msg: `Error: ${error}`})
}

const invalidPath = (path) => !instance.FS.analyzePath(path).exists
const isDir = (path) => instance.FS.isDir(instance.FS.stat(path).mode)
const isFile = (path) => instance.FS.isFile(instance.FS.stat(path).mode)
const invalidParent = (path) => {
    const analysis = instance.FS.analyzePath(path)
    if (!analysis.parentExists)
        return false
    return !isDir(analysis.parentPath)
}
const listdir = (path) => instance.FS.readdir(path).filter(path => path !== "." && path !== "..")
function* recursiveDirWalk(dir) {
    for (const path of listdir(dir)) {
        const subpath = dir + "/" + path
        if (isDir(subpath))
            yield* recursiveDirWalk(subpath)
        else if (isFile(subpath))
            yield subpath
        else throw `path ${subpath} neither dir nor file`
    }
}
const zipPaths = async (paths) => {
    // create zip archive of paths
    var zipWriter = new ZipWriter(new BlobWriter())
    // add all paths to archive
    await Promise.all(
        paths.map(path =>
            zipWriter.add(path, new BlobReader(new Blob([instance.FS.readFile(path)])))
        )
    )
    // return completed archive
    return await zipWriter.close()
}

const helpMessage = `DuneCopasi UNIX-like shell
Commands:
help            - print this message
cd [DIR]        - change to home or given directory
ls [DIR]        - list current or given directory
mv PATH PATH    - move file into directory or rename file
cp PATH PATH    - copy file into directory or to a new path
mkdir DIR       - make a directory
rmdir DIR       - remove an empty directory
rm PATH         - remove file
edit PATH       - open file in editor
save PATH       - save current editor content to PATH
cat PATH        - print content of file
dune-copasi ... - run the dune-copasi executable, usage as on linux
download PATH   - download a single file directly or a directory as zip
upload          - upload files to the current directory
open PATH       - open file in the browser
`

var calls = {}

const getEditorText = () => {
    const id = Date.now()
    postMessage({cmd: "getEditorText", id})
    return new Promise((resolve, reject) => {calls[id] = {resolve, reject}})
}

const getUploadedFiles = () => {
    const id = Date.now()
    postMessage({cmd: "upload", id})
    return new Promise((resolve, reject) => {calls[id] = {resolve, reject}})
}

const commands = {
    help: (args) => print(helpMessage),
    echo: (args) => print(args.join(" ")),
    pwd: (args) => {
        print(instance.FS.cwd())
    },
    cd: (args) => {
        if (args.length > 1) {
            postError("usage: cd [DIR]")
            return
        }

        const path = args.length === 0 ? home : args[0]

        if (invalidPath(path) || !isDir(path)) {
            postError(`${path} is not a valid directory`)
            return
        }

        instance.FS.chdir(path)

        postMessage({cmd: "cwd", cwd: instance.FS.cwd()})
    }, 
    ls: (args) => {
        if (args.length > 1) {
            postError("usage: ls [DIR]")
            return
        }

        const path = args.length === 0 ? "." : args[0]

        if (invalidPath(path) || !isDir(path)) {
            postError(`${path} is not a valid directory`)
            return
        }

        postMessage({
            cmd: "print",
            msg: instance.FS.readdir(path).filter((dir) => dir !== "." && dir !== "..").join(" ")
        })
    },
    mv: (args) => {
        if (args.length < 2) {
            postError("usage: mv PATH PATH")
        }

        const move = (src, dst) => {
            const file = instance.FS.readFile(src, {encoding: "binary"})
            instance.FS.unlink(src)
            instance.FS.writeFile(dst, file, {encoding: "binary"})
        }

        if (args.length === 2) {
            let [arg1, arg2] = args
            let src
            let dst

            if (invalidPath(arg1)) {
                postError(`source ${src} is not a valid path`)
            }

            if (isFile(arg1) && !invalidParent(arg2)) {
                src = arg1
                dst = arg2
            } else if (isFile(arg1) && isDir(arg2)) {
                src = arg1
                if (!arg2.endsWith("/"))
                    arg2 += "/"
                dst = arg2 + basename(arg1)
            } else {
                postError("unimplemented usage")
            }
            move(src,dst)
        } else {
            postError("unimplemented usage")
        }
    },
    cp: (args) => {
        if (args.length < 2) {
            postError("usage: cp PATH... DIR | cp PATH PATH")
        }

        const copy = (src, dst) => {
            const file = instance.FS.readFile(src, {encoding: "binary"})
            instance.FS.writeFile(dst, file, {encoding: "binary"})
        }

        if (args.length === 2) {
            let [arg1, arg2] = args
            let src
            let dst

            if (invalidPath(arg1)) {
                postError(`source ${src} is not a valid path`)
            }

            if (isFile(arg1) && !invalidParent(arg2)) {
                src = arg1
                dst = arg2
            } else if (isFile(arg1) && isDir(arg2)) {
                src = arg1
                if (!arg2.endsWith("/"))
                    arg2 += "/"
                dst = arg2 + basename(arg1)
            } else {
                postError("unimplemented usage")
            }
            copy(src,dst)
        }
    },
    mkdir: (args) => {
        if (args.length !== 1) {
            postError("usage: mkdir PATH")
            return
        }

        const path = args[0]

        if (invalidParent(path)) {
            postError(`parent of ${path} is not a valid directory`)
            return
        }

        if (instance.FS.analyzePath(path).exists) {
            postError(`${path} exists already`)
            return
        }

        instance.FS.mkdir(path)
    },
    rmdir: (args) => {
        if (args.length !== 1) {
            postError("usage: rmdir DIR")
            return
        }

        const path = args[0]

        if (invalidPath(path) || !isDir(path)) {
            postError(`${path} is not a valid directory`)
            return
        }

        if (listdir(path).length !== 0) {
            postError(`${path} is not empty`)
            return
        }

        instance.FS.rmdir(path)
    },
    rm: (args) => {
        if (args.length !== 1) {
            postError("usage: rm FILE")
            return
        }

        const path = args[0]

        if (invalidPath(path) || !isFile(path)) {
            postError(`${path} is not a valid file`)
            return
        }

        instance.FS.unlink(path)
    },
    edit: (args) => {
        if (args.length !== 1) {
            postError("usage: edit FILE")
            return
        }

        const path = args[0]

        if (invalidPath(path) || !isFile(path)) {
            postError(`${path} is not a valid file`)
            return
        }

        const text = instance.FS.readFile(path, {encoding: "utf8"})
        postMessage({cmd: "edit", text})
    },
    save: async (args) => {
        if (args.length !== 1) {
            postError("usage: save FILE")
            return
        }

        const path = args[0]

        if (invalidParent(path)) {
            postError(`${path} is not a valid file`)
            return
        }

        const text = await getEditorText()
        instance.FS.writeFile(path, text)
    },
    cat: (args) => {
        if (args.length !== 1) {
            postError("usage: cat FILE")
            return
        }

        const path = args[0]

        if (invalidPath(path)) {
            postError(`${path} is not a valid file`)
            return
        }

        print(instance.FS.readFile(path, {encoding: "utf8"}))
    },
    "dune-copasi": (args) => {
        instance.callMain(args)
    },
    download: async (args) => {
        if (args.length !== 1) {
            postError("usage: download PATH")
            return
        }

        const path = args[0]

        if (invalidPath(path)) {
            postError(`${path} is not a valid file`)
            return
        }

        let filename
        let blob
        if (isFile(path)) {
            filename = basename(path)
            blob = new Blob([instance.FS.readFile(path, {encoding: "binary"})])
        } else if (isDir(path)) {
            filename = basename(path) + ".zip"
            const paths = Array.from(recursiveDirWalk(path))
            blob  = await zipPaths(paths)
        } else {
            postError(`cannot download ${path}`)
        }
        postMessage({cmd: "download", blob, filename})
    },
    upload: async (args) => {
        const files = await getUploadedFiles()
        for (let file of files) {
            const buf = new Uint8Array(await file.arrayBuffer())
            instance.FS.writeFile(file.name, buf, {encoding: "binary"})
        }
    },
    open: async (args) => {
        if (args.length !== 1) {
            postError("usage: open PATH")
            return
        }
        for (let path of args) {
            const buf = instance.FS.readFile(path, {encoding: "binary"})
            postMessage({cmd: "open", buf, name: basename(path)})
        }
    },
    instantiate: async (args) => {
        if (!(args.version in versions)) {
            postError(`version ${args.version} not in ${Object.keys(versions)}`)
            return
        }
        instance = await versions[args.version](wasm_options)
    
        instance.callMain(["--version"])

        for (const url in args.files) {
            const resp = await fetch(url)
            const blob = await resp.blob()
            const buf = new Uint8Array(await blob.arrayBuffer())
            instance.FS.writeFile(args.files[url], buf, {encoding: "binary"})
        }
    }
}

onmessage = async (e) => {
    if (e.data.cmd === "response") {
        if ("resolve" in e.data)
            calls[e.data.id].resolve(e.data.resolve)
        if ("reject" in e.data)
            calls[e.data.id].resolve(e.data.reject)
        return
    }

    const {cmd, args} = e.data
    if (cmd in commands) {
        try {
            await commands[cmd](args)
        } catch (error) {
            postError(`WasmWorker: unexpected js error: ${error}`)
        }
    } else {
        postError(`WasmWorker: unknown command '${cmd}'. Type 'help' for a list of commands.`)
    }
    postMessage({cmd: "exit"})
}