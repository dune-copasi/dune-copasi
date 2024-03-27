import {BlobWriter, ZipWriter, BlobReader} from "@zip.js/zip.js"
import wasm from "dune-copasi-wasm-git"
// import wasm from "dune-copasi-wasm-local"

const home = "/dunecopasi"

// setup file system
const preRun = (instance) => {
    instance.FS.mkdir(home)
    instance.FS.chdir(home)
}

const print = (msg) => {
    postMessage({cmd: "print", msg: msg})
}

const instance = await wasm({print,preRun})

function basename(path) {
    return path.split('/').reverse()[0];
}

function* recursiveDirWalk(dir, ignore) {
    for (const path of instance.FS.readdir(dir)) {
        const mode = instance.FS.stat(dir + "/" + path)["mode"]
        if (path == "." || path == ".." || ignore.includes(basename(path)))
            continue
        if (instance.FS.isDir(mode))
            yield* recursiveDirWalk(instance, dir + "/" + path, ignore)
        else if (instance.FS.isFile(mode))
            yield dir + "/" + path
        else throw `path ${dir}/${path} neither dir nor file`
    }
}

const postError = (error) => {
    postMessage({cmd: "print", msg: `Error: ${error}`})
    postMessage({cmd: "exit"})
}

const invalidPath = (path) => !instance.FS.analyzePath(path).exists
const invalidParent = (path) => !instance.FS.analyzePath(path).parentExists
const isDir = (path) => instance.FS.isDir(instance.FS.stat(path).mode)

const helpMessage = `DuneCopasi shell: limited bash like UNIX shell
help     - print this message
cd [DIR] - change to home or given directory
ls [DIR] - list current or given directory
`

const commands = {
    help: (args) => print(helpMessage.replace(/\n/g, "\r\n")),
    echo: (args) => print(args.join(" ")),
    pwd: (args) => {
        print(instance.FS.cwd())
    },
    cd: (args) => {
        if (args.length > 1) {
            postError({cmd: "print", msg: "usage: cd [PATH]"})
            return
        }

        const path = args.length === 0 ? home : args[0]

        console.log()
        if (invalidPath(path) || !isDir(path)) {
            postError(`${path} is not a valid directory`)
            return
        }

        instance.FS.chdir(path)

        postMessage({cmd: "cwd", cwd: path})
    }, 
    ls: (args) => {
        if (args.length > 1) {
            postError({cmd: "print", msg: "usage: ls [PATH]"})
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
}

onmessage = (e) => {
    console.log("wasmworker received:", e.data)
    const {cmd, args} = e.data
    if (cmd in commands)
        commands[cmd](args)
    else
        postError(`wasmworker: ${cmd} not implemented`)
    postMessage({cmd: "exit"})
}

/*
onmessage = (e) => {
    wasm({print}).then(async instance => {
        // create and switch to working directory
        instance.FS.mkdir("/dune-copasi")
        instance.FS.chdir("/dune-copasi")

        // write config file to virtual FS
        instance.FS.writeFile('config.ini', configText)

        // write data files
        if (dataFiles.length !== 0) {
            const length = dataFilesPath.length
            if (dataFilesPath[length-1] == "/")
                dataFilesPath = dataFilesPath.slice(0, -1)

            if (dataFilesPath == "")
                dataFilesPath = "."

            const folders = dataFilesPath.split("/")
            var path = "."
            for (const folder of folders) {
                if (folder == "")
                    continue
                path = path + "/" + folder
                instance.FS.mkdir(path)
            }

            for (const dataFile of dataFiles) {
                const data = new Uint8Array(await dataFile.arrayBuffer())
                instance.FS.writeFile(dataFilesPath + "/" + dataFile.name, data)
            }
        }

        // call main with hardcoded paths
        instance.callMain(["--config=config.ini", "--model.writer.vtk.path=output"])

        print("\n")

        // gather all output file paths
        const ignore = [...Array.from(dataFiles, file => file.name), "config.ini"]
        const paths = Array.from(recursiveDirWalk(instance, ".", ignore))

        // create zip archive of paths
        var zipWriter = new ZipWriter(new BlobWriter())
        Promise.all(
            paths.map(path =>
                zipWriter.add(path, new BlobReader(new Blob([instance.FS.readFile(path)])))
            )
        ).then(() => {
            zipWriter.close().then(zipBlob => postMessage({zipBlob}))
        })
    }).catch((reason) => {
        postMessage({error: reason})
    })
}
*/