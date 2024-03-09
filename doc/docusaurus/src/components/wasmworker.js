import {BlobWriter, ZipWriter, BlobReader} from "@zip.js/zip.js"
import wasm from "dune-copasi-wasm-git"
// import wasm from "dune-copasi-wasm-local"

function basename(path) {
    return path.split('/').reverse()[0];
}

function* recursiveDirWalk(instance, dir, ignore) {
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

// pass stdout back to main thread which updates terminal accordingly
const print = (text) => {
    postMessage({printText: text})
}

onmessage = (e) => {
    const {configText, dataFiles, dataFilesPath} = e.data
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
