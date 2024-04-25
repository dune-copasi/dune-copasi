# DuneCopasi paper for Journal of Open Source Software

This folder contains the documents to generate the actual paper for JOSS

# Build the PDF

To build the PDF run the following ommand on this folder:

```bash
docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara
```

# Edit the document

Note that the document is formatted with an special flavor of Markdown for JOSS. A reference on the possible commands can be found here:

* [Markdown Format](https://www.markdownguide.org/)
* [JOSS Format](https://joss.readthedocs.io/en/latest/submitting.html#how-should-my-paper-be-formatted)
