# Utils

Here we describe some utils and how to process data with InCliniGene. Please follow the [configuration tutorial and demo](TUTORIAL_explore.md) for database configuration and query.
All files here referred, both Java files and their sample data, are placed into the [utils](utils) folder, except for the query script located in the custom folder [query_graphdb](query_graphdb).

## Step 1: metadata parsing

In order to process the data of a sparse matrix, you need to add and parse sample metadata through a file, formatted as tab-separated file format (TSV) and with the following columns: 

ADD REQUIRED COLUMNS

To import the metadata file, you need to execute the following parser.

```bash
javac MetadataParser.java
java MetadataParser sample_metadata_file.tsv
```

EXPECTED OUTPUT 

## Step 2: metadata parsing

To convert a general-purpose sparse matrix to a Jason file format, to be later processed for data import, you need need to run the file converter _sparseMatrix2Jason.java_ with the following statements:

```
javac SparseMatrix2Jason.java
java SparseMatrix2Jason sample_sparse_ISMatrix.tsv
```

The script will return a file (with extension _.clu_) for each column, indeed a file per sample.

## Step 3: database import

Once moved in the folder [query_graphdb](query_graphdb), run the following command:

```
javac TableCreator.java
java TableCreator table sample_file_clu ## sample_file_clu is a file sample name just created.
```

This statement will import all data in the specified table of the graph DB.
