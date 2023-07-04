# Utils

Here we describe some utils and how to process data with InCliniGene. Please follow the [configuration tutorial and demo](TUTORIAL_explore.md) for database configuration and query.
All files here referred, both Java files and their sample data, are placed in the folders [utils](utils) folder, except for the query script located in the custom folder [query_graphdb](query_graphdb).

## Step 1: metadata parsing

In order to process the data of a sparse matrix, you need to add and parse sample metadata through a file, formatted as tab-separated file format (TSV). An example of the file is provided here: [sample_metadata](utils/sample_metadata_file.tsv). 
The required columns are:
- "CompleteAmplificationID": the ID of the sample (corresponding to a column in the sparse matrix of IS.
- "UniqueID": an extra ID for the sample (both columns can have the same name/ID).

To import the metadata file, you need to execute the following command line [MetadataParser.java](utils/MetadataParser.java) with the input file [sample_metadata](utils/sample_metadata_file.tsv):

```bash
javac MetadataParser.java
java MetadataParser sample_metadata_file.tsv
```

NB: before running the command, please change UID and password.

## Step 2: sparse-matrix parsing

To convert a general-purpose sparse matrix to a Jason file format, to be later processed for data import, you need need to run the file converter [SparseMatrix2Jason.java](utils/SparseMatrix2Jason.java) on an input file such as the sample matrix ([sample_sparse_ISMatrix.tsv](utils/sample_sparse_ISMatrix.tsv)) with the following statements:

```
javac SparseMatrix2Jason.java
java SparseMatrix2Jason sample_sparse_ISMatrix.tsv
```

The script will return a file (with extension _.clu_) for each column, indeed a file "cluster" per sample.

If you already have the output from gTRIS, you can skip this step.

## Step 3: cluster file parsing and database import

Having all the cluster files ready for the import phase, we now need to parse each ".clu" file. An example is with the following statement (NB: update the code with your UID and password) using the Java script [CluParser.paper.java](import_data_graphdb/CluParser.paper.java):

```
javac CluParser.paper.java
for k in $( ls *.clu ); do
  echo "Parsing file ${k}"
  java CluParser.paper ${k}
done
```

This statement will import all data in the database, crating all the links and nodes.


## Step 4: database export into a sparse matrix file (TSV)

To export all desired data stored in the graph database, we will need to run the Java script [TableCreator.paper.java](query_graphdb/TableCreator.paper.java) with the explicit table name and the read count as quantification of each IS (here expressed with the name _clu_) as follows:

```
javac TableCreator.java
java TableCreator table clu ## 
```

This statement will recreate a sparse matrix file with IS in rows, sample in columns, and read count as the values. 
