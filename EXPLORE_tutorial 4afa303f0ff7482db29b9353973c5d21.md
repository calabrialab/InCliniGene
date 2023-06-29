# EXPLORE_tutorial

## Neo4j download

Download from the [Neo4j Download Center](https://neo4j.com/download-center/?utm_medium=cpc&utm_source=google&utm_campaign=emea-search-offers&utm_adgroup=desktop-download&utm_content=desktop-download&gclid=Cj0KCQjwyOuYBhCGARIsAIdGQRPdf7kUzGDg120L9eEsdj8U8rrrFIvwXfKpui8aq5Pnr0PH8SEbJokaApkiEALw_wcB#community) the Neo4j community edition in *.tar* format:
Currently, the latest version is the 4.4.21

```bash
$ wget https://neo4j.com/artifact.php?name=neo4j-community-4.4.21-unix.tar.gz
```

Then open the .tar package with the following command

```bash
$ tar -xf artifact.php?name=neo4j-community-4.4.21-unix.tar.gz
```

a folder named neo4j-community-4.4.21 should be created. Within that folder it is possible to start Neo4j.

## APOC download and configuration

APOC (Awesome Procedures on Cypher) is an add-on library for Neo4j that provides hundreds of procedures and functions adding a lot of useful functionality. More information on the APOC page on the neo4j [website](https://neo4j.com/labs/apoc/). It is required to run the function which import the database, which is a file in .graphml format.

```bash
$ cd neo4j-community-4.4.21/plugins/
$ wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/4.4.0.18/apoc-4.4.0.18-all.jar

```

In order to import the graph DB it is necessary to change some settings in the neo4j configuration file, adding line ‘apoc.import.file.enabled=true’ in the last row.

```bash
$ nano neo4j-community-4.4.21/conf/neo4j.conf
```

```bash
#********************************************************************
# Other Neo4j system properties
#********************************************************************
apoc.import.file.enabled=true
```

## Database download

The clonal expansion database is available on the Github repository of [InCliniGene](https://github.com/calabrialab/InCliniGene/raw/main/graph_db_export/clonalexp.graphml.gz). It can be downloaded int the import folder of neo4j-community-4.4.21 with the following commands.

```bash
$ cd neo4j-community-4.4.21/import
$ wget https://github.com/calabrialab/InCliniGene/raw/main/graph_db_export/clonalexp.graphml.gz
$ gunzip clonalexp.graphml.gz
```

## Neo4j starting

Finally Neo4j can be started and the clonal expansion database can be imported. 
First the database is started, this process can take several seconds, even after the output message has been printed. The query language for neo4j is Cypher, the full manual is available at this [page](https://neo4j.com/docs/cypher-manual/current/introduction/).
After that you can open the cypher-shell using the initial default credentials:
username: neo4j
password: neo4j
You will be asked to set a new password that you will use from now on to access the shell.

```bash
$ cd neo4j-community-4.4.21
$ bin/neo4j start

Directories in use:
home:         /home/neo4j-community-4.4.21
config:       /home/neo4j-community-4.4.21/conf
logs:         /home/neo4j-community-4.4.21/logs
plugins:      /home/neo4j-community-4.4.21/plugins
import:       /home/neo4j-community-4.4.21/import
data:         /home/neo4j-community-4.4.21/data
certificates: /home/neo4j-community-4.4.21/certificates
licenses:     /home/neo4j-community-4.4.21/licenses
run:          /home/neo4j-community-4.4.21/run
Starting Neo4j.
WARNING: Max 4096 open files allowed, minimum of 40000 recommended. See the Neo4j manual.
Started neo4j (pid:22518). It is available at http://localhost:7474
There may be a short delay until the server is ready.

$ bin/cypher-shell -u neo4j -p neo4j
Password change required
new password: 
confirm password: 
Connected to Neo4j using Bolt protocol version 4.4 at neo4j://localhost:7687 as user neo4j.
Type :help for a list of available commands or :exit to exit the shell.
Note that Cypher queries must end with a semicolon.
**neo4j@neo4j>**
```

Once inside the cypher-shell, you can query neo4j. The first command to run is to import the clonalexp.graphml database using a function from the APOC library.

```bash
**neo4j@neo4j>** CALL apoc.import.graphml("clonalexp.graphml", {readLabels: true});
+-----------------------------------------------------------------------------------------------------------------------------------+
| file                | source | format    | nodes  | relationships | properties | time  | rows | batchSize | batches | done | data |
+-----------------------------------------------------------------------------------------------------------------------------------+
| "clonalexp.graphml" | "file" | "graphml" | 916630 | 1056294       | 4811655    | 34505 | 0    | -1        | 0       | TRUE | NULL |
+-----------------------------------------------------------------------------------------------------------------------------------+

1 row
```