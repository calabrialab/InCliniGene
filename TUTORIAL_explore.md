# Tutorial: Configure the database and run a query

## Neo4j Docker download

DockerHub hosts an **[official Neo4j image](https://hub.docker.com/_/neo4j)** that provides a standard, ready-to-run package of Neo4j Community Edition and Enterprise Edition for a variety of versions. In this study the Community Edition 4.4.21 has been used, therefore the this tutorial aims to reproduce the same environment.

```bash
sudo docker pull neo4j:4.4.21-community
sudo docker images

REPOSITORY   TAG                IMAGE ID       CREATED       SIZE
neo4j        4.4.21-community   bf473beffa62   13 days ago   551MB
```

The `docker run` command creates and starts a container. On the next line, `--name testneo4j` defines the name we want to use for the container as `testneo4j`. This avoids us having to reference the container by its generic id, making our container easier to reference and to remember.
More complete info about how to run Neo4j in Docker are available [**here**](https://neo4j.com/developer/docker-run-neo4j/).

```bash
sudo docker run \
--name testneo4j \
--restart always \
--publish=7474:7474 \
--publish=7687:7687 \
--env NEO4J_AUTH=neo4j/yourpassword \
bf473beffa62
```

Then we can interrupt this process with Ctrl+c, and simply re-starting the container created. 
It’s possible to start the active container with `docker start` . 
This allows us to continue working on the same bash session without having to open a new one.

```bash
sudo docker start testneo4j
```

## Neo4j configuration

The Neo4j query language is Cypher, the full manual is available at this [page](https://neo4j.com/docs/cypher-manual/current/introduction/).
In order to run Cypher directly in our container, we need to first access our container. We will need to use the command below in order to run any commands in a running container. In this case, we are telling docker to run bash within our container, allowing us to interact with our container using Linux bash commands. For a full list of options, check out **[Docker’s info](https://docs.docker.com/engine/reference/commandline/exec/)** on the `exec` command.

```bash
sudo docker exec -it testneo4j bash
```

Neo4j uses the following directories:

```bash
Directories in use:
home:         /home/neo4j
config:       /home/neo4j/conf
logs:         /home/neo4j/logs
plugins:      /home/neo4j/plugins
import:       /home/neo4j/import
data:         /home/neo4j/data
certificates: /home/neo4j/certificates
licenses:     /home/neo4j/licenses
run:          /home/neo4/run
```

## APOC download and configuration

APOC (Awesome Procedures on Cypher) is an add-on library for Neo4j that provides hundreds of procedures and functions adding a lot of useful functionality. More information on the APOC page on the neo4j [website](https://neo4j.com/labs/apoc/). 
It is required to run the function which import the database, which is a file in .graphml format, therefore we have to donwload it in the plugins directory of Neo4j.

```bash
cd plugins/
wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/4.4.0.18/apoc-4.4.0.18-all.jar
cd ..
```

In order to import the graph DB it is necessary to change some settings in the neo4j configuration file. 
Since nano (or any other text editor) is not installed in the docker, it is necessary to install it to edit the Neo4j configuration file. It is possible with the following two commands.

```bash
apt-get update
apt-get install nano
```

Then open the configuration file (neo4j.conf) with the following command, adding the line ‘apoc.import.file.enabled=true’ in the last row, as in the example below

```bash
nano conf/neo4j.conf
```

```bash
#********************************************************************
# Other Neo4j system properties
#********************************************************************

dbms.memory.pagecache.size=512M

dbms.default_listen_address=0.0.0.0
dbms.directories.logs=/logs

apoc.import.file.enabled=true
```

## Database download

The Large Scale Validation Assay used to benchmark the performances in our study is a file is > 6 Gb of data, thus it cannot be uploaded into this repository. Therefore, a smaller dataset with clonal expansion data has been uploaded on the Github repository of [InCliniGene](https://github.com/calabrialab/InCliniGene/raw/main/graph_db_export/clonalexp.graphml.gz). It can be downloaded into the import directory of neo4j-community-4.4.21 with the following commands.

```bash
cd neo4j/import
wget https://github.com/calabrialab/InCliniGene/raw/main/graph_db_export/clonalexp.graphml.gz
gunzip clonalexp.graphml.gz
cd ..
```

At this point, in order for the configurations you have made to take effect, you must restart Neo4j. You need to leave the container with **exit** or **ctrl-D.**
It’s possible to view the list of containers with the `docker ps` command, and then stop the active container with `docker stop`.

```bash
sudo docker ps

CONTAINER ID   IMAGE          COMMAND                  CREATED          STATUS          PORTS                                                                                            NAMES
4e2d5bcfbd69   bf473beffa62   "tini -g -- /startup…"   30 minutes ago   Up 30 minutes   0.0.0.0:7474->7474/tcp, :::7474->7474/tcp, 7473/tcp, 0.0.0.0:7687->7687/tcp, :::7687->7687/tcp   testneo4j

sudo docker stop testneo4j
```

Then, re-start the container and enter it with the command previously used

```bash
sudo docker start testneo4j
sudo docker exec -it testneo4j bash
```

## Cypher shell

We can now access Cypher shell by running the `cypher-shell` command, which is shown below. Notice that we also need to specify the username (`-u neo4j`) and password (`-p password`) in order to access the database, using the authentication values we set up when we created the container.

```bash
cypher-shell -u neo4j -p yourpassword

Connected to Neo4j using Bolt protocol version 4.4 at neo4j://localhost:7687 as user neo4j.
Type :help for a list of available commands or :exit to exit the shell.
Note that Cypher queries must end with a semicolon.
neo4j@neo4j>
```

Once inside the cypher-shell, you can query Neo4j. The first command to run is to import the clonalexp.graphml database using a function from the APOC library.

```bash
neo4j@neo4j> CALL apoc.import.graphml("clonalexp.graphml", {readLabels: true});
+-----------------------------------------------------------------------------------------------------------------------------------+
| file                | source | format    | nodes  | relationships | properties | time  | rows | batchSize | batches | done | data |
+-----------------------------------------------------------------------------------------------------------------------------------+
| "clonalexp.graphml" | "file" | "graphml" | 916630 | 1056294       | 4811655    | 34505 | 0    | -1        | 0       | TRUE | NULL |
+-----------------------------------------------------------------------------------------------------------------------------------+

1 row
```

## Query

At this point you can run queries against the database, the outputs will be printed directly below the query.
Now, for example, let's try querying the database with query #02 from the query table used in the benchmarks.

```bash
neo4j@neo4j> match (s:sample)--(u:subg)
             return s.DNAnumber, count(u);
+---------------------------+
| s.DNAnumber    | count(u) |
+---------------------------+
| "VA2020-mix01" | 311      |
| "VA2020-mix02" | 10134    |
| "VA2020-mix03" | 12410    |
| "VA2020-mix04" | 16721    |
| "VA2020-mix05" | 12938    |
| "VA2020-mix06" | 11848    |
| "VA2020-mix07" | 6003     |
| "VA2020-mix08" | 3494     |
| "VA2020-mix09" | 946      |
| "VA2020-mix10" | 659      |
| "VA2020-mix11" | 334      |
| "VA2020-mix12" | 179      |
| "VA2020-mix13" | 379      |
| "VA2020-mix14" | 1327     |
| "VA2020-mix15" | 1065     |
| "VA2020-mix16" | 2063     |
| "VA2020-mix17" | 5114     |
| "VA2020-CEM37" | 136      |
+---------------------------+

18 rows
```
