
# Chlamytina #
- New *Chlamydomonas reinhardtii* chromatin states
- Additional down-stream analysis for differential proteins
- Epigenome-browser 

## Workflow ##
<p align=center>
<img src=Timeline_Workflows/Workflow_Chlamytina.jpg />
</p>

## Vignette ##

### 0. Purpose ###

Chlamytina is a small project focused in the well known green-algae model *Chlamydomonas reinhardtii* trying to answer a common question in some proteomic/transcriptomic studies: \
\
                                                  **Are my molecules of interest epigenetically regulated?** \
\
To fill this gap, we collected all epigenectic files published until the date and developed a new **chromatin states model** including 6mA, 5mC and nucleosome-profile for the first time. Additionally,
an **epigenome-browser** was conducted focusing on the site-specific approach. This tool engage the link-up between proteomic/transcriptomic changes and epigenetic patterns, thus displaying the 
*Chlamydomonas reinhardtii* epi-proteogenomic landscape.             

### 1. Inputs ###

The data has to be input as a table (first column = Phytozome accessions, other columns = different treatment abundance/expression) in .xlsx, .xlx or tab-separated .txt format.

### 2. Installation ###

#### Via Github #### 

At least a Linux Ubuntu 19.04 OS is required. Download this repository as zip or clone with the command line

```
git clone https://github.com/RocesV/Chlamytina
```
You need a user with root permissions. Change directory to Chlamytina and run 0_ChlamytinaEssential.sh to install all dependecies and software needed.

```
cd Chlamytina/  # IMPORTANT: all processes need to be done inside Chlamytina directory NOT from Chlamytina subdirectories
bash Code/0_ChlamytinaEssential.sh
```
Using this via the following steps (1_DataPrepare and 2_EnrichmentsLOLA) will require more time (~ 1.5 h) for installing the R packages needed.  

#### Via Dockerhub (epigenome-browser) ####

Install Docker or Docker Desktop as follow https://docs.docker.com/ . This is the only way to get access to **epigenome-browser**. Pull rocesv/chlamytina image from dockerhub


```
docker pull rocesv/chlamytina
``` 

Build the container using the image pulled. Because the jbrowse inside the container is running in apache2 server, an empy port from the host (8080) need to be connected to container's 80 port. In order to 
share data between host and the container it is advisable to define a volume (-v) linking a host directory to /home/rocesv/Documents/Transfer folder. 

```
docker run -t -i -d --name chlamytina_rocesv -p 8080:80 -v <ABSOLUTE PATH TO HOST SHARED DIRECTORY>:/home/rocesv/Documents/Transfer rocesv/chamytina bash
```
Check docker container is running 

```
docker ps -a
```

Get inside the container

```
docker exec -i -t chlamytina_rocesv bash
```
Now you should see your container user like root@6a32e10fc951:/home#. Change directory to Chlamytina

```
cd Chlamytina/
```

Brief docker tutorial:

```
docker stop chlamytina_rocesv # Stop the container
docker start chlamytina_rocesv # Start the container
exit # Get outside the container 
```
Using this via the time required for the following steps is minimum.

### 3. Data Prepare ###

1_DataPrepare.R script performs differential expression analysis, unwanted variation correction, Phytozome Accession liftover to the last version (5.5), intersection between proteins of different treatments and bed
data generation (possible backgrounds and inputs for 4. LOLA Enrichment).   

```
Rscript --vanilla Code/1_DataPrepare.R -h

Usage: 1_DataPrepare.R [file] [condition] [file] [condition] ... [options]


Options:
	-A CHARACTER, --file1=CHARACTER
		Dataset1 file path. First column CreIDs. Other columns quantification data.

	--condition1=CHARACTER
		Dataset1 Condition vector. It representes replicates for each treatment, separated by - 
 	 Condition vector must contain all your replicates. For example (9 samples): 1) 3-6 will set a contrast between the first three replicates and the last six 
 	 2) 3-3-3 will set all possible two by two contrasts between the three treatments 


	-B CHARACTER, --file2=CHARACTER
		Dataset2 file path

	--condition2=CHARACTER
		Dataset2 Condition vector

	-C CHARACTER, --file3=CHARACTER
		Dataset3 file path

	--condition3=CHARACTER
		Dataset3 Condition vector

	-D CHARACTER, --file4=CHARACTER
		Dataset4 file path

	--condition4=CHARACTER
		Dataset4 Condition vector

	-E CHARACTER, --file5=CHARACTER
		Dataset5 file path

	--condition5=CHARACTER
		Dataset5 Condition vector

	-d DIFFERENTIAL, --differential=DIFFERENTIAL
		If true, differential expression limma based test is performed [default = TRUE]

	-s SVA, --sva=SVA
		If true, sva removing unwanted variation is performed. Only for n>10-15 samples datasets. [default = FALSE]

	-i CHARACTER, --intersect=CHARACTER
		CreIDs intra-inter group specific discrimination [default = TRUE]

	-o CHARACTER, --out=CHARACTER
		Output directory [default = ./Outputs/]

	-c CHARACTER, --chromosome=CHARACTER
		If true, non-chromosome mapped (scaffolds ...) proteins are not taked into account [default = TRUE]

	-n CHARACTER, --normalization=CHARACTER
		Normalization metric used. Options: normalizeQuantiles (limma), none 
 	 It is advisable to set this argument as none and preprocess the data with other pkgs like Processomics [default = none]

	-h, --help
		Show this help message and exit
```

Example:

```
Rscript --vanilla Code/1_DataPrepare.R -A <PATH TO DATASET>/Dataset1.xlsx --condition1 3-15 -B <PATH TO DATASET>/Dataset2.xlsx --condition2 4-20 -C <PATH TO DATASET>/Dataset3.xlsx --condition3 3-6
-s TRUE -n normalizeQuantiles
```

### 4. LOLA Enrichment ###

To decipher the potential epigenetic regulation of the dataset, 2_EnrichmentsLOLA.R script performs enrichment analysis based on genomic regions overlap using LOLA package and plots a heatmap. This analysis needs three components: \
\
a) Query set or input, as genomic regions (.bed outputs from 1_DataPrepare.R) \
b) Universe or background, set of regions that could potentially have been included in the query set. This depends on the biological question, see **FAQ** (.bed outputs from 1_DataPrepare.R or Data/DB/Universe.bed) \
c) region data base sets that are to be tested for overlap with the input (Data/regionDB/Chlamytina) 

```
Rscript --vanilla Code/2_EnrichmentsLOLA -h 

Usage: 2_EnrichmentsLOLA.R [file] [file] [file] [background] [database] ... [options]


Options:

	-A CHARACTER, --file1=CHARACTER

		First BED file path. Any BED file with chr, start and end or DataPrepare output

	-B CHARACTER, --file2=CHARACTER

		Second BED file path

	-C CHARACTER, --file3=CHARACTER

		Third BED file path

	-D CHARACTER, --file4=CHARACTER

		Fourth BED file path

	-E CHARACTER, --file5=CHARACTER

		Fifth file path

	-F CHARACTER, --file6=CHARACTER

		Sixth file path

	-G CHARACTER, --file7=CHARACTER

		Seventh file path

	-H CHARACTER, --file8=CHARACTER

		Eighth file path

	-I CHARACTER, --file9=CHARACTER

		Nineth file path

	-J CHARACTER, --file10=CHARACTER

		Tenth file path

	-b CHARACTER, --background=CHARACTER

		Background BED file path. The set of regions tested for enrichments

	-l CHARACTER, --list=CHARACTER

		If true, the rest of args are ignored and list all the files for one regionDB

	-o CHARACTER, --out=CHARACTER

		Output directory [default = ./Outputs/]

	-r CHARACTER, --database=CHARACTER

		regionDB used. Options: Marks (epigenetic marks by original conditions), MMarks (merged Marks wo conditions) or CS_Control, CS_N, CS_S (Ngan et al., Nat.Plants 2015, Chromatin States !Nitrogen !Sulfur) or CS_Chlamytina (Updated Chromatin states with 5mC, 6mA and MNase) [default = MMarks]

	-c CHARACTER, --cores=CHARACTER

		Number of cores [default = 1]

	-h, --help

		Show this help message and exit
```
Example:

```
Rscript --vanilla Code/2_EnrichmentsLOLA.R -A Query1.bed -B Query2.bed -C Query3.bed -D Query4.bed -b Data/DB/Universe.bed -r CS_Chlamytina
```

### 5. Epigenome browser ###

Once you got into the docker container (2. Installation - Via Dockerhub) you need to start the apache2 server 
```
service apache2 start
```
Now you can enjoy the epigenome browser at http://localhost:8080/jbrowse . The browser will be available while the container is running so as long as ```docker stop chlamytina_rocesv``` is not executed you 
can acces to jbrowse. We recommend to always select refseq track and one of the .Genes tracks (Nuclear, Mitochondrion, Chloroplast). The epigenomic tracks can be displayed by condition (Control, light ...)
or merged (M-). 

### 6. FAQ ###

**(Q) What type of outputs can I obtain from Chlamytina?**

- It depends on the arguments selected. For a more detailed description, see Outputs/README.md 

**(Q) What type of inputs can I use?**

- Both transcripts and proteins can be used as long as the accessions comes from Phytozome

**(Q) Which background/universe should I use?**

- Using the same subset of proteins/transcripts as query you can obtain different results depending on the universe/background selected. Enrichments are computed using the background/universe as
reference so all potential enrichments already present in your background will be deprecated. To obtain a general view of your data we recommend use first
whole proteome/coding-transcriptome background (Data/DB/Universe.bed) because your global background may have some enrichments. Possible background in Chlamytina: Universe.bed, Global_background.bed 
(total set of your proteins), Differential_background.bed (total set of your differential proteins), file1.bed (single file set of proteins). Any .bed file that include your query can be used as background
but you have to be careful with query:background size.      

**(Q) Why are there different regionDBs?**

- To be able to elucidate the epigenetic regulation at different levels (marks with conditions, merged marks, chromatin states) we decided to maintain several regionDBs so you can choose the one that
best suits your questions. Data/DB/CS_Chlamytina_interpretation heatmaps would make easier the biological interpretation of CS_Chlamytina regionDB results  

