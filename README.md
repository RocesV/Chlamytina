
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

Chlamytina is a small project focused in the well known green-algae model *Chlamydomonas reinhardtii* that tries to answer a common questions in some proteomic/transcriptomic studies: \
\
**Are my molecules of interest epigenetically regulated?** \
\
To fill this gap, we collected all epigenectic files published until the date and developed a new **chromatin states model** including 6mA, 5mC and nucleosome-profile for the first time. Additionally,
an **epigenome-browser** was conducted focusing on the site-specific approach. This tool engage the link-up between proteomic/transcriptomic changes and epigenetic patterns, thus displaying the 
*Chlamydomonas reinhardtii* epi-proteogenomic/epi-transcriptomic landscape.             

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
cd Chlamytina/
bash Code/0_ChlamytinaEssential.sh # IMPORTANT: all processes need to be done inside Chlamytina directory NOT from Chlamytina subdirectories
```
Using this via the following steps (1_DataPrepare and 2_EnrichmentsLOLA) will require more time (~ 1.5h) for installing the R packages needed.  

#### Via Dockerhub (epigenome-browser) ####

Install Docker or Docker Desktop as follow https://docs.docker.com/ . This is the only way to get access to **epigenome-browser**. Pull rocesv/chlamytina image from dockerhub


```
docker pull rocesv/chlamytina
``` 

Build the container using the image pulled. Because the jbrowse inside the container is running in apache2 server, an empy port from the host (8080) need to be connected to container's 80 port. In order to 
share data between host and the container it is recommended to define a volume (-v) linking a host directory to /home/rocesv/Documents/Transfer folder. 

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

 


