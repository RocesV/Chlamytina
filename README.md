
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
 

#### Via Dockerhub (epigenome-browser) ####
