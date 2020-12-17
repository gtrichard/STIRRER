# STIRRER
STIRRER Project Code Repository


## Project Description

### Biological material

### Available Datasets

#### RNA-seq

#### small RNA-seq

#### BS-seq

#### Hi-C

## How to run the analyses of this project using this git repository

### Snakepipes installation to your favorite cluster

#### SlurmEasy installation

```
mkdir ~/bin
wget https://github.com/dpryan79/Misc/blob/master/MPIIE_internal/SlurmEasy -O ~/bin/SlurmEasy
export PATH="~/bin/:$PATH"
```
Now we need to modify some SlurmEasy parameters so it actually works on your cluster.

Line 23 you need to change the value of QUEUE to the name of the default slurm queue used by your cluster.
Line 29 you need to change the value of QUEUEs to the name of the available queues in your cluster.
Line 31 you need to change the value of MAINTEMP to the directory your cluster uses to store temproray files.

#### Miniconda installation

Install the latest miniconda version following these steps on your cluster:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

You might need to export the miniconda location to your PATH
```
export PATH="/path/to/miniconda3/bin:$PATH"
```

#### SnakePipes installation

```
conda create -n snakePipes2 -c mpi-ie -c conda-forge -c bioconda snakePipes==2.4.1
```

To verify that the installation worked:
```
conda activate snakePipes2
snakemake --help
snakePipes --help
```



#### Creating the necessary conda environments for SnakePipes

```
snakePipes createEnvs
```

### Other needed conda environments installation

Bedtools

### Running the data analyses using the RMarkdown notebooks

#### My cluster has a RStudio server running which has access to the cluster files

Consider you lucky, this will be easy. Simply clone this repository on your cluster and open the RMarkdown (.Rmd) files in each folder to run the pipelines and downstream analyses.

```
git clone https://github.com/gtrichard/STIRRER
```

#### I run RStudio locally and perform computationnaly intense tasks on a cluster

If you have a Linux machine locally, skip the Windows step.

##### I am unlucky and have a Windows machine for bioinformatics

You need to install WSL and follow this tutorial to get an RStudio server working on your linux subsystem (so it has access to sshfs):

https://jmgirard.com/rstudio-wsl2/

##### Clone the repository

Clone this git repo on your local machine and on the cluster where you want to perform the analyses (the cloning on your machine is only useful to get the Rmd files).

```
git clone https://github.com/gtrichard/STIRRER
```

##### Mounting your cluster filesystem on your local machine

Install sshfs and create a directory where you want to mount the cluster filesystem and mount it.
```
apt-get update && apt-get install sshfs
sudo mkdir /mnt/cluster
sudo sshfs -o allow_other,IdentityFile=/path/to/my/id_rsa name@cluster.host.address:/ /mnt/cluster
```

This last sudo line must be launched everytime you reboot your system. It might be convenient to put it as an alias in your ~/.bashrc

#### 



