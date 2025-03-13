#!/bin/bash
#SBATCH --job-name=rscript_cart # Job name
#SBATCH -o slurm.%j.out    # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err    # STDERR (%j = JobId)
#SBATCH --mail-type=FAIL   # notifications for job done & fail
##SBATCH --mail-user=hnatri@tgen.org # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G
##SBATCH --part=hmem
#SBATCH -t 24:00:00

export SIMG_FILE_NAME=rstudio-4.3.0-4-with_modules.sif

#check if there is a sif file specified on the sbatch, if not use default
#if [[ -z "${SIMG_FILE_NAME}" ]]; then
#  export SIMG_FILE_NAME=Rstudio-WithModules-WithJuypter.4.2.0.sif
#else
#  echo -e "The image we are using is $SIMG_FILE_NAME\n"
#fi


export RPRIHOME=$HOME/R/$SIMG_FILE_NAME/
#if [ -d /scratch ]
#then
#  export RScratch=/scratch/$USER/R/$SIMG_FILE_NAME/
#else
#  export RScratch=$HOME/R/$SIMG_FILE_NAME/
#fi

export RScratch=$HOME/R/$SIMG_FILE_NAME/

if [ ! -d ${RScratch}/tmp ]
then
  mkdir -p ${RScratch}/tmp
fi

if [ ! -d ${RPRIHOME} ]
then
  mkdir -p ${RPRIHOME}
fi

if [ ! -d ${HOME}/R/$SIMG_FILE_NAME/libs ]
then
  mkdir -p $HOME/R/$SIMG_FILE_NAME/libs
fi

if [ ! -d ${HOME}/R/$SIMG_FILE_NAME/.rstudio ]
then
  mkdir -p $HOME/R/$SIMG_FILE_NAME/.rstudio
fi

if [ ! -d ${HOME}/R/$SIMG_FILE_NAME/.config ]
then
  mkdir -p $HOME/R/$SIMG_FILE_NAME/.config
fi

# User-installed R packages go into their home directory
if [ ! -e ${RPRIHOME}/.Renviron ]
then
  printf '\nNOTE: creating ~/.Renviron file'
  echo "R_LIBS_USER=$HOME/R/$SIMG_FILE_NAME/libs" > ${RPRIHOME}/.Renviron
fi

if [ ! -e ${RPRIHOME}/.Rprofile ]
then
  printf '\nNOTE: creating ~/.Rprofile file\n\n'
  echo "myPaths <- .libPaths()" > ${RPRIHOME}/.Rprofile
  echo "myPaths <- c(myPaths[3], myPaths[2],myPaths[1])" >> ${RPRIHOME}/.Rprofile
  echo ".libPaths(myPaths)" >> ${RPRIHOME}/.Rprofile
fi

mkdir -p $RScratch/rstudio-server/lib
mkdir -p $RScratch/rstudio-server/run
mkdir -p $RScratch/rstudio-server/tmp
export RSLIB=$RScratch/rstudio-server/lib
export RSRUN=$RScratch/rstudio-server/run
export RSTMP=$RScratch/rstudio-server/tmp

# This example bind mounts the /project directory on the host into the Singularity container.
# By default the only host file systems mounted within the container are $HOME, /tmp, /proc, /sys, and /dev.

#change this before deploying this to point at our repo
#this was pulled using singularity build singularity-rstudio4.R4.2.simg docker://rocker/rstudio

if [ -d /packages/containers/RStudio/ ]
then
  export SIMG_IMAGE=/packages/containers/RStudio/$SIMG_FILE_NAME
elif [ -d /packages/containers/Rstudio ]
then
  export SIMG_IMAGE=/packages/containers/Rstudio/$SIMG_FILE_NAME
elif [ -d /opt/RStudio/ ]
then
  export SIMG_IMAGE=/opt/RStudio/$SIMG_FILE_NAME
else
  echo 'cant fine a simg file that matches'
  exit
fi


module load singularity;

RSTUDIO_PASSWORD=$PASSWORD singularity exec -B /tgen_labs,/opt,$RSTMP:/tmp,$RSLIB:/var/lib/rstudio-server/,$RSRUN:/var/run/rstudio-server/,${RPRIHOME}/.Renviron:$HOME/.Renviron,${RPRIHOME}/.rstudio:$HOME/.rstudio,${RPRIHOME}/.config:$HOME/.config $SIMG_IMAGE R CMD BATCH 13384_GBM_CellChat_compare.R


#RSTUDIO_PASSWORD=$PASSWORD singularity exec -B /labs/banovich -B /scratch/${USER} -B /packages/petagene/protect_1.3.21-1/bin -B /opt,$RSTMP:/tmp,$RSLIB:/var/lib/rstudio-server/,$RSRUN:/var/run/rstudio-server/ $SIMG_IMAGE R --vanilla CMD BATCH CART_tumor_scGSVA.R

