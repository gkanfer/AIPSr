# <p>  <b>AI-photoswitchable screening (AI-PS) </b> </p>
<img src="https://github.com/gkanfer/AI-PS/raw/master/logoAIPS.png" width="250" title="cellpose" alt="cellpose" align="right" vspace = "10">

This project permits the development and utilization of a pooled optical screen by making use of photoactivation and sorting. AI-PS algorithms provide a structure for segmenting cells, training machine learning algorithms to categorize cellular phenotypes, and applying these for executing a screen on the fly. To learn more, take a look at Kanfer et al., 2021.

### About photoswitchable Crispr screening 
  
AI-PS is a platform that uses machine learning and deep learning algorithms to facilitate pooled genetic screening for subcellular image phenotypes. This method reduces time, cost and complexity compared to standard screening methods. AI-PS can detect and label cells according to subcellular protein localization, abundance, size and shape. It is compatible with adherent tissue culture cells and is accessible with no need for specialized flow instrumentation. Results of this platform have been validated by identifying PINK1 as the only known reported hit required for Parkin translocation to damaged mitochondria. It has also been used to explore a completely different protein translocation process undetectable via FACS.


### CITATION
Kanfer, G., Sarraf,S., Mamman, K., Baldwin, H., Johnson, K., Kampmann, M., Ward, M., Lippincott-Schwartz, J., Youle, R. Image-based pooled whole genome CRISPR screening for intracellular phenotypes – Parkin and TFEB subcellular localization. J Cell Biol. 2021 Feb 1;220(2):e202006180. doi: 10.1083/jcb.202006180

### Installation
```
library(devtools)
install_github("gkanfer/AIPSr")
```
#### System requirements
This software is supported for running on Linux, Windows and Mac OS. It requires a Mac OS later than Yosemite for the graphical interface and 8GB of RAM to run, with 16GB-32GB being recommended for larger images and 3D volumes. If there are any installation issues, please open an issue.

#### Dependencies
```
library(EBImage)
library(e1071)
library(data.table)
library(plyr)
library(dplyr)
library(gtools)
library(yaImpute)
library(outliers)
library(ggplot2)
library(ggbiplot)
library(progress)
library(parallel)
library(reticulate)
use_condaenv("tf-keras")
library(keras)
library(tensorflow)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(caret)
library(ROCR)
```

### Model training examples from AIPS paper
 - **support vector machine (SVM) classification:** A support vector machine (SVM) classification model is trained on images of cells with either cytosolic or mitochondrial GFP-Parkin, utilizing 18 features computed from 2,500 single-cell images . The features were calculated by EBImage, a R image processing and analysis package (Pau et al., 2010), and it was determined that this binary switch in the Parkin location was suitable for detection by this model.
**Deploy:**
[PARKIN screen deploy code](https://github.com/gkanfer/AI-PS/tree/master/Parkin_screen)



- **Convolutional neural network (CNN) classification:** A CNN model is created using TensorFlow and Keras to classify TFEB localization. 107,226 example images of GFP-TFEB used, 80% for training, 15% for validation and 5% for testing. The image input size is 150x150 pixels, with three steps of convolution and max pooling at a learning rate of 1e-4. 50 epochs and a batch size of 200 were used for training, with the model weights being saved after each epoch to prevent overfitting. Brightness augmentation (10%-90%) was also applied to the training data set to account for variations in fluorescence signal intensity.
**Deploy:**
[TFEB screen deploy code](https://github.com/gkanfer/AI-PS/tree/master/TFEB_screen)

### Targted sgRNA  mapping and read count
```
#!/bin/sh
module load bcl2fastq
module load umitools
module load trimmomatic
module load fastxtoolkit
module load bowtie/1.1.2
module load samtools
module load umitools
input_path=$(pwd)
mkdir tmp
bcl2fastq --no-lane-splitting --runfolder-dir ./190917_NB552201_0006_AHMMNVAFXY/ --output-dir tmp &> tmp/bcl2fastq.log
bcl2fastq --no-lane-splitting --runfolder-dir $input_path --output-dir tmp &> tmp/bcl2fastq.log
cd tmp
gunzip Undetermined_S0_R1_001.fastq.gz
umi_tools extract --stdin=Undetermined_S0_R1_001.fastq --extract-method='regex' --bc-pattern="(?P<umi_1>.{3}T{1}.{3}T{1}.{3}T{1}.{8})(?P<discard_1>GCACAAAAGGAAAC.*AGTATCCCTTGGAGAACCACCTTGTTGG){s<=2}(.*)(?P<discard_2>GTTTAAGAGCTAAGCT.*){s<=2}" --log=processed.log --stdout processed.fastq.gz
java -jar $TRIMMOJAR SE -phred33 processed.fastq.gz trimmed_processed.fastq.gz SLIDINGWINDOW:4:15 MINLEN:19
gunzip trimmed_processed.fastq.gz
#fastq-sample -n 69000000 trimmed_processed.fastq -o trimmed_processed_sample
fastx_trimmer -f 1 -l 19 -i trimmed_processed.fastq -o trimmed_processed_fastx.fastq
#gzip trimmed_processed_fastx.fastq
path="~Crispri_h"

#gunzip trimmed_processed_fastx.fastq.gz

bowtie --threads 10 --tryhard -n 0 -k 1 -l 19 $path -q trimmed_processed_fastx.fastq --chunkmbs 200 --sam | samtools view -Sb - > Post_Map_e.bam

samtools sort Post_Map_e.bam -o Sorted_Post_Map_e.bam

samtools index Sorted_Post_Map_e.bam

umi_tools count --per-contig -I Sorted_Post_Map_e.bam > e_UMI_countList_test_e.txt
```

### Demo
 <img src="https://github.com/gkanfer/AI-PS/raw/master/Video_tfeb_2.gif" width="250" align="right" hspace = 30 vspace = "0">
<img src="https://github.com/gkanfer/AI-PS/raw/master/Parkin_deploy.gif" width="350" title="cellpose" alt="cellpose" align="left" vspace = "0">
<br/><br/>
<br/><br/>
<br/><br/>
<br/><br/>
<br/><br/>

<br>
<br>
<br>
<br>
<br>

