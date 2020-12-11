# COMIRI

__*Please note: this tool is not finished. This README will be updated as the tool is developped.*__

## Introduction:

**COMIRI** (**CO**sensus **M**etacaller for **I**mmune **R**eceptor **I**dentification) automates the extraction of T-cell and B-cell CDR3 sequences from bulk RNA sequencing. The program can run the following algorithms:  
* [CATT](https://github.com/GuoBioinfoLab/CATT)  
* [MiXCR](https://mixcr.readthedocs.io/en/master/index.html)  
* [TRUST3](https://bitbucket.org/liulab/ng-bcr-validate/src/master/TRUST3/)  
* [TRUST4](https://github.com/liulab-dfci/TRUST4)  
* [V'DJer](https://github.com/mozack/vdjer)  
<!-- * [ImRep](https://github.com/mandricigor/imrep)   -->

## Overview:


![](images/COMIRI_flowchart.jpg?raw=true)

### Features:

<img align="right" width="200" src='images/logos.jpg'>

* Single command to run all tools and dependencies, with adjustable parameters for each tool, filtering step and consensus method  
* Possibility of running with Docker or Singularity containers, for increased portability and reproducibility  
* Parallelization with [Parsl](https://parsl.readthedocs.io/en/stable/index.html), a flexible and scalable parallel programming library for Python. Parsl allows for a modular allocation of resources to run each step concurrently, while respecting any dependencies between programs    


## Installation:

```
pip install comiri
```

## Usage:

```
comiri -f FASTQ_PATH -o OUTPUT_DIR <OPTIONS>
```

`-f, --fastq_path FASTQ_PATH`    Path of the input fastq files  
`-o, --output_dir OUTPUT_DIR`    Path of the output directory  

Options:  
`-h, --help`    Show this help message and exit  
`-g, --genome_version {hg19,hg38}`    Genome version to use  
`-p, --run_program {TRUST3,TRUST4,MiXCR,VDJer,CATT}`    List of programs to run  
`-m, --run_mode {docker,singularity,local}`    Run mode for pipeline  
`-r, --receptor {TCR,BCR-heavy,BCR-light}`    List of receptors to extract  
`-c, --config CONFIG`     Config file to use for running Parsl  
`-j, --json JSON`    JSON file to use for running Parsl  
`-s, --single_end`    Run in single-end mode  



## Poster:

[CSHL Biological Data Science Conference (2020)](https://meetings.cshl.edu/posters/data20/images/viewer.html?file=data_20_142.pdf)

