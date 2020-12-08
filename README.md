# COMIRI

__*Please note: this tool is not finished. This README will be updated as the tool is developped*__

## Introduction:

**COMIRI** (**CO**sensus **M**etacaller for **I**mmune **R**eceptor **I**dentification) automates the extraction of T-cell and B-cell CDR3 sequences from bulk RNA sequencing. The program can run the following algorithms:  
* [CATT](https://github.com/GuoBioinfoLab/CATT)  
* [ImRep](https://github.com/mandricigor/imrep)  
* [MiXCR](https://mixcr.readthedocs.io/en/master/index.html)  
* [TRUST3](https://bitbucket.org/liulab/ng-bcr-validate/src/master/TRUST3/)  
* [TRUST4](https://github.com/liulab-dfci/TRUST4)  
* [V'DJer](https://github.com/mozack/vdjer)  

## Overview:


![Flowchart](images/COMIRI_flowchart.jpg?raw=true)

### Features:

<a href=''><img src='images/Parsl_logo.jpg' align="right" width="150" /></a> <a href=''><img src='images/Docker_logo.jpg' align="right" width="150" /></a> <a href=''><img src='images/Singularity_logo.jpg' align="right" width="150" /></a>


* Single command to run all tools and dependencies, with adjustable parameters for each tool, filtering step and consensus method  
* Possibility of running with Docker or Singularity containers, for increased portability and reproducibility  
* Parallelization with [Parsl](https://parsl.readthedocs.io/en/stable/index.html), a flexible and scalable parallel programming library for Python. Parsl allows for a modular allocation of resources to run each step concurrently, while respecting any dependencies between programs    


## Installation:



## Usage:





## 

[CSHL Biological Data Science Conference (2020)](https://meetings.cshl.edu/posters/data20/images/viewer.html?file=data_20_142.pdf)

