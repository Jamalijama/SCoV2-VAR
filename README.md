# SCoV2-VAR
## A light-weighted, customizable, and open-source database of 12 million SARS-CoV-2 genomes.

This package provides an implementation of the inference pipeline of 
Extraction information from GB file, Fast sequence alignment, Database construction, 
Database application, Gene fragment retrieval program, and Annotation information statistics.

### Running your database construction

Before establishing the database, it is necessary to download the original files containing 
sequence information from the corresponding websites, which will consume a significant 
amount of storage space. It is recommended to perform this task with sufficient memory 
available.

You will need a machine running python environment with any operating systems, 
that need install the dependencies.

Note: You may optionally wish to
    create a
    [Python Virtual Environment](https://docs.python.org/3/tutorial/venv.html)
    to prevent conflicts with your system's Python environment.

    ```bash
    pip3 install -r requirements.txt
    ```

#### If you want to construct a database from scratch, please follow these steps:

1.  In this project, we provide a program to extract sequence information from GB files, 
if needed, please run `ExtractFromGB_0.py`

    ```bash
    python ExtractFromGB_0.py
    ```
    
2.  Run `DBConstuction_1.py` to construct a new SCoV2-VAR database, 
you need to select the corresponding code based on whether the downloaded data 
has already been aligned with the reference sequence. 
If it has been aligned, run the first segment of code in the main function;
otherwise, run the second segment of code. 
Additionally, you also need to modify the sequence file name at the beginning of the program.

    ```bash
    python DBConstuction_1.py
    ```
### Running your database application

Run `Application_2.py` to utilize the already established SCoV2-VAR database, 
input the corresponding numbers and files according to the prompts.

Note: This program currently does not support undo operations, so if an error is made during input, you will need to re-run the program.

```bash
python Application_2.py
```

### Run your gene fragment retrieval and annotation information statistics
In this project, we provide programs to generate corresponding gene fragments, 
protein fragments, and SNVs (Single Nucleotide Variants) within sequences 
based on the SCoV2-VAR database, if needed, please run `GeneRetrieval_3.py`
    
```bash
python GeneRetrieval_3.py
```

In this project, we also provide a program for compiling statistics 
on the annotation files in the database, if needed, please run `AnnotationStatistics_3.py`

```bash
python AnnotationStatistics_3.py
```