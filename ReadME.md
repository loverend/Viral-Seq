# Viral-Seq Single Cell / Bulk RNA seq Analysis

- Refactored by Lauren Overend

- Adapted from Viral-Track https://www.nature.com/articles/s41597-019-0116-4, Author Pierre Bost)

-------------------------------------------------


**Rescomp Dependencies @LaurenOverend:** 

```
module load UMI-tools/1.0.1-foss-2020a-Python-3.8.2 
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load SAMtools/1.10-GCC-9.3.0
module load STAR/2.7.3a-GCC-9.3.0
module load Subread/2.0.1-GCC-9.3.0
```
## R dependencies 
```




```


The first step consists in creating a **STAR** index that include both host and virus reference genomes.
38
To do so first download the **ViruSite**  [genome reference database](http://www.virusite.org/index.php?nav=download). Host genome has also to be downloaded from the [ensembl website](https://www.ensembl.org/info/data/ftp/index.html). 
39
This can take some time and requires large amount of memory and storage space : please check that you have at least 32 GB of RAM and more than 100GB of avaible memory.
40
​
41
The  **STAR** index can now be built by typing :
42
​
No file chosen
Attach files by dragging & dropping, selecting or pasting them.
Styling with Markdown is supported
@loverend
Commit changes
Commit summary
Create README.md
Optional extended description
Add an optional extended description…
 Commit directly to the main branch.
 Create a new branch for this commit and start a pull request. Learn more about pull requests.
 
Footer
© 2023 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
Editing Viral-Seq/README.md at main · loverend/Viral-Seq
