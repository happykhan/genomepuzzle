# genomepuzzle


Uses reseq for read simulation, 

```
docker pull quay.io/biocontainers/reseq:1.1--py38h7ce28ed_4
```


bowtie2-build my_reference.fa my_reference
bowtie2 -p 32 -X 2000 -x my_reference -1 my_data_1.fq -2 my_data_2.fq | samtools sort -m 10G -@ 4 -T _tmp -o my_mappings.bam -

reseq illuminaPE -j 32 -r my_reference.fa -b my_mappings.bam -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq

reseq illuminaPE -j 32 -r my_reference.fa -b my_mappings.bam --statsOnly
reseq illuminaPE -j 32 -s my_mappings.bam.reseq --stopAfterEstimation
reseq illuminaPE -j 32 -R my_reference.fa -s my_mappings.bam.reseq --ipfIterations 0 -1 my_simulated_data_1.fq -2 my_simulated_data_2.fq


To run a simulation with tiles the tile information needs to stay in the read names after the mapping. This means there must not be a space before it, like it is often the case for read archive data. To replace the space on the fly with an underscore the reseq-prepare-names.py script is provided. In this case run the mapping like this:

bowtie2 -p 32 -X 2000 -x my_reference -1 <(reseq-prepare-names.py my_data_1.fq my_data_2.fq) -2 <(reseq-prepare-names.py my_data_2.fq my_data_1.fq) | samtools sort -m 10G -@ 4 -T _tmp -o my_mappings.bam -






## Kleb dataset

from S. David, V. Cohen, S. Reuter, A.E. Sheppard, T. Giani, J. Parkhill, , , G.M. Rossolini, E.J. Feil, H. Grundmann, D.M. Aanensen, Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread among Klebsiella pneumoniae, Proc. Natl. Acad. Sci. U.S.A.
117 (40) 25043-25054,

https://doi.org/10.1073/pnas.2003407117 (2020). 



## outbreak module 

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1592-1