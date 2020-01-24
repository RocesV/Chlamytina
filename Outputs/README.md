# OUTPUTS FORMAT #
It depends on the arguments defined but we can classify the outputs in 4 types:

- *.txt - Example file1_Treatment2-Treatment1.txt : All the files finished with .txt comes from differential expression tests. Syntax: [file used by position in arguments]_[Treatments comparison from condition vector]
- *_background.bed - Possible Universes you can use to run LOLA enrichment analysis. File background is named without _background. Example: file1.bed
- Uniq_***.bed - bed input for LOLA that comes from whole files intersection
- Diff_uniq_***.bed - bed input for LOLA that comes from differential files intersection
- *.pdf - LOLA enrichment results
- *.svg - nVennR algorithm intersection results



