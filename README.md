Steps taken to produce output (assumes you have python 3.9x and augustus installed):

1. download https://genome.jgi.doe.gov/portal/Aspni_NRRL3_1/Aspni_NRRL3_1.download.ftp.html
2. extract the files somewhere. Examples in this guide assume the contents are extracted under ~/Aspni_NRRL3
3. run augustus. this generates a new gff file under ~/Aspni_NRRL3_augustus. you may need to create the directory:
```augustus --species=aspergillus_nidulans ~/Aspni_NRRL3/Aspni_NRRL3_1_AssemblyScaffolds.fasta > ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gff```
4. run the following command:
```find_model_diffs.py ~/Aspni_NRRL3/Aspni_NRRL3_1_GeneCatalog_genes_20140311.gff ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gff```
5. run the following two commands. these create .gbk files of the supplied gff files
    - ```gff.py ~/Aspni_NRRL3/Aspni_NRRL3_1_GeneCatalog_genes_20140311.gff```
    - ```gff.py ~/Aspni_NRRL3_augstus/Aspni_NRRL3_1_AssemblyScaffolds.gff```
6. for the first set of results, run
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out [your@email.com] include_ids.txt 0.9```
this may take a long time (~10 hours). there could be errors in the meantime if you are rate limited by ncbi. it is important to replace the [your@email.com] with your own email
7. check if any of the output xml files under blastp_out lack results
8. run ```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out output``` (replace the gemcapy and macse path with your paths)
9. run ```compare_results.py output inexact_matches.csv ~/Aspni_NRRL3/Aspni_NRRL3_1_GeneCatalog_genes_20140311.gbk ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk base_results```
10. check results in ```base_results/results.csv```

other result sets:

90% coverage, 40% identity:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_4 [your@email.com] include_ids.txt 0.9 0.4```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_4 output_9_4```

90% coverage, 60% identity:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_6 [your@email.com] include_ids.txt 0.9 0.6```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_6 output_9_6```

90% coverage, 80% identity:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_8 [your@email.com] include_ids.txt 0.9 0.8```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_8 output_9_8```

90% coverage, 60% identity, realign aggressive:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_6 [your@email.com] include_ids.txt 0.9 0.6```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_6 output_9_6_ra -ra```

90% coverage, 60% identity, realign short:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_6 [your@email.com] include_ids.txt 0.9 0.6```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_6 output_9_6_rs -rs```

90% coverage, 60% identity, realign aggressive+short:
```python retrieve_model_informants.py ~/Aspni_NRRL3_augustus/Aspni_NRRL3_1_AssemblyScaffolds.gbk gbk_out_9_6 [your@email.com] include_ids.txt 0.9 0.6```
```./run_gemcapy.sh ~/gemcapy/gemcapy.py ~/macse_path/ gbk_out_9_6 output_9_6_ra_rs -ra -rs```
