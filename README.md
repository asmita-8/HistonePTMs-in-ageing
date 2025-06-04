# HistonePTMs-in-ageing
1. NucPosSimulator
   a. subsetbedfile_genewise_onePTMatatime.sh : This shell script subsets the main bed file(PTM data) downloded from NIH Roadmap Epigenomics Project gene wise using Bedtools intersect. This script runs one PTM bed file at a time. Download Bedtools before running this script.
   b. run_nucpossimulator.sh : after subsetting the downloaded bed file gene wise, run NucPosSimulator. In one parent directory there will be one directory per gene, in which the output of NucPosSimulator will be saved. Download NucPosSimulator before running this script. Also, after downloading a paramerter "NucLength" needs to changed to 200 from 147(default) in params.txt file. This script too runs one PTM at a time.
   c. round_off_result_bed_file.py : Round off the output of result.bed file of NucPosSimultor to nearest multiple of 200.
   
