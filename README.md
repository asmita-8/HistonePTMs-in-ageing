# HistonePTMs-in-ageing
1. NucPosSimulator
   a. subsetbedfile_genewise_onePTMatatime.sh : This shell script subsets the main bed file(PTM data) downloded from NIH Roadmap Epigenomics Project gene wise using Bedtools intersect. This script runs one PTM bed file at a time. Download Bedtools before running this script.
   b. run_nucpossimulator.sh : after subsetting the downloaded bed file gene wise, run NucPosSimulator. In one parent directory there will be one directory per gene, in which the output of NucPosSimulator will be saved. Download NucPosSimulator before running this script. Also, after downloading a paramerter "NucLength" needs to changed to 200 from 147(default) in params.txt file. This script too runs one PTM at a time.
   c. round_off_result_bed_file.py : Round off the output of result.bed file of NucPosSimultor to nearest multiple of 200.
   d. discretise_rounded_values.py : Asign value one to the nucleosome where histone PTM is present and 0 where its absent. Input file is bins_for_all_genes.csv where each genes length has been divided into 200 base pair regions(1 nucleosome). Runs one PTM at a time.
   e. concatenate.py : Finally concatenate all the PTM discretised values into one csv file.
   
2. Bit Error Rate
   a. ber.py : Calculates the BER between two samples and plots histogram for all genes.
   b. ber_LV.py : Selects those genes that have BER<0.05 between 3 year and 34 year old sample of heart left ventricle
   c. ber_RV.py : Selects those genes that have BER<0.05 between 3 year and 34 year old sample of heart right ventricle
   d. ber_3.py : Selects those genes that have BER>0.75 between heart left and right ventricle of 3 year old sample. PTMs are grouped into activating, repressive and neutral.
   e. ber_34.py : Selects those genes that have BER>0.75 between heart left and right ventricle of 34 year old sample. PTMs are grouped into activating, repressive and neutral.

3. CNN
   a. cnn_h3k4me3_kfold_34yr.py : CNN model to predict the PTM pattern in 5 nucleosome positions upstream of TSS(1 being the farthest and 5 being the closest to TSS) based on PTM patterns in the gene body regions across 5 fold cross validation for 34 year old sample of hear left ventricle taking 15,338 genes(input file: morethan50_nopromoter_nogenebody_cond_34yr_h3k4me3_heartLV.csv---selected genes having gene body length(nucleosomes) 50 and more).
   b. run_captum.py : This plot reveals which nucleosomes in the gene body region are most important in predicting the PTM pattern in the promoter region. This script uses the saved model from one k-fold.
   c. promoterplusgenebodysignal.py : Plots the signal of PTM in promotyer and gene body region for the two age groups.
