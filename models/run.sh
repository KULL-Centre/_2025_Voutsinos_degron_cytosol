
# conda activate tensorflow
# python PAP.py -m  cnn2w1 -s 0 --skip-nonnatural all_data.fasta > all_data_pap.txt
# python PAP.py -m human30 -s 0 --skip-nonnatural all_data.fasta > all_data_h30.txt
# python PAP.py -m human25 -s 0 --skip-nonnatural all_data.fasta > all_data_h25.txt

# Combine library and scores and output training and test data
Rscript model_data.r > model_data.out

# Do lasso analysis
Rscript regression_analysis.r > regression_analysis.out

# Train PAP model
conda activate tensorflow
python train_cnn2w.py

