
# Combine library and scores and output training and test data
Rscript split_data.r > split_data.out

# Do lasso analysis
Rscript regression_analysis.r > regression_analysis.out

# Train PAP model
conda activate tensorflow
python train_cnn2w.py

