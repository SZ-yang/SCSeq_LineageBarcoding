# I somehow lost how Larry_41093_2000_norm_log_cleaned.h5ad or scBaseEncoderFeat_Z_bs260_tau0.5.csv were made...

import numpy as np
import pandas as pd

# Load the .npy file
data = np.load('/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/Larry_scvi_full_embeddings.npy')

# Convert the data to a DataFrame (if it's 2D or higher dimensional)
df = pd.DataFrame(data)

# Save as a .csv file
df.to_csv('/home/users/kzlin/kzlinlab/projects/scContrastiveLearn/out/kevin/Writeup5/Larry_scvi_full_embeddings.csv', index=False)