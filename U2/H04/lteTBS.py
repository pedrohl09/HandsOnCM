import numpy as np
# Load the .npz file
data = np.load("TBS_data.npz")
loaded_array = data["TBS"]

print(len(loaded_array))