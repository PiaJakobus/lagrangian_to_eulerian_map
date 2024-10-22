import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd

# Load the data
file_path = 'out.dat'
data = pd.read_csv(file_path, delim_whitespace=True, header=None)

# Extract columns
x = data.iloc[:, 2]  # Column 3
y = data.iloc[:, 3]  # Column 4
z = data.iloc[:, 4]  # Column 4 for color

cb = plt.scatter(x,y,c=z,s=1,vmin=-1,vmax=1)
plt.colorbar(cb)
plt.show()
