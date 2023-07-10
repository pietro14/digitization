import h5py
import matplotlib.pyplot as plt
import numpy as np

def plot_image_from_h5(file_path, dataset_name):
    # Open the HDF5 file
    with h5py.File(file_path, 'r') as file:
        # Access the dataset
        dataset = file[dataset_name]

        # Convert the dataset to a NumPy array
        image = dataset[:]

        # Plot the image with origin at the bottom left
        plt.imshow(image, cmap="hot", alpha=1, origin="lower", vmin=0, vmax=np.max(image))

        plt.xticks(np.arange(0, 2305, 256))
        plt.yticks(np.arange(0, 2305, 256))
        plt.xlabel('X')
        plt.ylabel('Y')


        plt.axis('off')
        plt.savefig('image.png')

# Provide the file path and dataset name
file_path = './histograms_Run00001.h5'

dataset_name = 'pic_run1_ev0'

# Call the function to plot the specified image with origin at the bottom left and save it as 'image.png'
plot_image_from_h5(file_path, dataset_name)

