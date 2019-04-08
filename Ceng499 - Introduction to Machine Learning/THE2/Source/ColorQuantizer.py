#!/usr/bin/env python
# -*- coding: utf-8 -*-
from MyKMeans import MyKMeans
import numpy as np
from sklearn.datasets import load_sample_image
import os
from PIL import Image
from scipy import misc
from sklearn.utils import shuffle
from sklearn.metrics import pairwise_distances_argmin


class ColorQuantizer:
    """Quantizer for color reduction in images. Use MyKMeans class that you implemented.
    
    Parameters
    ----------
    n_colors : int, optional, default: 64
        The number of colors that wanted to exist at the end.
    random_state : int, RandomState instance or None, optional, default: None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Read more from:
    http://scikit-learn.org/stable/auto_examples/cluster/plot_color_quantization.html
    """
    
    def __init__(self, n_colors=64, random_state=None):
        self.n_colors = n_colors
        self.random_state = random_state

    
    def read_image(self, path):
        """Reads jpeg image from given path as numpy array. Stores it inside the
        class in the image variable.
        
        Parameters
        ----------
        path : string, path of the jpeg file
        """
                        
        self.image = misc.imread(path)              # read image and save into numpy array
        
        width, height, depth = self.image_shape = tuple(self.image.shape)           # save shape of the array

        self.image = np.array(self.image, dtype = np.float64) / 255

        self.image = np.reshape(self.image, (width * height, depth))


        
    def recreate_image(self, path_to_save):
        """Reacreates image from the trained MyKMeans model and saves it to the
        given path.
        
        Parameters
        ----------
        path_to_save : string, path of the png image to save
        """
        image_samples = shuffle(self.image, random_state = self.random_state)[:2000]

        # train MyKMeans model
        kmeans = MyKMeans(init_method = "random", n_clusters = self.n_colors, random_state = self.random_state)
        kmeans.initialize(image_samples)
        kmeans.fit(image_samples)

        #print kmeans.labels

        self.cluster_centers = kmeans.cluster_centers
        #print kmeans.cluster_centers

        predicted_labels = kmeans.predict(self.image)

        #print predicted_labels

        d = kmeans.cluster_centers.shape[1]

        recreated_image = np.zeros(self.image_shape)

        label_index = 0

        for i in range(self.image_shape[0]):
            for j in range(self.image_shape[1]):
                recreated_image[i][j] = kmeans.cluster_centers[predicted_labels[label_index]]
                label_index += 1

        #plt.imshow(recreated_image)
        #plt.show()

        misc.imsave(path_to_save, recreated_image)


    def export_cluster_centers(self, path):
        """Exports cluster centers of the MyKMeans to given path.

        Parameters
        ----------
        path : string, path of the txt file
        """
        np.savetxt(path, self.cluster_centers, delimiter = "     ", fmt = "%.18g", newline = os.linesep)
        
        
    def quantize_image(self, path, weigths_path, path_to_save):
        """Quantizes the given image to the number of colors given in the constructor.
        
        Parameters
        ----------
        path : string, path of the jpeg file
        weigths_path : string, path of txt file to export weights
        path_to_save : string, path of the output image file
        """
        self.read_image(path)
        self.recreate_image(path_to_save)
        self.export_cluster_centers(weigths_path)
        

if __name__ == "__main__":

    quantizer = ColorQuantizer(n_colors = 64)

    quantizer.quantize_image("metu.jpg", "metu_cluster_centers.txt", "my_quantized_metu.png")

    quantizer.quantize_image("ankara.jpg", "ankara_cluster_centers.txt", "my_quantized_ankara.png")

