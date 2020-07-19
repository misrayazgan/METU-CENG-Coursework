import cv2
from scipy import signal
import numpy as np
import operator
import itertools
from matplotlib import pyplot as plt

epsilon = 0.00000001
dataset_directory = "dataset"

name_dict = {}              # Dictionary keeping (key:image name, value: index in RGB_images array)
image_names = []
RGB_images = []
grayscale_histograms = []
color_histograms = []
gradient_histograms = []
query_images = []

class Image:
    def __init__(self, data, name):
        self.data = data
        self.name = name
        self.size = data.shape
        self.matches = {}

    def resize_image(self):
        print("Image is resized, before: ", self.size)

        if self.size[0] == 640:
            self.data = cv2.resize(self.data, (480, 640))       #(width, height) --> (#columns, #rows)
        elif self.size[1] == 640:
            self.data = cv2.resize(self.data, (640, 480))

        self.size = self.data.shape
        print(self.size)

    def grayscale_histogram(self, data, bin_count):
        # Convert the image to grayscale

        grayscale_data = cv2.cvtColor(data, cv2.COLOR_BGR2GRAY)
        flat_grayscale = grayscale_data.ravel()

        step_size = 256 // bin_count
        hist = np.zeros(bin_count, float)

        bins = np.arange(0, 256 + step_size, step_size)         # Array keeping bin values
        digitized = np.digitize(flat_grayscale, bins) - 1           # Array keeping indices

        bin_indices, counts = np.unique(digitized, return_counts=True)

        # Assign counts to bins
        for i, bin_index in enumerate(bin_indices):
            hist[bins[bin_index] // step_size] = counts[i]

        hist = hist / hist.sum(axis=0)

        return hist

    def color_histogram(self, data, bin_count):
        step_size = 256 // bin_count

        # Each channel has bin_count many bins
        # Combine them --> bin_count^3
        hist = np.zeros(bin_count ** 3, float)

        # Prepare index triplets to combine bins from R, G and B
        combined_indices = {key: i for i, key in enumerate(itertools.product(np.arange(bin_count), repeat=3))}
        # returns {(0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2, (0, 1, 1): 3, (1, 0, 0): 4, (1, 0, 1): 5, (1, 1, 0): 6, (1, 1, 1): 7} for 2 bins

        r_column = (data[:,:,0]).ravel()
        g_column = (data[:,:,1]).ravel()
        b_column = (data[:,:,2]).ravel()

        r_id = r_column // step_size
        g_id = g_column // step_size
        b_id = b_column // step_size

        rgb_id = list(zip(r_id, g_id, b_id))

        rgb_bin_indices = [combined_indices[triplet] for triplet in rgb_id]

        np.add.at(hist, rgb_bin_indices, 1)
        hist = hist / hist.sum(axis=0)
        return hist

    def hog(self, data, bin_count):
        # Convert the image to grayscale
        grayscale_data = cv2.cvtColor(data, cv2.COLOR_BGR2GRAY)

        v_filter = np.array([[1, 1, 1], [0, 0, 0], [-1, -1, -1]])
        h_filter = np.array([[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]])

        horizontal_filtered = signal.convolve2d(grayscale_data, h_filter, mode="same")
        vertical_filtered = signal.convolve2d(grayscale_data, v_filter, mode="same")

        orientation = np.degrees(np.arctan(vertical_filtered / (horizontal_filtered + epsilon))).astype(int).ravel()       # Keep orientation in degrees
        orientation = np.where(orientation >= 0, orientation, orientation + 180)                     # If orientation is negative, add 180 degrees
        magnitude = np.sqrt(vertical_filtered ** 2 + horizontal_filtered ** 2)

        magnitude = magnitude.ravel()

        # Histogram range is 0 - 180 degrees
        hist = np.zeros(bin_count, float)
        step_size = 180 // bin_count
        bins = np.arange(0, 180, step_size)         # Array keeping bin values

        indices = (orientation // step_size).astype(int)
        np.add.at(hist, indices, magnitude)

        hist = hist / hist.sum(axis = 0)

        return hist

    def split_image(self, level):
        # Keep grids as 3d numpy arrays in a list
        grids = []
        row_col_count = 2 ** (level - 1)

        row_step = self.size[0] // row_col_count
        col_step = self.size[1] // row_col_count

        h_grids = np.hsplit(self.data, row_col_count)
        for h in h_grids:
            vs = np.vsplit(h, row_col_count)
            for v in vs:
                grids.append(v)

        return grids



    def grid_based(self, type, bin_count, level):
        # If size of the image is different, resize it
        if self.size != (640, 480, 3) and self.size != (480, 640, 3):
            self.resize_image()

        # If level == 1, 1 grid
        # If level == 2, 4 grids
        # If level == 3, 16 grids
        grids = []
        if level != 1:
            grids = self.split_image(level)
        else:
            grids = [self.data]

        hist_concat = np.array([], float)
        if type == "grayscale":
            for grid in grids:
                hist = self.grayscale_histogram(grid, bin_count)
                hist_concat = np.concatenate((hist_concat, hist), axis=0)

        elif type == "color":
            for grid in grids:
                hist = self.color_histogram(grid, bin_count)
                hist_concat = np.concatenate((hist_concat, hist), axis=0)

        elif type == "gradient":
            for grid in grids:
                hist = self.hog(grid, bin_count)
                hist_concat = np.concatenate((hist_concat, hist), axis=0)

        return hist_concat


    def find_matches(self, type):
        # type: can be one of "grayscale", "color", "gradient"

        global grayscale_histograms, color_histograms, gradient_histograms, name_dict, image_names
        self_index = name_dict[self.name]

        if type == "grayscale":
            # Calculate euclidean distance between self.data's histogram and histograms of the other images in the dataset
            distances = np.linalg.norm(grayscale_histograms[self_index] - grayscale_histograms, axis=1)
            self.matches = dict(zip(image_names, distances))
        elif type == "color":
            distances = np.linalg.norm(color_histograms[self_index] - color_histograms, axis=1)
            self.matches = dict(zip(image_names, distances))
        elif type == "gradient":
            distances = np.linalg.norm(gradient_histograms[self_index] - gradient_histograms, axis=1)
            self.matches = dict(zip(image_names, distances))


    def write_output(self, filename):
        # Sort the match images in ascending order
        self.matches = dict(sorted(self.matches.items(), key=operator.itemgetter(1)))

        with open(filename, "a") as output:
            output.write(self.name + ": ")

            for key in self.matches:
                # Write euclidean distance and name of the match image
                output.write(str(self.matches[key]) + " " + key + " ")

            output.write("\n")



def read_images(filename):
    global RGB_images, name_dict, image_names

    # Read image names from the input file
    with open(filename, "r") as input:
        file = input.read()
        image_names = file.split("\n")

        # Remove the last blank lines
        while len(image_names[-1]) == 0:
            image_names = image_names[:-1]

        for i, name in enumerate(image_names):
            name_dict[name] = i
            data = cv2.imread(dataset_directory + "/" + name)
            RGB_image = Image(data, name)
            RGB_images.append(RGB_image)



def read_query_images(filename):
    global query_images, name_dict, RGB_images

    with open(filename, "r") as input:
        file = input.read()
        query_names = file.split("\n")

        # Remove the last blank lines
        while len(query_names[-1]) == 0:
            query_names = query_names[:-1]

        for name in query_names:
            query_image = RGB_images[name_dict[name]]
            query_images.append(query_image)


def calculate_grayscale_histograms():
    global grayscale_histograms

    for i, image in enumerate(RGB_images):
        #print("Image", i)
        grayscale_histograms.append(image.grid_based("grayscale", 16, 1))



def calculate_color_histograms():
    global color_histograms

    for i, image in enumerate(RGB_images):
        #print("Image", i)
        color_histograms.append(image.grid_based("color", 16, 1))


def calculate_hog():
    global gradient_histograms

    for i, image in enumerate(RGB_images):
        #print("Image", i)
        gradient_histograms.append(image.grid_based("gradient", 45, 1))


def get_query_results():
    output = open("output_color_l1.txt","w+")

    for image in query_images:
        image.find_matches("color")
        image.write_output("output_color_l1.txt")


if __name__ == "__main__":
    read_images("src/images.dat")
    read_query_images("src/test_queries.dat")
    #calculate_grayscale_histograms()
    #calculate_color_histograms()
    #calculate_hog()
    color_histograms = np.load("colorhist/colorhist_16_l1.npy")
    get_query_results()
