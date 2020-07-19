import os
import sys
import cv2
import numpy as np
import pickle
from cyvlfeat.sift import *
from sklearn.cluster import KMeans
from sklearn.metrics import confusion_matrix

all_bagoffeatures = []


class Image:
    def __init__(self, data, name=""):
        self.data = data            # Grayscale image data is stored
        self.size = data.shape      # Size of the image
        self.name = name            # Name of the image for test dataset
        self.descriptor = []        # Will be initialized when sift descriptors are calculated

    # Calculate the descriptor of the image object according to the method specified (SIFT/dense SIFT)
    def calculateDescriptor(self, method, step_size=0):
        if method == "sift":
            sift = cv2.xfeatures2d.SIFT_create(1000)
            key_points, descriptor = sift.detectAndCompute(self.data, None)

            if len(key_points) == 0:
                # Sample keypoint manually
                keypoint = cv2.KeyPoint(x = 15, y = 15, _size = 1)
                descriptor = cv2.sift.compute(self.data, [keypoint])

            self.descriptor = descriptor

        elif method == "densesift":
            key_points, descriptor = dsift(self.data, step=step_size)

            self.descriptor = descriptor

        return self


class Classifier:
    def __init__(self):
        self.train_directory = "the2_data/train"
        self.validation_directory = "the2_data/validation"
        self.test_directory = "test"
        self.unique_labels = []         # Store unique class labels
        self.train_images = []          # Store grayscale images for train set
        self.train_labels = []          # Store true labels for each image in train set
        self.validation_images = []     # Store grayscale images for validation set
        self.validation_labels = []     # Store true labels for each image in validation set
        self.predicted_labels = []      # Store predicted labels for images in validation set
        self.test_images = []           # Store grayscale images for test dataset
        self.kmeans = None              # Will be initialized

    def readImages(self):
        # Read images in the train dataset
        for filename in os.listdir(self.train_directory):
            self.unique_labels.append(filename)
            class_directory = os.path.join(self.train_directory, filename)

            for image_name in os.listdir(class_directory):
                data = cv2.imread(self.train_directory + "/" + filename + "/" + image_name)
                # Convert image to grayscale
                gray_data = cv2.cvtColor(data, cv2.COLOR_BGR2GRAY)
                self.train_images.append(Image(gray_data))
                self.train_labels.append(filename)

        # Read images in the validation dataset
        for filename in os.listdir(self.validation_directory):
            class_directory = os.path.join(self.validation_directory, filename)

            for image_name in os.listdir(class_directory):
                data = cv2.imread(self.validation_directory + "/" + filename + "/" + image_name)
                # Convert image to grayscale
                gray_data = cv2.cvtColor(data, cv2.COLOR_BGR2GRAY)
                self.validation_images.append(Image(gray_data))
                self.validation_labels.append(filename)

    def readTestImages(self):
        # Reset validation_images to store test images for this time
        self.validation_images = []
        # Read images in the test dataset
        for image_name in os.listdir(self.test_directory):
            data = cv2.imread(self.test_directory + "/" + image_name)
            # Convert image to grayscale
            gray_data = cv2.cvtColor(data, cv2.COLOR_BGR2GRAY)
            self.validation_images.append(Image(gray_data, image_name))

    # Detect keypoints of the images by applying SIFT or dense SIFT algorithms
    def detectKeypoints(self, method, step_size):
        # sift/dsift returns descriptor arrays with shape (#key_points, 128)
        self.train_descriptors = np.ndarray((0, 128))

        for i, image in enumerate(self.train_images):
            image.calculateDescriptor(method, step_size)
            #print("Image", i, "is done.")
            self.train_descriptors = np.concatenate((self.train_descriptors, image.descriptor))

    # Calculate the Bag of Features representation for each image in the training set
    def calculateBOF(self, kmeans_k, pickle_name):
        global all_bagoffeatures

        # Perform k-means clustering
        kmeans = KMeans(n_clusters=kmeans_k, random_state=0).fit(self.train_descriptors)

        #pickle.dump(kmeans, open(pickle_name, "wb"))

        # Calculate BOF for each image
        for image in self.train_images:
            # Create histogram with bin_count = kmeans_k
            histogram = np.zeros(kmeans_k, float)
            labels = kmeans.predict(image.descriptor)

            np.add.at(histogram, labels, 1)

            # Normalize histogram
            histogram = histogram / histogram.sum(axis=0)
            image.bag_of_features = histogram

            all_bagoffeatures.append(histogram)

    def predict(self, method, step_size, kmeans_k, knn_k, pickle_name):
        #self.kmeans = pickle.load(open(pickle_name, "rb"))

        for i, image in enumerate(self.validation_images):
            image.calculateDescriptor(method, step_size)

            histogram = np.zeros(kmeans_k, float)
            labels = self.kmeans.predict(image.descriptor)

            np.add.at(histogram, labels, 1)

            # Normalize histogram
            histogram = histogram / histogram.sum(axis=0)
            image.bag_of_features = histogram

            result = self.knn(knn_k, image)
            self.predicted_labels.append(result)


    # Return predicted label for test_image
    def knn(self, k, test_image):
        global all_bagoffeatures

        # Get distance between the features of the test_image and train images
        distances = np.linalg.norm(test_image.bag_of_features - all_bagoffeatures, axis=1)
        distances_labels = list(zip(list(distances), self.train_labels))

        # Sort with respect to distances
        sorted_distances = sorted(distances_labels, key = lambda x: x[0])

        # Find labels for k nearest neighbors
        nn_labels = []

        for i in range(k):
            nn_labels.append(sorted_distances[i][1])

        # Find the majority label
        result_label = max(set(nn_labels), key=nn_labels.count)
        return result_label

    def accuracy(self):
        correct = 0

        for i, label in enumerate(self.predicted_labels):
            if label == self.validation_labels[i]:
                correct += 1

        accuracy = correct / len(self.validation_labels) * 100
        return accuracy


    def printConfusionMatrix(self):
        print(confusion_matrix(self.validation_labels, self.predicted_labels, labels=self.unique_labels))


    def writeResults(self):
        with open("results.txt", "w+") as output:
            for i, image in enumerate(self.validation_images):
                output.write(image.name + ": " + self.predicted_labels[i] + "\n")



if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Please run as python3 the2.py sift/densesift step_size kmeans_k knn_k")
        sys.exit()

    sift_type = sys.argv[1]
    step_size = int(sys.argv[2])
    kmeans_k = int(sys.argv[3])
    knn_k = int(sys.argv[4])

    classifier = Classifier()
    classifier.readImages()
    #classifier.readTestImages()
    classifier.detectKeypoints(sift_type, step_size)
    classifier.calculateBOF(kmeans_k, "kmeans128_densesift3.pkl")
    classifier.predict(sift_type, step_size, kmeans_k, knn_k, "kmeans128_densesift3.pkl")
    #classifier.writeResults()
    #print("Accuracy is:", classifier.accuracy())
    #classifier.printConfusionMatrix()
