#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import distance

EPSILON = 0.0001


class MyKMeans:
	"""K-Means clustering similar to sklearn 
	library but different.
	https://goo.gl/bnuM33

	But still same.

	Parameters
	----------
	n_clusters : int, optional, default: 3
		The number of clusters to form as well as the number of
		centroids to generate.
	init_method : string, optional, default: 'random'
		Initialization method. Values can be 'random', 'kmeans++'
		or 'manual'. If 'manual' then cluster_centers need to be set.
	max_iter : int, default: 300
		Maximum number of iterations of the k-means algorithm for a
		single run.	
	random_state : int, RandomState instance or None, optional, default: None
		If int, random_state is the seed used by the random number generator;
		If RandomState instance, random_state is the random number generator;
		If None, the random number generator is the RandomState instance used
		by `np.random`.
	cluster_centers : np.array, used only if init_method is 'manual'.
		If init_method is 'manual' without fitting, these values can be used
		for prediction.
	"""

	def __init__(self, init_method="random", n_clusters=3, max_iter=300, random_state=None, cluster_centers=[]):
		self.init_method = init_method
		self.n_clusters = n_clusters
		self.max_iter = max_iter
		self.random_state = np.random.RandomState(random_state)
		if init_method == "manual":
			self.cluster_centers = cluster_centers
		else:
			self.cluster_centers = []

	def fit(self, X):
		"""Compute k-means clustering.
		Parameters
		----------
		X : array-like, shape=(n_samples, n_features)
			Training instances to cluster.
		Returns
		----------
		self : MyKMeans
		"""
		previous_centers = np.zeros(self.cluster_centers.shape)
		difference = np.linalg.norm(self.cluster_centers - previous_centers, None)
		self.labels = []				                	# store the clusters that each sample belongs to
															# labels for the clusters are indexes in self.cluster_centers e.g. 0,1,2,..
		number_of_iterations = 1                                                                
		
		while difference > EPSILON and number_of_iterations <= self.max_iter:
			# for each training sample find the closest clusters and update self.labels
			for i in range(X.shape[0]):
				distances_to_centroids = np.linalg.norm(X[i] - self.cluster_centers, axis = 1)
				min_distance_cluster = np.argmin(distances_to_centroids)
				self.labels.append(min_distance_cluster)

			previous_centers = self.cluster_centers
			print "miso"
			# find the new cluster centers
			for i in range(self.n_clusters):
				points = []
				for j in range(len(self.labels)):
					if self.labels[j] == i:
						points.append(X[j])

				self.cluster_centers[i] = np.mean(points, axis = 0)
			
			
			difference = np.linalg.norm(previous_centers - self.cluster_centers, None)
			
			number_of_iterations += 1

		
		return self

	def initialize(self, X):
		""" Initialize centroids according to self.init_method
		Parameters
		----------
		X : array-like, shape=(n_samples, n_features)
			Training instances to cluster.
		Returns
		----------
		self.cluster_centers : array-like, shape=(n_clusters, n_features)
		"""
		shape = (0, X.shape[1])

		if self.init_method == "random":
			np_index = self.random_state.permutation(X.shape[0])[:self.n_clusters]

			self.cluster_centers = X[np_index]
			
		elif self.init_method == "kmeans++":
			self.cluster_centers = np.empty(shape, float)
             
			first_index = self.random_state.randint(X.shape[0])									# randomly select the first index
			first_centroid = X[first_index]                                                     # randomly select the first centroid			

			self.cluster_centers = np.append(self.cluster_centers, [first_centroid], axis = 0)    # append first centroid to the numpy array
			
			# find the centroids after randomly selecting the first one
			for i in range(self.n_clusters - 1):                
				max_distance = 0
				new_centroid = np.empty(shape, float)

				# find the new centroid with the highest total euclidean distance to the previous centroids
				for np_x in X:                  
					total_distance_to_previous = 0

					for centroid in self.cluster_centers:
						total_distance_to_previous += distance.euclidean(np_x, centroid)

					if total_distance_to_previous > max_distance:
						max_distance = total_distance_to_previous
						new_centroid = np_x

				self.cluster_centers = np.append(self.cluster_centers, [new_centroid], axis = 0)

		return self.cluster_centers


	def predict(self, X):
		"""Predict the closest cluster each sample in X belongs to.
		In the vector quantization literature, `cluster_centers` is called
		the code book and each value returned by `predict` is the index of
		the closest code in the code book.
		Parameters
		----------
		X : array-like, shape = [n_samples, n_features]
			New data to predict.
		Returns
		-------
		labels : array, shape [n_samples,]
			Index of the cluster each sample belongs to.
		"""
		labels = []                  # store the clusters that each sample belongs to
		
		for i in range(len(X)):
			distances_to_centroids = np.linalg.norm(X[i] - self.cluster_centers, axis = 1)
			labels.append(np.argmin(distances_to_centroids))

		return labels

	def fit_predict(self, X):
		"""Compute cluster centers and predict cluster index for each sample.
		Convenience method; equivalent to calling fit(X) followed by
		predict(X).
		Parameters
		----------
		X : {array-like, sparse matrix}, shape = [n_samples, n_features]
			New data to transform.
		Returns
		-------
		labels : array, shape [n_samples,]
			Index of the cluster each sample belongs to.
		"""
		self.fit(X)
		return self.predict(X)


if __name__ == "__main__":
	if __name__ == "__main__":
		X = np.array([[1, 2], [1, 4], [1, 0],
					  [4, 2], [4, 4], [4, 0]])
		kmeans = MyKMeans(n_clusters=2, random_state=0, init_method='kmeans++')
		print kmeans.initialize(X)
		#test = KMeans(n_clusters = 2, init = "k-means++", n_init = 1, random_state = 0, algorithm = "full")
		#test.fit(X)
		# [[4. 4.]
		#  [1. 0.]]
		kmeans = MyKMeans(n_clusters=2, random_state=0, init_method = 'random')
		print kmeans.initialize(X)
		# [[4. 0.]
		#  [1. 0.]]
		kmeans.fit(X)
		print kmeans.labels
		# array([1, 1, 1, 0, 0, 0])
		print kmeans.predict([[0, 0], [4, 4]])
		# array([1, 0])
		print kmeans.cluster_centers
		# array([[4, 2],
		#       [1, 2]])