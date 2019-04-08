#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import distance
import copy

class MyKMedoids:
	"""KMedoids implementation parametric with 'pam' and 'clara' methods.

	Parameters
	----------
	n_clusters : int, optional, default: 3
		The number of clusters to form as well as the number of medoids to
		determine.
	max_iter : int, default: 300
		Maximum number of iterations of the k-medoids algorithm for a
		single run.
	method : string, default: 'pam'
		If it is pam, it applies pam algorithm to whole dataset 'pam'.
		If it is 'clara' it selects number of samples with sample_ratio and applies
			pam algorithm to the samples. Returns best medoids of all trials
			according to cost function.
	sample_ratio: float, default: .2
		It is used if the method is 'clara'
	clara_trials: int, default: 10,
		It is used if the method is 'clara'
	random_state : int, RandomState instance or None, optional, default: None
		If int, random_state is the seed used by the random number generator;
		If RandomState instance, random_state is the random number generator;
		If None, the random number generator is the RandomState instance used
		by `np.random`.

	Examples

	"""

	def __init__(self, n_clusters=3, max_iter=300, method='clara', sample_ratio=.2, clara_trials=10, random_state=0):
		self.n_clusters = n_clusters
		self.max_iter = max_iter
		self.method = method
		self.sample_ratio = sample_ratio
		self.clara_trials = clara_trials
		self.random_state = np.random.RandomState(random_state)
		self.best_medoids = []
		self.min_cost = float('inf')

	def fit(self, X):
		"""Compute k-medoids clustering. If method is 'pam'
		Parameters
		----------
		X : array-like, shape=(n_samples, n_features)
			Training instances to cluster.
		Returns
		----------
		self : MyKMedoids
		"""
		self.X = X
		#self.random_state = np.random.RandomState(self.random_state)

		if self.method == "clara":

			for i in range(self.clara_trials):

				samples = self.sample()
				best_medoids, min_cost = self.pam(samples)

				clusters = self.generate_clusters(best_medoids, X)			# generate clusters for the whole dataset

				cost = self.calculate_cost(best_medoids, clusters)			# calculate cost for the entire dataset
				
				if cost < self.min_cost:
					self.min_cost = cost
					self.best_medoids = best_medoids


		elif self.method == "pam":
			self.best_medoids, self.min_cost = self.pam(X)

		return self


	def sample(self):
		"""Samples from the data with given sample_ratio.

		Returns
		-------
		X : {array-like, sparse matrix}, shape = [n_samples, n_features]
			New data to transform.
		"""
		#random_state = np.random.RandomState(self.random_state)
		sample_indexes = self.random_state.permutation(self.X.shape[0])[:int(self.sample_ratio * self.X.shape[0])]
		
		samples = self.X[sample_indexes]

		return samples


	def pam(self, X):
		"""
		kMedoids - PAM
		See more : http://en.wikipedia.org/wiki/K-medoids
		The most common realisation of k-medoid clustering is the Partitioning Around Medoids (PAM) algorithm and is as follows:[2]
		1. Initialize: randomly select k of the n data points as the medoids"""
		shape = (0, X.shape[1])

		#random_state = np.random.RandomState(self.random_state)		
		medoid_indexes = self.random_state.permutation(X.shape[0])[:int(self.n_clusters)]      # select unique medoids randomly		

		medoids = X[medoid_indexes]

		sample_best_medoids = medoids.copy()		

		#2. Associate each data point to the closest medoid. ("closest" here is defined using any valid distance metric, most commonly Euclidean distance, Manhattan distance or Minkowski distance)
		clusters = self.generate_clusters(medoids, X)
		
		sample_min_cost = self.calculate_cost(medoids, clusters)

		"""3. For each medoid m
				For each non-medoid data point o
				  Swap m and o and compute the total cost of the configuration
		4. Select the configuration with the lowest cost.
		5. repeat steps 2 to 4 until there is no change in the medoid.

		"""
		
		iteration = 0
		stop = False

		for i, medoid in enumerate(medoids):
			for j, nonmedoid in enumerate(clusters[i]):

				new_medoids = medoids.copy()
				new_clusters = copy.deepcopy(clusters)

				new_medoids[i], new_clusters[i][j] = new_clusters[i][j].copy(), new_medoids[i].copy()

				new_cost = self.calculate_cost(new_medoids, new_clusters)
				
				if new_cost < sample_min_cost:
					sample_best_medoids = new_medoids.copy()
					sample_min_cost = new_cost


				iteration += 1

				if iteration >= 300:
					stop = True
					break
			
			if stop:
				break


		"""
		Parameters
		----------
		X : {array-like, sparse matrix}, shape = [n_samples, n_features]
			New data to transform.
		Returns
		-------
		best_medoids, min_cost : tuple, shape [n_samples,]
			Best medoids found and the cost according to it.
		"""
		return (sample_best_medoids, sample_min_cost)


	def generate_clusters(self, medoids, samples):
		"""Generates clusters according to distance to medoids. Order
		is same with given medoids array.
		Parameters
		----------
		medoids: array_like, shape = [n_clusters, n_features]
		samples: array-like, shape = [n_samples, n_features]
		Returns
		-------
		clusters : array-like, shape = [n_clusters, elements_inside_cluster, n_features]
		"""

		#shape = (self.n_clusters, samples.shape[0], samples.shape[1])
		clusters = []
		clusters_dict = {}
		
		for i in range(self.n_clusters):				# create a dictionary with keys: cluster numbers, values: empty lists
			clusters_dict.setdefault(i, [])

		for sample in samples:
			if not np.any((medoids[:] == sample).all(1)):			# if sample is not in medoids (to add non-meodids to clusters)
				
				distances_to_medoids = np.linalg.norm(sample - medoids, axis = 1)
				closest_medoid = np.argmin(distances_to_medoids)

				clusters_dict[closest_medoid].append(sample)		# add the sample to the closest cluster

		for i in range(self.n_clusters):				 # dict to list
			clusters.append(clusters_dict[i])

		return clusters
			

	def calculate_cost(self, medoids, clusters):
		"""Calculates cost of each medoid's cluster with squared euclidean function.
		Parameters
		----------
		medoids: array_like, shape = [n_clusters, n_features]
		clusters: array-like, shape = [n_clusters, elements_inside_cluster, n_features]
		Returns
		-------
		cost : float
			total cost of clusters
		"""
		total_cost = 0.0

		for i in range(len(clusters)):
			for j in range(len(clusters[i])):
				# find squared euclidean distance of each element to its cluster center(medoid)
				squared_distance = distance.sqeuclidean(clusters[i][j], medoids[i])
				total_cost += squared_distance

		return total_cost



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
		labels = []

		for i in range(len(X)):
			distances_to_medoids = np.linalg.norm(X[i] - self.best_medoids, axis = 1)
			labels.append(np.argmin(distances_to_medoids))

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
	X = np.array([np.array([2., 6.]),
				  np.array([3., 4.]),
				  np.array([3., 8.]),
				  np.array([4., 7.]),
				  np.array([6., 2.]),
				  np.array([6., 4.]),
				  np.array([7., 3.]),
				  np.array([7., 4.]),
				  np.array([8., 5.]),
				  np.array([7., 6.])
				  ])

	kmedoids = MyKMedoids(n_clusters=2, random_state=0, sample_ratio = 1)
	print kmedoids.fit_predict(X)
	# [1 1 1 1 0 0 0 0 0 0]
	print kmedoids.best_medoids
	# [array([7., 4.]), array([2., 6.])]
	print kmedoids.min_cost
	# 28.0


	#test = MyKMedoids(method="clara", n_clusters = 2, sample_ratio=0.5, clara_trials= 10, random_state= 0)
	#print test.fit_predict(X)
	#print "FINISH", test.best_medoids
	#print test.min_cost	