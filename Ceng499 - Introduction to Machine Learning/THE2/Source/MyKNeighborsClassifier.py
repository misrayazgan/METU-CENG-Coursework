#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.spatial import distance

class MyKNeighborsClassifier:
	"""Classifier implementing the k-nearest neighbors vote similar to sklearn 
	library but different.
	https://goo.gl/Cmji3U

	But still same.
	
	Parameters
	----------
	n_neighbors : int, optional (default = 5)
		Number of neighbors to use by default.
	method : string, optional (default = 'classical')
		method for voting. Possible values:
		- 'classical' : uniform weights.  All points in each neighborhood
		  are weighted equally.
		- 'weighted' : weight points by the inverse of their distance.
		  in this case, closer neighbors of a query point will have a
		  greater influence than neighbors which are further away.
		- 'validity' weights are calculated with distance and multiplied
		  of validity for each voter.  
		Note: implementing kd_tree is bonus.
	norm : {'l1', 'l2'}, optional (default = 'l2')
		Distance norm. 'l1' is manhattan distance. 'l2' is euclidean distance.
	Examples
	--------
	"""
	def __init__(self, n_neighbors=5, method='classical', norm='l2'):
		self.n_neighbors = n_neighbors
		self.method = method
		self.norm = norm
		self.labels = []
		self.distinct_labels = []
		self.data = []

	def fit(self, X, y):
		"""Fit the model using X as training data and y as target values
		Parameters
		----------
		X : array-like, shape (n_query, n_features),
			Training data. 
		y : array-like, shape = [n_samples] 
			Target values.
		"""

		self.data = X
		self.labels = y

		self.distinct_labels = sorted(list(set(y)))					# get distinct classes

		if self.method == "validity":
			self.validities = []

			# calculate validities for each sample			
			for sample_i, sample in enumerate(self.data):
				distances = []

				if self.norm == "l1":
					# apply Manhattan distance
					for i, x in enumerate(self.data):
						if x != sample:
							print x, sample
							distances.append([distance.cityblock(x, sample), self.labels[i]])			# store [distance, label] pairs in distances list

				elif self.norm == "l2":
					# apply Euclidean distance
					for i, x in enumerate(self.data):
						if x != sample:
							print x, sample
							distances.append([distance.euclidean(x, sample), self.labels[i]])			# store [distance, label] pairs in distances list

				# get n nearest neighbors
				nearest_neighbors = sorted(distances, key = lambda x: x[0])[:self.n_neighbors]		# sort wrt distance

				label_weights = [0.0] * len(self.distinct_labels)						# store label weights wrt neighbors for each label

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					label_weights[label_index] += (1.0 / (neighbor[0] + 1e-15))

				sample_label = self.labels[sample_i]
				validity = label_weights[self.distinct_labels.index(sample_label)] / sum(label_weights)

				self.validities.append(validity)


	def predict(self, X):
		"""Predict the class labels for the provided data
		Parameters
		----------
		X : array-like, shape (n_query, n_features),
			Test samples.
		Returns
		-------
		y : array of shape [n_samples]
			Class labels for each data sample.
		"""
		if len(self.labels) == 0:
			raise ValueError("You should fit first!")
			
		y = []			# store labels for each data

		for x_i, x in enumerate(X):
			# find distance of x to all training data
			distances = []

			if self.norm == "l1":				
				for i, data in enumerate(self.data):
					distances.append([distance.cityblock(x, data), self.labels[i], i])			# store [distance, label, data_index] pairs in distances list

			elif self.norm == "l2":
				for i, data in enumerate(self.data):
					distances.append([distance.euclidean(x, data), self.labels[i], i])			# store [distance, label, data_index] pairs in distances list

			# get n nearest neighbors
			nearest_neighbors = sorted(distances, key = lambda x: x[0])[:self.n_neighbors]		# sort wrt distance
				

			votes = [0] * self.n_neighbors 							# store vote(label) of each nearest neighbor

			if self.method == "classical":

				for i, neighbor in enumerate(nearest_neighbors):
					votes[i] = neighbor[1]

				y.append(max(votes, key = votes.count))

			elif self.method == "weighted":
					
				weights = []										# store weight of each neighbor (1/(distance + 1))

				for i, neighbor in enumerate(nearest_neighbors):
					weights.append(1.0 / (neighbor[0] + 1e-15))

				total_weights = [0] * len(self.distinct_labels)					# store total weights for each label wrt indexes in self.distinct_labels

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					total_weights[label_index] += weights[i]

				y.append(self.distinct_labels[total_weights.index(max(total_weights))])

			elif self.method == "validity":

				validities = []										# store weight of each neighbor (1/(distance + 1))

				for i, neighbor in enumerate(nearest_neighbors):
					validities.append((1.0 / (neighbor[0] + 1e-15)) * self.validities[neighbor[2]])						
				
				total_validities = [0] * len(self.distinct_labels)					# store total validity*weight values for each label wrt indexes in self.distinct_labels

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					total_validities[label_index] += validities[i]

				y.append(self.distinct_labels[total_validities.index(max(total_validities))])

		return y

		
	def predict_proba(self, X, method=None):
		"""Return probability estimates for the test data X.
		Parameters
		----------
		X : array-like, shape (n_query, n_features),
			Test samples.
		method : string, if None uses self.method.
		Returns
		-------
		p : array of shape = [n_samples, n_classes]
			The class probabilities of the input samples. Classes are ordered
			by lexicographic order.
		"""
		if method == None:
			method = self.method

		p = [[0] * len(self.distinct_labels)] * len(X)					# store probabilities of each class for each sample
		
		for x_i, x in enumerate(X):			
			# find distance of x to all training data
			distances = []

			if self.norm == "l1":				
				for i, data in enumerate(self.data):
					distances.append([distance.cityblock(x, data), self.labels[i], i])			# store [distance, label, data_index] pairs in distances list

			elif self.norm == "l2":
				for i, data in enumerate(self.data):
					distances.append([distance.euclidean(x, data), self.labels[i], i])			# store [distance, label, data_index] pairs in distances list

			# get n nearest neighbors
			nearest_neighbors = sorted(distances, key = lambda x: x[0])[:self.n_neighbors]		# sort wrt distance
				

			classes = [0] * len(self.distinct_labels) 					# store how many neighbors there are in each class
			
			if method == "classical":

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					classes[label_index] += 1
				
				for class_i in range(len(classes)):
					p[x_i][class_i] = float(classes[class_i]) / self.n_neighbors
				

			elif method == "weighted":

				weights = []										# store weight of each neighbor (1/(distance + 1))

				for i, neighbor in enumerate(nearest_neighbors):
					weights.append(1.0 / (neighbor[0] + 1e-15))

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					classes[label_index] += weights[i]

				for class_i in range(len(classes)):
					p[x_i][class_i] = float(classes[class_i]) / self.n_neighbors
				

			elif method == "validity":

				validities = []										# store validity*weight of each neighbor (1/(distance + 1))*validity

				for i, neighbor in enumerate(nearest_neighbors):
					validities.append((1.0 / (neighbor[0] + 1e-15)) * self.validities[neighbor[2]])						

				class_validities = [0] * len(self.distinct_labels)					# store total validity*weight values for each label wrt indexes in self.distinct_labels

				for i, neighbor in enumerate(nearest_neighbors):
					label_index = self.distinct_labels.index(neighbor[1])
					class_validities[label_index] += validities[i]

				for valid_i in range(len(class_validities)):
					p[x_i][valid_i] = float(class_validities[valid_i]) / self.n_neighbors

			# normalize probability list for each sample
			norm = sum(p[x_i])			
			p[x_i] = [round(j / norm, 8) for j in p[x_i]]

		return p


if __name__=='__main__':
	X = [[0], [1], [2], [3]]
	y = [0, 0, 1, 1]
	neigh = MyKNeighborsClassifier(n_neighbors=3, method="validity")
	neigh.fit(X, y)
	
	print neigh.predict([[1.1]]) #, [6.7], [5], [1.9], [0]])

	n = 0.9
	print(neigh.predict_proba([[n]], method='classical'))
	# [[0.66666667 0.33333333]]
	print(neigh.predict_proba([[n]], method='weighted'))
	# [[0.92436975 0.07563025]]
	print(neigh.predict_proba([[n]], method='validity'))
	# [[0.92682927 0.07317073]]