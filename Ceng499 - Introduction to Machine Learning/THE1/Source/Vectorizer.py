#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nltk
import numpy as np
import math
from collections import Counter

class Vectorizer:
	def __init__(self, min_word_length = 3, max_df = 1.0, min_df = 0.0):
		self.min_word_length = min_word_length
		self.max_df = max_df
		self.min_df = min_df
		self.term_df_dict = {}
		self.vocabulary = []
		self.term_idf_dict = {}

	def fit(self, raw_documents):
		#print " fitting started"
		"""Generates vocabulary for feature extraction. Ignores words shorter than min_word_length and document frequency
		not between max_df and min_df.

		:param raw_documents: list of string for creating vocabulary
		:return: None
		"""
		self.document_count = len(raw_documents)
		# TODO: Implement this method
		vocab = set()

		#doc_number = 0
		for document in raw_documents:
			words = document.split()            # get each word in the document as a list
			set_words = set(words)				# set contains unique words, eliminates the duplicates in the document

			for word in set_words:
				if len(word) >= self.min_word_length:
					if word not in vocab:
						vocab.add(word)
						self.term_df_dict[word] = 1.0
					else:
						self.term_df_dict[word] += 1.0

			#doc_number += 1
			#print "fitted DOC ", doc_number

		for word in vocab:			
			document_frequency = self.term_df_dict[word] / self.document_count     # find document frequency
			
			if document_frequency >= self.min_df and document_frequency <= self.max_df:
				self.vocabulary.append(word)
				idf = math.log((1.0 + self.document_count)/(1.0 + self.term_df_dict[word])) + 1.0          # find idf
				self.term_idf_dict[word] = idf

			
		#print "fitting DONE for all docs"


	def _transform(self, raw_document, method):
		"""Creates a feature vector for given raw_document according to vocabulary.

		:param raw_document: string
		:param method: one of count, existance, tf-idf
		:return: numpy array as feature vector
		"""
		# TODO: Implement this method
		size = len(self.vocabulary)		

		words = raw_document.split()        # get each word in the document as a list
		counter = Counter(words)			# store how many times each word appear in the document

		if method == 'count':
			feature_vector = np.zeros((size,), int)
			for i in range(size):
				feature_vector[i] = counter[self.vocabulary[i]]

			return feature_vector

		elif method == 'existance':
			feature_vector = np.zeros((size,), int)
			for i in range(size):
				if counter[self.vocabulary[i]] > 0:
					feature_vector[i] = 1

			return feature_vector

		elif method == 'tf-idf':
			feature_vector = np.zeros((size,), float)
			i = 0
			for word in self.vocabulary:                				# for each word in a document
				tf = counter[word]             							# find tf

				tf_idf = tf * self.term_idf_dict[word]                	# find tf-idf
				
				feature_vector[i] = tf_idf

				i += 1

			norm = np.linalg.norm(feature_vector)
			if norm > 0.0:
				normalized = feature_vector / norm
				return normalized
			else:
				return feature_vector


	def transform(self, raw_documents, method = "tf-idf"):
		#print "transform started"
		"""For each document in raw_documents calls _transform and returns array of arrays.

		:param raw_documents: list of string
		:param method: one of count, existance, tf-idf
		:return: numpy array of feature-vectors
		"""
		# TODO: Implement this method
		size = len(raw_documents)
		feature_vectors = np.empty((0, len(self.vocabulary)), int)

		for i in range(size):
			feature_vectors = np.append(feature_vectors, [self._transform(raw_documents[i], method)], axis = 0)
			#print "DOC ", i, " is transformed."

		#print "transform ended"
		return feature_vectors

	def fit_transform(self, raw_documents, method = "tf-idf"):
		"""Calls fit and transform methods respectively.

		:param raw_documents: list of string
		:param method: one of count, existance, tf-idf
		:return: numpy array of feature-vectors
		"""
		# TODO: Implement this method
		self.fit(raw_documents)
		return self.transform(raw_documents, method)

	def get_feature_names(self):
		"""Returns vocabulary.

		:return: list of string
		"""
		try:
			self.vocabulary
		except AttributeError:
			print "Please first fit the model."
			return []
		return self.vocabulary

	def get_term_dfs(self):
		"""Returns number of occurances for each term in the vocabulary in sorted.

		:return: array of tuples
		"""
		return sorted(self.term_df_dict.iteritems(), key = lambda (k, v): (v, k), reverse = True)

if __name__=="__main__":
	v = Vectorizer(min_df = 0.25, max_df = 0.75)
	contents = [
	 "this is the first document",
	 "this document is the second document",
	 "and this is the third one",
	 "is this the first document",
	]

	v.fit(contents)
	print v.get_feature_names()
	existance_vector = v.transform(contents, method = "existance")        
	print existance_vector
	count_vector = v.transform(contents, method = "count")        
	print count_vector
	tf_idf_vector = v.transform(contents, method = "tf-idf")
	print tf_idf_vector
