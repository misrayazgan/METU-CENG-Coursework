#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nltk
import os
from nltk.corpus import stopwords
import codecs
import errno
import string
import re


class Preprocessor:
    def __init__(self, dataset_directory="Dataset", processed_dataset_directory= "ProcessedDataset"):
        self.dataset_directory = dataset_directory
        self.processed_dataset_directory=processed_dataset_directory
        nltk.download("stopwords")
        nltk.download("punkt")
        self.stop_words = set(stopwords.words('english'))

    def _remove_puncs_numbers_stop_words(self, tokens):
        """Remove punctuations in the words, words including numbers and words in the stop_words list.

        :param tokens: list of string
        :return: list of string with cleaned version
        """
        # TODO: Implement this method
        numbers_removed = [token for token in tokens if not any(c.isdigit() for c in token)]
        punc_removed = []
        for word in numbers_removed:
            new_word = re.sub(r'[^\w\s]', '', word)
            #new_word = word.translate(None, string.punctuation)
            #tokenizer = nltk.tokenize.RegexpTokenizer(r'\w+')
            #new_word = tokenizer.tokenize(word)                 # returns list such as ['greetings']
            if len(new_word) != 0:
                punc_removed.append(new_word)
        	#    punc_removed = punc_removed + new_word          # merge new word with the previous ones

        result = [word for word in punc_removed if word not in self.stop_words]

        return result

    def _tokenize(self, sentence):
        """Tokenizes given string.

        :param sentence: string to tokenize
        :return: list of string with tokens
        """
        # TODO: Implement this method
        tokens = nltk.word_tokenize(sentence.lower())        

        return tokens

    def _stem(self, tokens):
        """Stems the tokens with nltk SnowballStemmer

        :param tokens: list of string
        :return: list of string with words stems
        """
        # TODO: Implement this method
        result = []
        stemmer = nltk.SnowballStemmer('english')
        for token in tokens:
        	token = stemmer.stem(token)
        	if token != '':
        		result.append(token) #.encode('UTF-8', 'strict'))
        
        return result

    def preprocess_document(self, document):
        """Calls methods _tokenize, _remove_puncs_numbers_stop_words and _stem respectively.

        :param document: string to preprocess
        :return: string with processed version
        """
        # TODO: Implement this method
        tokens = self._tokenize(document)
        words = self._remove_puncs_numbers_stop_words(tokens)
        stemmed_words = self._stem(words)
        result_string = ""
        for word in stemmed_words:
            result_string = result_string + word + ' '
        result_string = result_string[:-1]          # remove the last ' ' from the string
        
        return result_string

    def preprocess(self):
        """Walks through the given directory and calls preprocess_document method. The output is
        persisted into processed_dataset_directory by keeping directory structure.

        :return: None
        """
        for root, dirs, files in os.walk(self.dataset_directory):
            if os.path.basename(root) != self.dataset_directory:
                print "Processing", root, "directory."
                dest_dir = self.processed_dataset_directory+"/"+root.lstrip(self.dataset_directory+"/")
                if not os.path.exists(dest_dir):
                    try:
                        os.makedirs(dest_dir)
                    except OSError as exc:
                        if exc.errno != errno.EEXIST:
                            raise
                for file in files:
                    file_path = root + "/" + file
                    with codecs.open(file_path, "r", "ISO-8859-1") as f:
                        data = f.read().replace("\n", " ")
                    processed_data = self.preprocess_document(data)
                    output_file_path = dest_dir + "/" + file
                    with codecs.open(output_file_path, "w", "ISO-8859-1") as o:
                        o.write(processed_data)

if __name__=="__main__":
    text =  """ Greetings,
                shall I sit or stand?
                - Tell us.
                - Tell us.
                I'll tell. We bought the goods
                from Black Faik.
                We reloaded the truck
                in Karabuk.
                I was driving the truck
                till Adana.
                - What are you talking about?
                - And you?!
                You've abducted me,
                you'll do the talking.
                I'm confused anyway.
                - Aggressive.
                - Aggressive.
                Yeah, aggressive.
                Is that it?"""

    text2 = "I can't bridge. No...don't kill me"
    p = Preprocessor()
    '''token =  p._tokenize(text)
                print token
                removed =  p._remove_puncs_numbers_stop_words(token)
                print removed
                stemmed =  p._stem(removed)
                print stemmed
                print p.preprocess_document(text)
                '''
    print p.preprocess_document(text)
    
