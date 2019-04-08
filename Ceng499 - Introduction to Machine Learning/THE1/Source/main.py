from ExperimentSuite import ExperimentSuite
from Preprocessor import Preprocessor
from Vectorizer import Vectorizer
import tensorflow as tf

DEFAULT_EPOCH = 75
DEFAULT_LAYERS = (512,)
DEFAULT_ACTIVATION = tf.nn.relu
DEFAULT_LOSS = "categorical_hinge" # "categorical_crossentropy", categorical_hinge


if __name__ == "__main__":
    # TODO: Make your experiments here
    
    es = ExperimentSuite()
    
    #pre = Preprocessor()
    #pre.preprocess()
    #print "preprocess is DONE"
    vec = Vectorizer(max_df = 0.97, min_df = 0.1)
    #print "you called the Vectorizer, starting to fit"

    # only fit with train contents
    vec.fit(es.train_contents)

    train_x = vec.transform(raw_documents = es.train_contents, method = 'tf-idf')
    #print train_x
    test_x = vec.transform(raw_documents = es.test_contents, method = 'tf-idf')
    #print test_x

    #with open("testsave.npy", 'w') as outf:
    #	np.save(outf, train_x)
    #with open("testsave.npy", 'r') as inpf:
    #	loaded = np.load(inpf)

    #print "fit_transform is DONE"

    callback = tf.keras.callbacks.TensorBoard(log_dir='callback_logs', histogram_freq = 0, write_graph = True, write_images = True)

    #print "training started"
    es.train_model(layers = DEFAULT_LAYERS, train_x = train_x, train_y = es.train_y, 
    	test_x = test_x, test_y = es.test_y, tbCallBack = callback, 
    	epoch = DEFAULT_EPOCH, activation = DEFAULT_ACTIVATION, loss = DEFAULT_LOSS)
    #print "training ended"
