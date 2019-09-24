import sys
import numpy as np
import h5py
import tensorflow as tf
print(tf.test.is_gpu_available())

from keras.models import Sequential
from keras.models import model_from_json
from keras.layers import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers import Convolution1D, MaxPooling1D, ZeroPadding1D
from keras.utils import np_utils
from helper_cab import *


INPUTDATA = sys.argv[1]
FILELABEL = sys.argv[2]

ARCHITECTURE = "../models/dv_4c_10_4001_hetplusWT_talmojan_cnn_arch.json"
WEIGHTS = "../models/dv_4c_10_4001_hetplusWT_talmojan_cnn_weights.h5"


model = model_from_json(open(ARCHITECTURE).read())
model.load_weights(WEIGHTS)
model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])

f = h5py.File(INPUTDATA,'r') 
data = f.get('data') 
data = np.transpose(np.array(data)) # For converting to numpy array

L = 4001 # total length of each window
nb_hist = (L-1)/2 

datalen = data.shape[0]-nb_hist*2;

out_p=np.zeros((data.shape[0],4))

for i in range(0,10):
		print(i)
		print(data.shape)
		lim = round((data.shape[0]-L)/10)
		idx=np.arange(lim)+lim*i+nb_hist
		print(idx.shape)
		idx = np.intp(idx)
		idx = np.array(idx)
		print(idx.shape)
		print(idx[0])
		print(idx[-1])
		X = time_embed_predict(data, nb_hist, batch_idx=idx)
		X = X.reshape(X.shape[0], X.shape[1], 1)
		X = X.astype('float32')
		out_p[idx,:] = model.predict(X)
		print(out_p[-1,:])
		np.save(FILELABEL + "_out_p", out_p)