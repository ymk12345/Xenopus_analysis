import csv
import numpy as np
import h5py
import math
import random
from helper_cab import *

# 1. Set up
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers import Convolution1D, MaxPooling1D, ZeroPadding1D
from keras.utils import np_utils
from keras.utils.io_utils import HDF5Matrix

##### Training/Model parameters

batch_size = 128 #16 or 32 recommended (smaller for less data; smaller batches give "better" results; larger batches make training faster)... number of training examples to consider at once
nb_epoch = 6 # # of times to iterate over the entire training set
nb_classes = 4 #nosong, female, male, overlap

# number of convolutional filters to use
nb_filters = 16 #less_filters #32 usual
# size of pooling area for max pooling
nb_pool = 2  # 2x2 pooling
# convolution kernel size
nb_conv = 9  
nb_pad = 4
#nb_hist = np.intp(1000 / 5)
nb_samples = 18997050;
#samples_per_epoch = np.intp(np.floor(nb_samples / batch_size) * batch_size)


L = 4001 # total length of each window
train_prop = .1
nb_hist = (L-1)/2 # # of samples before and after each timepoint to include in training


f = h5py.File('../preprocessed_data/preprocessed_fortraining_4001_hetplusWT.mat','r') 
outputModelName = 'dv_4c_10_4001_hetplusWT_talmojan_cnn'

data = f.get('data') 
classes = f.get('class')
data = np.transpose(np.array(data))
classes = np.transpose(np.array(classes))

print (data.shape)
print (classes.shape)

labels = np_utils.to_categorical(classes, nb_classes)
labels = np.reshape(labels, (labels.shape[0], labels.shape[1]))

#nb_samples = np.intp(np.floor(data.shape[0] / batch_size) * batch_size)
samples_per_epoch = np.intp(np.floor(data.shape[0] * train_prop / batch_size) * batch_size)
# samples_per_epoch = np.intp(np.floor(10000 / batch_size) * batch_size)

##### Calculating weights for training
class_weights = {0:1, 1:1, 2:1, 3:1}
print (class_weights)

#### Build Keras convolutional net architecture

model = Sequential()

#model.add(ZeroPadding1D(padding=nb_pad, input_shape=(L, 1)))
model.add(Convolution1D(nb_filters, nb_conv, border_mode='valid', subsample_length=2, input_shape=(L,1)))
model.add(Activation('relu')) 
#model.add(Dropout(0.25))
print(model.output_shape)
#model.add(ZeroPadding1D(padding=nb_pad))
model.add(Convolution1D(nb_filters, nb_conv, subsample_length=2))
model.add(Activation('relu'))
print(model.output_shape)

model.add(MaxPooling1D(pool_length=nb_pool))#model.output_shape[1]
print(model.output_shape)

#model.add(Dropout(0.25))
#model.add(ZeroPadding1D(padding=nb_pad))
model.add(Convolution1D(nb_filters, nb_conv, subsample_length=2))
model.add(Activation('relu'))
print(model.output_shape)

#model.add(MaxPooling1D(pool_length=nb_pool))
#print(model.output_shape)

model.add(Convolution1D(nb_filters, 15))
print(model.output_shape)

model.add(Activation('relu'))
model.add(Flatten())
print(model.output_shape)
model.add(Dense(128)) #64 smaller
model.add(Activation('relu'))
print(model.output_shape)
model.add(Dropout(0.4))
model.add(Dense(nb_classes))
model.add(Activation('softmax'))
print(model.output_shape)

model.compile(loss='categorical_crossentropy',
              optimizer='adam',
              metrics=['accuracy'])

#data_gen = generate(data, labels, nb_classes, nb_hist, samples_per_epoch, batch_size=32, train_size=0.6, shuffle=True, balance=False)

data_gen = generate_cab(data, labels, nb_classes, nb_hist, nb_epoch, batch_size, samples_per_epoch, shuffle=True, balance=False)

model.fit_generator(data_gen, samples_per_epoch=samples_per_epoch, nb_epoch=nb_epoch, verbose=1, class_weight=class_weights)

#model.fit(X_train, Y_train, class_weight=class_weights, 
#		  batch_size=batch_size, nb_epoch=nb_epoch,
#          verbose=1, validation_data=(X_test, Y_test))
#score = model.evaluate(X_test, Y_test, verbose=0)

#save_model(model, '../models/' + outputModelName + )

#print(score)
#print('Test score:', score[0])
#print('Test accuracy:', score[1])

save_model(model, '../models/' + outputModelName)
model.save_weights('../models/' + outputModelName + '_weights.h5')

json_string = model.to_json()
open('../models/' + outputModelName + '_architecture.json', 'w').write(json_string)
#model.save_weights('../models/' + outputModelName + '_weights.h5')

#X, resp = time_embed(data, labels, nb_hist)
#pred = model.predict(X)

#with h5py.File('../predict/1001_talmojan_pred.h5', 'w') as f:
#    f.create_dataset("pred", data=pred, compression="gzip")
#    f.create_dataset("resp", data=resp, compression="gzip")

