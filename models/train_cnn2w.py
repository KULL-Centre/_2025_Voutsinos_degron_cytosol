import os,csv,pickle,time
# os.environ['CUDA_VISIBLE_DEVICES'] = '2'

import numpy as np
import sklearn
import tensorflow as tf
print( tf.config.experimental.list_physical_devices() )
    
from tensorflow import keras
import keras.backend as K
from keras import layers
from keras import regularizers

@tf.autograph.experimental.do_not_convert
def pearson(y_true, y_pred):
    # assert len(y_true.shape) == 1
    assert y_true.shape == y_pred.shape
    cent_true = y_true - tf.math.reduce_mean(y_true, axis=0)
    cent_pred = y_pred - tf.math.reduce_mean(y_pred, axis=0)
    dot2_true = tf.math.reduce_sum(tf.math.multiply(cent_true,cent_true))
    dot2_pred = tf.math.reduce_sum(tf.math.multiply(cent_pred,cent_pred))
    return( tf.math.reduce_sum(cent_true*cent_pred)/tf.math.sqrt(dot2_true*dot2_pred) )

def calc_mse(y_true, y_pred):
    return( np.mean(np.square(y_true-y_pred)) )

    
aa2i = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9,
        'M':10, 'N':11, 'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19}

seqs = []
scores = []
input_filename = "../github/models/degron_train.csv"
with open(input_filename) as csv_file:
    ci_aa = 1
    ci_score = 6  # FACS transformed
    csv_reader = csv.reader(csv_file, delimiter=',')
    header = next(csv_reader)
    print("Reading file '%s' with sequence in column '%s' and score in column '%s'" % (input_filename,header[ci_aa],header[ci_score]))
    for row in csv_reader:
        assert len(row) == len(header)
        if all( [aa in aa2i.keys() for aa in row[ci_aa]]):
            # name, aa, n_rep, cells, score, std, score_facs, std_facs
            seqs.append(row[ci_aa])
            scores.append(row[ci_score])

# labels
scores = np.array(scores, dtype="float32")

def str2mat(s):
    n = len(s)
    v = np.zeros((n,len(aa2i)))
    for (i,aa) in enumerate(s):
        v[i,aa2i[aa] ] = 1
    return v

seqs_mat = np.array([str2mat(seq) for seq in seqs])

# Settings
settings = {}
settings['name']        = 'cnn2w1'
settings['cnn_win']     =    1
settings['filters_int'] =   25
settings['n_ct_res']    =    5
settings['filters_ct']  =   25
settings['dense_dim']   =   50
settings['activation']  = 'elu'
settings['ct_l1weight']   = None
settings['l2weight']      = 1e-6
settings['dropout_rate']  = 0.2
settings['batch_size']    = 1024
settings['learning_rate']    = 1e-4
settings['learn_rate_step3'] = 1e-5
settings['epochs']           = 2000
settings['stop_patience']    = 100
settings['validation_split'] = 0.2
print(settings)

t_begin = time.time()

# Step 1, fit internal layer
K.clear_session()
inputs    = layers.Input(shape=(None, len(aa2i)), dtype="float32")
internal  = layers.Conv1D(settings['filters_int'], settings['cnn_win'], strides=1, padding="valid", name='int_cnn', activation=settings['activation'],
                          kernel_regularizer=regularizers.l2(settings['l2weight']))(inputs)
internal  = layers.GlobalAveragePooling1D(data_format='channels_last', keepdims=False)(internal)
internal  = layers.Dropout(settings['dropout_rate'])(internal)
internal  = layers.Dense(settings['dense_dim'], activation=settings['activation'], name='int2', kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)
internal  = layers.Dropout(settings['dropout_rate'])(internal)
internal  = layers.Dense(settings['dense_dim'], activation=settings['activation'], name='int3', kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)
output    = layers.Dense(1, name='latent_int', activation=None, kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)

opt = keras.optimizers.legacy.Adam(learning_rate=settings['learning_rate'])
stop_callback = keras.callbacks.EarlyStopping(monitor='val_loss', verbose=1, start_from_epoch=1, patience=settings['stop_patience'])

cnn_int = keras.models.Model(inputs=inputs, outputs=output)
cnn_int.summary()
cnn_int.compile(loss='mse', optimizer=opt, metrics=['mse','mae',pearson])

# train internal weights on internal part only and fix after this
seqs_mat_int = seqs_mat[:,0:25,:]
fit_int = cnn_int.fit(x=seqs_mat_int, y=scores, epochs=settings['epochs'], batch_size=settings['batch_size'],
                      validation_split=settings['validation_split'], callbacks=[stop_callback], verbose=2)

seconds_int = time.time()-t_begin
epochs_int = len(fit_int.epoch)
val_mse_int = min(fit_int.history['val_mse'])
val_rp_int = max(fit_int.history['val_pearson'])

# predict all data considering all positions for later comparison
pred_int1 = cnn_int.predict(seqs_mat, verbose=2).flatten()
mse_int1 = calc_mse(scores, pred_int1)
rp_int1 = pearson(scores, pred_int1)
print("Done step 1: Time %d s, epochs %d, best validation MSE %.4f Pearson %.4f, full data MSE %.4f Pearson %.4f" %
      (seconds_int,epochs_int,val_mse_int,val_rp_int,mse_int1,rp_int1))

# Step 2, add CT layer and fix parameters of internal model
K.clear_session()
input_int = layers.Input(shape=(None, len(aa2i)), name='input_int', dtype="float32")
input_ct  = layers.Input(shape=(settings['n_ct_res']*len(aa2i),), name='input_ct', dtype="float32")

internal  = layers.Conv1D(settings['filters_int'], settings['cnn_win'], strides=1, padding='valid', name='int_cnn', activation=settings['activation'],
                          trainable=False)(input_int)
internal  = layers.GlobalAveragePooling1D(data_format='channels_last', keepdims=False)(internal)
# internal  = layers.GlobalMaxPooling1D(data_format='channels_last', keepdims=False)(internal)
internal  = layers.Dense(settings['dense_dim'], name='int2', activation=settings['activation'], trainable=False)(internal)
internal  = layers.Dense(settings['dense_dim'], name='int3', activation=settings['activation'], trainable=False)(internal)
internal  = layers.Dense(1, name='latent_int', activation=None, trainable=False)(internal)

cterm = layers.Dense(settings['filters_ct'], name='ct1', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(input_ct)
# cterm = layers.Dense(settings['filters_ct'], name='ct1', activation=settings['activation'], kernel_regularizer=regularizers.l1(settings['ct_l1weight']))(input_ct)
cterm     = layers.Dropout(settings['dropout_rate'])(cterm)
cterm     = layers.Dense(settings['dense_dim'], name='ct2', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)
cterm     = layers.Dropout(settings['dropout_rate'])(cterm)
cterm     = layers.Dense(settings['dense_dim'], name='ct3', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)
# cterm     = layers.Dropout(settings['dropout_rate'])(cterm)
cterm     = layers.Dense(1, name='latent_ct', activation=None, kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)

output    = layers.Add()([internal,cterm])
cnn_ct    = keras.models.Model(inputs=[input_int,input_ct], outputs=output)
cnn_ct.summary()

for node_name in ['int_cnn','int2','int3','latent_int']:
    int_weights = cnn_int.get_layer(node_name).weights
    cnn_ct.get_layer(node_name).set_weights(int_weights)
    
cnn_ct.compile(loss='mse', optimizer=opt, metrics=['mse','mae',pearson])
seqs_mat_ct = seqs_mat[:,(30-settings['n_ct_res']):30,:]
seqs_oh_ct = seqs_mat_ct.reshape(*seqs_mat_ct.shape[:-2], -1)
fit_ct = cnn_ct.fit(x=[seqs_mat,seqs_oh_ct], y=scores, epochs=settings['epochs'], batch_size=settings['batch_size'],
                    validation_split=settings['validation_split'], callbacks=[stop_callback], verbose=2)

assert len(cnn_ct.get_layer('ct1').get_weights()[0][0]) == settings['filters_ct']
cut = 1e-5
ct_nz_weights = [np.sum(np.abs(cnn_ct.get_layer('ct1').get_weights()[0][:,i]) > cut) for i in range(settings['filters_ct'])]
print("Number of weights greater than %.1g in each of %d C-degron filters (mean %.2f)" % (cut, settings['filters_ct'], np.mean(ct_nz_weights)))
print(ct_nz_weights)

seconds_ct = time.time()-t_begin
epochs_ct = len(fit_ct.epoch)
val_mse_ct = min(fit_ct.history['val_mse'])
val_rp_ct = max(fit_ct.history['val_pearson'])

# predict all data
pred2 = cnn_ct.predict([seqs_mat,seqs_oh_ct], verbose=2).flatten()
mse2 = calc_mse(scores, pred2)
rp2 = pearson(scores, pred2)
print("Done step 2: Time %d s, epochs %d, MSE %.4f Pearson %.4f, best validation MSE %.4f Pearson %.4f" % (seconds_ct,epochs_ct,mse2,rp2,val_mse_ct,val_rp_ct))

# Step 3, refine 
K.clear_session()
input_int = layers.Input(shape=(None, len(aa2i)), name='input_int', dtype="float32")
input_ct  = layers.Input(shape=(settings['n_ct_res']*len(aa2i),), name='input_ct', dtype="float32")

internal  = layers.Conv1D(settings['filters_int'], settings['cnn_win'], strides=1, padding='valid', name='int_cnn', activation=settings['activation'],
                          trainable=False)(input_int)
internal  = layers.GlobalAveragePooling1D(data_format='channels_last', keepdims=False)(internal)
internal  = layers.Dropout(settings['dropout_rate'])(internal)
internal  = layers.Dense(settings['dense_dim'], name='int2', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)
internal  = layers.Dropout(settings['dropout_rate'])(internal)
internal  = layers.Dense(settings['dense_dim'], name='int3', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)
internal  = layers.Dense(1, name='latent_int', activation=None, kernel_regularizer=regularizers.l2(settings['l2weight']))(internal)

cterm     = layers.Dense(settings['filters_ct'], name='ct1', activation=settings['activation'], trainable=False)(input_ct)
cterm     = layers.Dropout(settings['dropout_rate'])(cterm)
cterm     = layers.Dense(settings['dense_dim'], name='ct2', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)
cterm     = layers.Dropout(settings['dropout_rate'])(cterm)
cterm     = layers.Dense(settings['dense_dim'], name='ct3', activation=settings['activation'], kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)
cterm     = layers.Dense(1, name='latent_ct', activation=None, kernel_regularizer=regularizers.l2(settings['l2weight']))(cterm)

output    = layers.Add()([internal,cterm])
cnn       = keras.models.Model(inputs=[input_int,input_ct], outputs=output)
cnn.summary()

for node_name in ['int_cnn','int2','int3','ct1','ct2','ct3','latent_int','latent_ct']:
    ct_weights = cnn_ct.get_layer(node_name).weights
    cnn.get_layer(node_name).set_weights(ct_weights)

opt_step3 = keras.optimizers.legacy.Adam(learning_rate=settings['learn_rate_step3'])

cnn.compile(loss='mse', optimizer=opt_step3, metrics=['mse','mae',pearson])
fit = cnn.fit(x=[seqs_mat,seqs_oh_ct], y=scores, epochs=settings['epochs'], batch_size=settings['batch_size'],
              validation_split=settings['validation_split'], callbacks=[stop_callback], verbose=2)

t_end = time.time()

# predict all data
pred = cnn.predict([seqs_mat,seqs_oh_ct], verbose=2).flatten()
int_model = keras.models.Model(inputs=cnn.get_layer('input_int').input, outputs=cnn.get_layer('latent_int').output)
pred_int = int_model.predict(seqs_mat, verbose=2).flatten()

seconds = t_end-t_begin
epochs = len(fit.epoch)

# evaluate full data
mse = calc_mse(scores, pred)
rp = pearson(scores, pred)
mse_int = calc_mse(scores, pred_int)
rp_int = pearson(scores, pred_int)

val_mse = min(fit.history['val_mse'])
val_rp = max(fit.history['val_pearson'])

print("Done training model: Time %d s, epochs %d, MSE %.4f Pearson %.4f, best validation MSE %.4f Pearson %.4f" % (seconds,epochs,mse,rp,val_mse,val_rp))
print("Evaluation as internal tiles MSE %.4f vs %.4f after step 1 (diff %.5f)" % (mse_int,mse_int1,mse_int-mse_int1))

# store it
cnn.save("%s_int.tf" % (settings['name']), save_format='tf')
cnn.save_weights("%s.keras" % (settings['name']), save_format='keras')
cnn_ct.save_weights("%s_step2.keras" % (settings['name']), save_format='keras')
with open("%s.pickle" % (settings['name']), 'wb') as pickle_file:
    pickle.dump((fit.history,pred,pred_int,settings), pickle_file)

print("Training done")
