#!/usr/bin/env python
"""
Program: train_esym_dnn.py
         Train the DNN that maps neutron-star data to Esym(rho)
         (the training of MR_EOS_v2.ipynb, as a script)

         python train_esym_dnn.py --input MR        # M-R input
         python train_esym_dnn.py --input MLambda   # M-Lambda input

         Network and training as in the paper: 12 dense layers of 100
         (ReLU, linear input and output layers), Adam with AMSgrad,
         learning rate 0.003, MSE loss, batch size 500, at most 2000
         epochs; here with early stopping on the validation loss.

         Output (in --out-dir, default models/), for tag = esym_<input>:
         <tag>.keras           model (best epoch)
         <tag>.weights.h5      weights
         <tag>_scaler.npz      input scaler (fitted on the training set)
         <tag>_test.npz        test set: raw and scaled inputs, outputs,
                               predictions, indices, EOS parameters
         <tag>_history.csv     loss per epoch
         <tag>_metrics.json    test errors and run settings
"""
import os
import sys
import json
import time
import argparse
import numpy as np

import esym_data as ed


def build_model(lr=0.003):
    """The network of MR_EOS_v2.ipynb (build_model)"""
    from tensorflow import keras
    model = keras.Sequential()
    model.add(keras.Input(shape=(100,)))
    model.add(keras.layers.Dense(100, activation='linear'))
    for _ in range(10):
        model.add(keras.layers.Dense(100, activation='relu'))
    model.add(keras.layers.Dense(100, activation='linear'))
    opt = keras.optimizers.Adam(learning_rate=lr, beta_1=0.9, beta_2=0.999, epsilon=1e-08, amsgrad=True)
    model.compile(optimizer=opt, loss='mse', metrics=['mae'])
    return model


def test_metrics(y_test, pred):
    """Absolute errors at 5 rho0 (last output point, as in the paper) and
    over all densities"""
    err = np.abs(y_test - pred)
    return {'mae_5rho0': float(np.mean(err[:, -1])),
            'std_5rho0': float(np.std(err[:, -1])),
            'mae_all': float(np.mean(err)),
            'max_5rho0': float(np.max(err[:, -1])),
            'mae_vs_rho': [float(v) for v in np.mean(err, axis=0)]}


def train(kind='MR', files=ed.NS_FILES, n_train=40000, epochs=2000, patience=100,
          batch_size=500, lr=0.003, seed=42, out_dir='models', tag=None, verbose=2):
    """Build the data sets, train with early stopping, save everything.
    Returns (model, sets, metrics)."""
    import tensorflow as tf
    from tensorflow import keras

    tag = tag or 'esym_' + kind
    os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(out_dir, tag)

    # --- reproducibility ---
    keras.utils.set_random_seed(seed)
    tf.config.experimental.enable_op_determinism()

    # --- data ---
    t0 = time.time()
    sets = ed.build_sets(kind, files, n_train, seed)
    d = sets['data']
    print('%s: %d sequences, %d pass the quality filter; train %d, validation %d, test %d'
          % (kind, d['n_total'], d['n_kept'], len(sets['idx_train']), len(sets['idx_val']),
             len(sets['idx_test'])), flush=True)
    sets['scaler'].save(base + '_scaler.npz')

    # --- model and training ---
    model = build_model(lr)
    callbacks = [
        keras.callbacks.EarlyStopping(monitor='val_loss', patience=patience,
                                      restore_best_weights=True, verbose=1),
        keras.callbacks.ModelCheckpoint(base + '_best.keras', monitor='val_loss',
                                        save_best_only=True),
        keras.callbacks.CSVLogger(base + '_history.csv'),
    ]
    f32 = lambda a: a.astype('float32')
    t1 = time.time()
    hist = model.fit(x=f32(sets['x_train']), y=f32(sets['y_train']),
                     epochs=epochs, batch_size=batch_size,
                     validation_data=(f32(sets['x_val']), f32(sets['y_val'])),
                     callbacks=callbacks, verbose=verbose)
    t_train = time.time() - t1

    # --- save the model (best weights restored) ---
    model.save(base + '.keras')
    model.save_weights(base + '.weights.h5')
    if os.path.exists(base + '_best.keras'):
        os.remove(base + '_best.keras')

    # --- test set ---
    pred = model.predict(f32(sets['x_test']), batch_size=batch_size, verbose=0)
    met = test_metrics(sets['y_test'], pred)
    vl = hist.history['val_loss']
    met.update({'input': kind, 'epochs_run': len(vl), 'best_epoch': int(np.argmin(vl)) + 1,
                'best_val_loss': float(np.min(vl)), 'patience': patience, 'max_epochs': epochs,
                'batch_size': batch_size, 'learning_rate': lr, 'seed': seed,
                'n_total': d['n_total'], 'n_kept': d['n_kept'], 'n_train': len(sets['idx_train']),
                'n_val': len(sets['idx_val']), 'n_test': len(sets['idx_test']),
                'files': list(files), 'train_time_s': t_train, 'total_time_s': time.time() - t0,
                'tensorflow': tf.__version__,
                'gpu': [x.name for x in tf.config.list_physical_devices('GPU')]})
    with open(base + '_metrics.json', 'w') as fo:
        json.dump(met, fo, indent=1)
    np.savez(base + '_test.npz', x_test_raw=sets['x_test_raw'], x_test=sets['x_test'],
             y_test=sets['y_test'], pred=pred, idx_test=sets['idx_test'],
             L=d['L'][sets['idx_test']], Ksym=d['Ksym'][sets['idx_test']],
             Jsym=d['Jsym'][sets['idx_test']], rhox=ed.RHOX)

    print('%s: stopped after %d epochs (best epoch %d, val_loss %.4g), %.0f s' %
          (kind, met['epochs_run'], met['best_epoch'], met['best_val_loss'], t_train))
    print('%s: test error at 5 rho0: %.2f +- %.2f MeV (max %.2f); all densities: %.2f MeV' %
          (kind, met['mae_5rho0'], met['std_5rho0'], met['max_5rho0'], met['mae_all']))
    return model, sets, met


def main():
    p = argparse.ArgumentParser(description='Train the Esym(rho) DNN on M-R or M-Lambda input')
    p.add_argument('--input', choices=['MR', 'MLambda'], required=True,
                   help='DNN input: M-R or M-Lambda sequences')
    p.add_argument('--ns-files', nargs='+', default=ed.NS_FILES,
                   help='NS files from tov_ml_cli.x (default: %s)' % ' '.join(ed.NS_FILES))
    p.add_argument('--n-train', type=int, default=40000, help='training set size (default: 40000)')
    p.add_argument('--epochs', type=int, default=2000, help='maximum number of epochs (default: 2000)')
    p.add_argument('--patience', type=int, default=100,
                   help='early stopping: epochs without improvement of val_loss (default: 100)')
    p.add_argument('--batch-size', type=int, default=500, help='batch size (default: 500)')
    p.add_argument('--lr', type=float, default=0.003, help='learning rate (default: 0.003)')
    p.add_argument('--seed', type=int, default=42, help='random seed (default: 42)')
    p.add_argument('--out-dir', default='models', help='output directory (default: models)')
    p.add_argument('--tag', default=None, help='file name prefix (default: esym_<input>)')
    a = p.parse_args()
    train(a.input, a.ns_files, a.n_train, a.epochs, a.patience, a.batch_size, a.lr,
          a.seed, a.out_dir, a.tag)


if __name__ == '__main__':
    main()
