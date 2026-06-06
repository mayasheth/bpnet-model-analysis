import tensorflow as tf
import bpnet.model.arch
import json
import os

folds = ["fold0", "fold1", "fold2", "fold3", "fold4"]
out_dir = "/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models_h5"
os.makedirs(out_dir, exist_ok = True)

keras_model_paths = [os.path.join("/oak/stanford/groups/engreitz/Users/sheth/EP300_BPNet/2025_0517_official_EP300_K562_model/models/release_run_1", f, "ENCSR000EGE/ENCSR000EGE_split000")
    for f in folds]
h5_models_out = [os.path.join(out_dir, f, "ENCSR000EGE_model.h5") for f in folds]

for f, model_in, model_out in zip(folds, keras_model_paths, h5_models_out):
    print(f)
    model = tf.keras.models.load_model(model_in)
    model.save(model_out)
    

    
