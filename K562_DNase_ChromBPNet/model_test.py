import tensorflow as tf
import bpnet.model.arch
import time, datetime, json, os, sys
import pandas as pd
import numpy as np
import pysam
import random
from scipy.special import logsumexp

def inspect_model_signature(model_path):
    model = tf.saved_model.load(model_path)
    print(f"Model type: {type(model)}")
    
    # Check if it has signatures
    if hasattr(model, 'signatures'):
        print("Available signatures:", list(model.signatures.keys()) if hasattr(model.signatures, 'keys') else model.signatures)
        
        # Try to access serving_default
        if 'serving_default' in model.signatures:
            sig = model.signatures['serving_default']
            print("serving_default signature:")
            
            # Handle inputs as list
            print("Inputs:")
            if hasattr(sig.inputs, 'items'):
                for name, tensor in sig.inputs.items():
                    print(f"  {name}: shape {tensor.shape}, dtype {tensor.dtype}")
            else:
                # It's a list
                for i, tensor in enumerate(sig.inputs):
                    print(f"  input[{i}]: shape {tensor.shape}, dtype {tensor.dtype}")
            
            # Handle outputs as list
            print("Outputs:")
            if hasattr(sig.outputs, 'items'):
                for name, tensor in sig.outputs.items():
                    print(f"  {name}: shape {tensor.shape}, dtype {tensor.dtype}")
            else:
                # It's a list
                for i, tensor in enumerate(sig.outputs):
                    print(f"  output[{i}]: shape {tensor.shape}, dtype {tensor.dtype}")
    
    # Test with sample data
    test_seq = np.random.rand(2, 2114, 4).astype('float32')
    
    # Try calling via signature
    if hasattr(model, 'signatures') and 'serving_default' in model.signatures:
        try:
            # Try different ways to call the signature
            sig = model.signatures['serving_default']
            
            # Method 1: positional argument
            try:
                result = sig(test_seq)
                print(f"Signature call (positional) result type: {type(result)}")
                if isinstance(result, (list, tuple)):
                    print(f"Result has {len(result)} elements:")
                    for i, item in enumerate(result):
                        if hasattr(item, 'shape'):
                            print(f"  [{i}]: shape {item.shape}, dtype {item.dtype}")
                        else:
                            print(f"  [{i}]: {type(item)}")
                elif hasattr(result, 'shape'):
                    print(f"Result shape: {result.shape}")
                else:
                    print(f"Result: {result}")
            except Exception as e:
                print(f"Positional call failed: {e}")
                
                # Method 2: Try with keyword argument 'sequence'
                try:
                    result = sig(sequence=test_seq)
                    print(f"Signature call (sequence=) successful")
                except Exception as e2:
                    print(f"Keyword call also failed: {e2}")
                    
        except Exception as e:
            print(f"Signature access failed: {e}")
    
    # Try direct call
    try:
        result = model(test_seq)
        print(f"Direct call result type: {type(result)}")
        if isinstance(result, (list, tuple)):
            print(f"Direct call - List with {len(result)} elements:")
            for i, item in enumerate(result):
                if hasattr(item, 'shape'):
                    print(f"  [{i}]: shape {item.shape}, dtype {item.dtype}")
                else:
                    print(f"  [{i}]: {type(item)} - {str(item)[:100]}")
    except Exception as e:
        print(f"Direct call failed: {e}")
    
    return model