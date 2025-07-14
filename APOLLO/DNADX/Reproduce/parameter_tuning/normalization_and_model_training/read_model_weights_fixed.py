# -*- coding: utf-8 -*-
import torch
import numpy as np
import pandas as pd
import os

def read_model_weights(model_path):
    """Read and analyze the model weights from a .pth file"""
    print("Reading model weights from: {}".format(model_path))
    
    if not os.path.exists(model_path):
        print("Error: Model file not found: {}".format(model_path))
        return None, None
    
    # Load the model weights
    try:
        weights = torch.load(model_path, map_location='cpu')
        print("Model weights loaded successfully!")
        print("Type of loaded object: {}".format(type(weights)))
        
        if isinstance(weights, dict):
            print("\nModel state dictionary contains {} layers:".format(len(weights)))
            print("-" * 80)
            
            layer_info = []
            total_params = 0
            
            for key, tensor in weights.items():
                shape = tuple(tensor.shape)
                num_params = tensor.numel()
                total_params += num_params
                
                # Handle different tensor types
                tensor_info = {
                    'Layer': key,
                    'Shape': shape,
                    'Parameters': num_params,
                    'Dtype': str(tensor.dtype),
                }
                
                # Only calculate statistics for floating point tensors
                if tensor.dtype.is_floating_point:
                    tensor_info.update({
                        'Min': float(tensor.min()),
                        'Max': float(tensor.max()),
                        'Mean': float(tensor.mean()),
                        'Std': float(tensor.std())
                    })
                else:
                    tensor_info.update({
                        'Min': 'N/A (non-float)',
                        'Max': 'N/A (non-float)',
                        'Mean': 'N/A (non-float)',
                        'Std': 'N/A (non-float)'
                    })
                
                layer_info.append(tensor_info)
                
                print("Layer: {:<35} Shape: {:<20} Params: {:>8} Type: {}".format(
                    key, str(shape), num_params, str(tensor.dtype)))
            
            print("-" * 80)
            print("Total parameters: {:,}".format(total_params))
            
            # Create detailed analysis
            df = pd.DataFrame(layer_info)
            print("\nDetailed Weight Statistics:")
            print(df.to_string(index=False))
            
            # Analyze weight distributions for floating point weights only
            print("\nWeight Distribution Analysis (Float tensors only):")
            print("=" * 80)
            
            for key, tensor in weights.items():
                if tensor.dtype.is_floating_point and 'weight' in key:
                    print("\n{}:".format(key))
                    print("  Shape: {}".format(tuple(tensor.shape)))
                    print("  Data type: {}".format(tensor.dtype))
                    print("  Min: {:.6f}, Max: {:.6f}".format(float(tensor.min()), float(tensor.max())))
                    print("  Mean: {:.6f}, Std: {:.6f}".format(float(tensor.mean()), float(tensor.std())))
                    
                    # Count zeros and near-zeros
                    zeros = (tensor == 0).sum()
                    near_zeros = (tensor.abs() < 1e-6).sum()
                    print("  Exact zeros: {} ({:.2f}%)".format(
                        int(zeros), float(zeros) / tensor.numel() * 100))
                    print("  Near-zeros (|x| < 1e-6): {} ({:.2f}%)".format(
                        int(near_zeros), float(near_zeros) / tensor.numel() * 100))
                    
                    # Check for dead neurons (all zeros) in linear layers
                    if len(tensor.shape) == 2:  # Linear layer weights
                        dead_input_neurons = (tensor.abs().sum(dim=0) == 0).sum()
                        dead_output_neurons = (tensor.abs().sum(dim=1) == 0).sum()
                        print("  Dead input neurons: {} out of {}".format(
                            int(dead_input_neurons), tensor.shape[1]))
                        print("  Dead output neurons: {} out of {}".format(
                            int(dead_output_neurons), tensor.shape[0]))
                        
                        # Weight magnitude analysis
                        weight_norms = tensor.norm(dim=1)  # L2 norm of each output neuron
                        print("  Output neuron weight norms - Min: {:.6f}, Max: {:.6f}, Mean: {:.6f}".format(
                            float(weight_norms.min()), float(weight_norms.max()), float(weight_norms.mean())))
            
            # Analyze batch norm parameters
            print("\nBatch Normalization Parameters:")
            print("=" * 60)
            for key, tensor in weights.items():
                if 'running_mean' in key or 'running_var' in key:
                    if tensor.dtype.is_floating_point:
                        print("\n{}:".format(key))
                        print("  Shape: {}".format(tuple(tensor.shape)))
                        print("  Min: {:.6f}, Max: {:.6f}".format(float(tensor.min()), float(tensor.max())))
                        print("  Mean: {:.6f}, Std: {:.6f}".format(float(tensor.mean()), float(tensor.std())))
                    else:
                        print("\n{}: (Non-float type: {})".format(key, tensor.dtype))
            
            # Show first few values of each layer for inspection
            print("\nSample Weight Values (first 5 elements):")
            print("=" * 60)
            for key, tensor in weights.items():
                if tensor.dtype.is_floating_point:
                    flat_tensor = tensor.flatten()
                    sample_size = min(5, len(flat_tensor))
                    sample_values = [float(flat_tensor[i]) for i in range(sample_size)]
                    print("{:<35}: {}".format(key, sample_values))
            
            # Save weights to CSV for detailed analysis
            print("\nSaving weight statistics to CSV...")
            df.to_csv('model_weights_analysis.csv', index=False)
            print("Weight analysis saved to: model_weights_analysis.csv")
            
            return weights, df
            
        else:
            print("Warning: Expected dictionary, got {}".format(type(weights)))
            return weights, None
            
    except Exception as e:
        print("Error loading model weights: {}".format(str(e)))
        import traceback
        traceback.print_exc()
        return None, None

def extract_specific_layer_weights(weights, layer_name):
    """Extract and display specific layer weights"""
    if layer_name in weights:
        tensor = weights[layer_name]
        print("\nDetailed analysis of layer: {}".format(layer_name))
        print("-" * 50)
        print("Shape: {}".format(tuple(tensor.shape)))
        print("Data type: {}".format(tensor.dtype))
        
        if tensor.dtype.is_floating_point:
            print("Min: {:.6f}".format(float(tensor.min())))
            print("Max: {:.6f}".format(float(tensor.max())))
            print("Mean: {:.6f}".format(float(tensor.mean())))
            print("Std: {:.6f}".format(float(tensor.std())))
            
            # Show actual weight matrix (if small enough)
            if tensor.numel() <= 100:  # Only show if small
                print("\nActual weights:")
                if len(tensor.shape) == 1:
                    print(tensor.numpy())
                else:
                    print(tensor.numpy())
            else:
                print("\nWeight matrix too large to display ({} elements)".format(tensor.numel()))
                print("First 10 elements: {}".format(tensor.flatten()[:10].tolist()))
        else:
            print("Non-floating point tensor, showing raw values:")
            if tensor.numel() <= 20:
                print(tensor.tolist())
            else:
                print("First 10 elements: {}".format(tensor.flatten()[:10].tolist()))
        
        return tensor
    else:
        print("Layer '{}' not found in model weights".format(layer_name))
        print("Available layers: {}".format(list(weights.keys())))
        return None

if __name__ == "__main__":
    print("Model Weight Analysis Tool")
    print("=" * 80)
    
    # Analyze the results_nn model (most recent one)
    model_path = './results_nn/models/best_model.pth'
    
    if os.path.exists(model_path):
        print("\nAnalyzing model file: {}".format(model_path))
        weights, df = read_model_weights(model_path)
        
        if weights is not None:
            print("\n" + "="*80)
            print("SPECIFIC LAYER ANALYSIS")
            print("="*80)
            
            # Analyze first layer weights in detail
            first_layer = extract_specific_layer_weights(weights, 'layers.0.weight')
            
            print("\n" + "="*80)
            print("ARCHITECTURE SUMMARY")
            print("="*80)
            
            # Reconstruct the network architecture from weights
            print("Neural Network Architecture:")
            linear_layers = [key for key in weights.keys() if 'weight' in key and 'layers.' in key and 'running' not in key]
            for i, layer_key in enumerate(linear_layers):
                weight_shape = weights[layer_key].shape
                if len(weight_shape) == 2:
                    print("  Layer {}: {} -> {} neurons".format(i+1, weight_shape[1], weight_shape[0]))
            
            # Output layer
            if 'output.weight' in weights:
                output_shape = weights['output.weight'].shape
                print("  Output: {} -> {} (final output)".format(output_shape[1], output_shape[0]))
        
    else:
        print("Model file not found: {}".format(model_path))
    
    print("\nAnalysis completed!") 