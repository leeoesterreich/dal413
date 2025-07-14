#!/usr/bin/env python3
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import numpy as np
import os
import argparse

# Define the DeepCNANet class (based on previous findings)
class DeepCNANet(nn.Module):
    """Deep Neural Network for CNA prediction"""
    def __init__(self, input_size, hidden_sizes=[256, 128, 64]): # Default hidden_sizes
        super(DeepCNANet, self).__init__()
        self.layers = nn.ModuleList()
        
        # Input layer
        self.layers.append(nn.Linear(input_size, hidden_sizes[0]))
        self.layers.append(nn.BatchNorm1d(hidden_sizes[0]))
        
        # Hidden layers
        for i in range(len(hidden_sizes)-1):
            self.layers.append(nn.Linear(hidden_sizes[i], hidden_sizes[i+1]))
            self.layers.append(nn.BatchNorm1d(hidden_sizes[i+1]))
        
        # Output layer
        self.output = nn.Linear(hidden_sizes[-1], 1)
        
        self.dropout = nn.Dropout(0.2) # Does not affect weight extraction directly
        
    def forward(self, x):
        for layer_idx, layer in enumerate(self.layers):
            if isinstance(layer, nn.Linear):
                x = layer(x)
                # Apply ReLU and Dropout only if it's not the last linear layer before batchnorm
                # and there's a subsequent batchnorm or it's not the output layer
                # This logic should match the original training script's forward pass for consistency
                # For weight extraction, only the structure matters, but good to be accurate.
                # The original forward pass in model_training_nn.py was:
                # if isinstance(layer, nn.Linear):
                #     x = F.relu(layer(x))
                #     x = self.dropout(x)
                # else:  # BatchNorm layer
                #     x = layer(x)
                # This means ReLU and dropout are applied *after* the linear transformation.
                # And batchnorm is applied to the result of that.
                
                # More accurate forward based on typical patterns:
                # Linear -> Activation -> Dropout (if applicable) -> BatchNorm (if applicable)
                # However, the provided structure from search was Linear -> BatchNorm for input and hidden.
                
                # Sticking to the simpler interpretation from previous `model_training_nn.py`
                # where activation and dropout follow linear, and batchnorm is separate.
                x = F.relu(x)
                x = self.dropout(x)
            elif isinstance(layer, nn.BatchNorm1d):
                x = layer(x)
        return self.output(x) # Output layer is applied last

def extract_feature_weights(model_path):
    """
    Loads a saved PyTorch model, extracts feature importances from the first linear layer,
    and saves them to a CSV file.
    """
    if not os.path.exists(model_path):
        print(f"Error: Model file not found at {model_path}")
        return

    print(f"Loading model from {model_path}...")
    try:
        # Try loading with map_location in case the model was saved on GPU and we are on CPU
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        print(f"Using device: {device}")
        model_info = torch.load(model_path, map_location=device)
    except Exception as e:
        print(f"Error loading model: {e}")
        return

    if not isinstance(model_info, dict):
        print("Error: Model file does not contain a dictionary (model_info).")
        return

    required_keys = ['architecture', 'model', 'feature_names']
    for key in required_keys:
        if key not in model_info:
            print(f"Error: Model info is missing required key: {key}")
            return
            
    architecture = model_info['architecture']
    state_dict = model_info['model']
    feature_names = model_info['feature_names']
    
    input_size = architecture.get('input_size')
    hidden_sizes = architecture.get('hidden_sizes', [256, 128, 64]) # Use default if not found

    if input_size is None:
        print("Error: 'input_size' not found in model architecture.")
        return

    if len(feature_names) != input_size:
        print(f"Warning: Number of feature names ({len(feature_names)}) does not match model input size ({input_size}).")
        # Decide if to proceed or exit. For now, proceed but this is a potential issue.

    print(f"Reconstructing model with input size: {input_size} and hidden sizes: {hidden_sizes}")
    model = DeepCNANet(input_size=input_size, hidden_sizes=hidden_sizes)
    model.load_state_dict(state_dict)
    model.to(device) # Move model to the device
    model.eval() # Set to evaluation mode

    # Extract weights from the first linear layer (self.layers[0])
    # self.layers[0] is nn.Linear(input_size, hidden_sizes[0])
    # Its weights have shape: (hidden_sizes[0], input_size)
    try:
        first_layer_weights = model.layers[0].weight.data
    except IndexError:
        print("Error: Could not access model.layers[0]. The model structure might be different than expected.")
        return
    except AttributeError:
        print("Error: model.layers[0] does not seem to be a linear layer with a 'weight' attribute.")
        return


    # Calculate importance: sum of absolute weights for each input feature across all neurons in the first hidden layer
    # This gives one importance value per input feature.
    # Weights: (num_hidden_neurons, num_input_features)
    # We want to sum along dim 0 (over hidden neurons) for each input feature (column)
    # Or, transpose and sum along dim 1.
    # W.abs().sum(dim=0) would be if features are columns.
    # If W is (out_features, in_features), then W.abs().sum(dim=0) gives per in_feature importance.
    importances_abs = first_layer_weights.abs().sum(dim=0).cpu().numpy()
    importances_net = first_layer_weights.sum(dim=0).cpu().numpy() # Net sum of weights

    if len(importances_abs) != len(feature_names):
        print(f"Error: Mismatch between number of absolute importances ({len(importances_abs)}) and feature names ({len(feature_names)}).")
        print(f"Weight matrix shape: {first_layer_weights.shape}")
        print(f"Expected input features (from feature_names): {len(feature_names)}")
        print(f"Expected input features (from input_size): {input_size}")
        return
    
    if len(importances_net) != len(feature_names):
        print(f"Error: Mismatch between number of net importances ({len(importances_net)}) and feature names ({len(feature_names)}).")
        return

    # Create a DataFrame
    df_weights = pd.DataFrame({
        'Segment': feature_names,
        'Importance_Abs': importances_abs,
        'Importance_Net': importances_net
    })

    # Sort by absolute importance
    df_weights = df_weights.sort_values(by='Importance_Abs', ascending=False)

    # Save to CSV
    output_dir = os.path.dirname(model_path)
    model_filename = os.path.basename(model_path)
    
    # Construct the output filename
    if model_filename.endswith('_final_model.pt'):
        base_name = model_filename[:-len('_final_model.pt')]
    else: # Fallback if the suffix is different
        base_name = os.path.splitext(model_filename)[0]
        
    output_csv_filename = f"{base_name}_feature_weights.csv"
    output_csv_path = os.path.join(output_dir, output_csv_filename)

    try:
        df_weights.to_csv(output_csv_path, index=False)
        print(f"Successfully saved feature weights to: {output_csv_path}")
    except Exception as e:
        print(f"Error saving CSV file: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract feature weights from a DeepCNANet PyTorch model.")
    parser.add_argument("model_path", type=str, help="Path to the saved .pt model file.")
    
    args = parser.parse_args()
    extract_feature_weights(args.model_path) 