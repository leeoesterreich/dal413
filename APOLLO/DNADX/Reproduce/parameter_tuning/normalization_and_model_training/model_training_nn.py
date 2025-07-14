# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
import pickle
from scipy.stats import pearsonr
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F

class CNADataset(Dataset):
    """Custom Dataset for CNA data"""
    def __init__(self, X, y):
        # Convert to float32 for PyTorch
        X = np.asarray(X, dtype=np.float32)
        y = np.asarray(y, dtype=np.float32)
        
        self.X = torch.FloatTensor(X)
        self.y = torch.FloatTensor(y)
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

class DeepCNANet(nn.Module):
    """Deep Neural Network for CNA prediction"""
    def __init__(self, input_size, hidden_sizes=[256, 128, 64]):
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
        
        # Dropout for regularization
        self.dropout = nn.Dropout(0.2)
        
    def forward(self, x):
        for layer in self.layers:
            if isinstance(layer, nn.Linear):
                x = F.relu(layer(x))
                x = self.dropout(x)
            else:  # BatchNorm layer
                x = layer(x)
        return self.output(x)

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("training_data/rna_signature_score_median_no_norm.pkl")
    segment_score = pd.read_pickle("training_data/cna_segment_score_mean_no_norm.pkl")
    
    print("\nSample name formats:")
    print("RNA signature score sample names (first 5):", list(signature_score.columns[:5]))
    print("CNA segment score sample names (first 5):", list(segment_score.columns[:5]))
    
    # Convert CNA sample names to match RNA format (replace . with -)
    segment_score.columns = segment_score.columns.str.replace('.', '-')
    
    # Find common samples using exact names
    common_samples = sorted(set(signature_score.columns) & set(segment_score.columns))
    print("\nFound {} common samples".format(len(common_samples)))
    if len(common_samples) > 0:
        print("First 5 common samples:", common_samples[:5])
    
    # Align by common samples
    signature_score = signature_score[common_samples]
    segment_score = segment_score[common_samples]
    
    return signature_score, segment_score

def clean_data(X):
    """Clean data by handling NaN values"""
    if isinstance(X, pd.DataFrame):
        X = X.values
    
    # Convert to float type first
    X = X.astype(float)
    
    # Replace inf values with NaN
    X = np.where(np.isinf(X), np.nan, X)
    
    # Impute NaN values with column means
    imputer = SimpleImputer(strategy='mean')
    X_cleaned = imputer.fit_transform(X)
    
    return X_cleaned

def train_model(model, train_loader, val_loader, criterion, optimizer, device, num_epochs=100):
    """Train the neural network model without early stopping"""
    train_losses = []
    val_losses = []
    
    # Create models directory in results_nn
    os.makedirs('results_nn/models', exist_ok=True)
    best_model_path = os.path.join('results_nn/models', 'best_model.pth')
    
    best_val_loss = float('inf')
    
    for epoch in range(num_epochs):
        # Training phase
        model.train()
        train_loss = 0.0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y.view(-1, 1))
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation phase
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X)
                val_loss += criterion(outputs, batch_y.view(-1, 1)).item()
        
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        
        train_losses.append(train_loss)
        val_losses.append(val_loss)
        
        print('Epoch [{}/{}], Train Loss: {:.4f}, Val Loss: {:.4f}'.format(
            epoch+1, num_epochs, train_loss, val_loss))
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), best_model_path)
            print('    New best validation loss: {:.4f}'.format(val_loss))
    
    # Load best model
    model.load_state_dict(torch.load(best_model_path))
    return model, train_losses, val_losses

def calculate_tercile_roc(y_true, y_pred):
    """Calculate ROC metrics for top tercile vs bottom 2/3 classification"""
    from sklearn.metrics import roc_auc_score
    
    # Calculate tercile threshold (66.67th percentile)
    tercile_threshold = np.percentile(y_true, 66.67)
    
    # Create binary labels (1 for top tercile, 0 for bottom 2/3)
    y_true_binary = (y_true >= tercile_threshold).astype(int)
    
    # Calculate AUC
    auc_score = roc_auc_score(y_true_binary, y_pred)
    
    return auc_score

def build_prediction_model(X, y, signature_name, feature_names, batch_size=32, learning_rate=0.001):
    """Build and evaluate the neural network model"""
    print("\nBuilding neural network model for signature: {}".format(signature_name))
    
    # Ensure y is in the correct shape (n_samples, 1)
    y = y.reshape(-1, 1)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Clean and normalize data
    X_train = clean_data(X_train)
    X_test = clean_data(X_test)
    
    # Scale the data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Create datasets and dataloaders
    train_dataset = CNADataset(X_train_scaled, y_train)
    test_dataset = CNADataset(X_test_scaled, y_test)
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    # Initialize model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = DeepCNANet(input_size=X_train.shape[1]).to(device)
    
    # Define loss function and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    
    # Train the model
    model, train_losses, val_losses = train_model(
        model, train_loader, test_loader, criterion, optimizer, device
    )
    
    # Make predictions
    model.eval()
    with torch.no_grad():
        X_train_tensor = torch.FloatTensor(X_train_scaled).to(device)
        X_test_tensor = torch.FloatTensor(X_test_scaled).to(device)
        y_train_pred = model(X_train_tensor).cpu().numpy()
        y_test_pred = model(X_test_tensor).cpu().numpy()
    
    # Calculate metrics
    train_correlation, _ = pearsonr(y_train.ravel(), y_train_pred.ravel())
    test_correlation, _ = pearsonr(y_test.ravel(), y_test_pred.ravel())
    train_mse = mean_squared_error(y_train.ravel(), y_train_pred.ravel())
    test_mse = mean_squared_error(y_test.ravel(), y_test_pred.ravel())
    
    # Calculate tercile ROC AUC
    train_auc = calculate_tercile_roc(y_train.ravel(), y_train_pred.ravel())
    test_auc = calculate_tercile_roc(y_test.ravel(), y_test_pred.ravel())
    
    # Create results directories
    os.makedirs('results_nn/plots', exist_ok=True)
    os.makedirs('results_nn/models', exist_ok=True)
    
    # Plot performance
    plt.figure(figsize=(12, 8))
    
    # Plot training history
    plt.subplot(2, 2, 1)
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Validation Loss')
    plt.title('Training History')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    
    # Plot predictions vs actual
    plt.subplot(2, 2, 2)
    plt.scatter(y_test.ravel(), y_test_pred.ravel(), alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    plt.title('Test Set Predictions\nCorr: {:.3f}, MSE: {:.3f}'.format(test_correlation, test_mse))
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    
    # Plot tercile classification
    plt.subplot(2, 2, 3)
    tercile_threshold = np.percentile(y_test.ravel(), 66.67)
    y_test_binary = (y_test.ravel() >= tercile_threshold).astype(int)
    plt.hist(y_test_pred[y_test_binary == 0], bins=30, alpha=0.5, label='Bottom 2/3', density=True)
    plt.hist(y_test_pred[y_test_binary == 1], bins=30, alpha=0.5, label='Top 1/3', density=True)
    plt.title('Tercile Classification\nAUC: {:.3f}'.format(test_auc))
    plt.xlabel('Predicted Values')
    plt.ylabel('Density')
    plt.legend()
    
    # Save plot
    plt.tight_layout()
    plt.savefig('results_nn/plots/{}_performance.png'.format(signature_name))
    plt.close()
    
    # Save model and metadata
    model_info = {
        'model': model.state_dict(),
        'scaler': scaler,
        'feature_names': feature_names,
        'train_correlation': train_correlation,
        'test_correlation': test_correlation,
        'train_mse': train_mse,
        'test_mse': test_mse,
        'train_auc_tercile': train_auc,
        'test_auc_tercile': test_auc,
        'architecture': {
            'input_size': X_train.shape[1],
            'hidden_sizes': [256, 128, 64]
        }
    }
    
    torch.save(model_info, 'results_nn/models/{}_final_model.pt'.format(signature_name))
    
    return model_info

def monte_carlo_cv(X, y, signature_name, n_splits=200, train_size=0.7):
    """Perform Monte Carlo cross validation with neural network model"""
    print("\nProcessing signature: {}".format(signature_name))
    print("Starting Monte Carlo cross-validation with {} splits".format(n_splits))
    
    # Ensure y is in the correct shape (n_samples,)
    y = y.ravel()
    
    cv_scores = []
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    for i in range(n_splits):
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, train_size=train_size, random_state=i
        )
        
        # Clean and normalize data
        X_train = clean_data(X_train)
        X_test = clean_data(X_test)
        
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        # Create datasets and dataloaders
        train_dataset = CNADataset(X_train_scaled, y_train)
        test_dataset = CNADataset(X_test_scaled, y_test)
        
        train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=32)
        
        # Initialize model
        model = DeepCNANet(input_size=X_train.shape[1]).to(device)
        
        # Train model
        criterion = nn.MSELoss()
        optimizer = optim.Adam(model.parameters(), lr=0.001)
        
        model, _, _ = train_model(
            model, train_loader, test_loader, criterion, optimizer, device,
            num_epochs=50  # Reduced epochs for CV
        )
        
        # Evaluate
        model.eval()
        with torch.no_grad():
            X_test_tensor = torch.FloatTensor(X_test_scaled).to(device)
            y_test_pred = model(X_test_tensor).cpu().numpy()
        
        # Calculate metrics
        test_correlation, _ = pearsonr(y_test.ravel(), y_test_pred.ravel())
        test_mse = mean_squared_error(y_test.ravel(), y_test_pred.ravel())
        
        cv_scores.append({
            'split': i,
            'correlation': test_correlation,
            'mse': test_mse
        })
        
        if (i + 1) % 20 == 0:
            print("Completed {}/{} splits".format(i + 1, n_splits))
    
    return pd.DataFrame(cv_scores)

if __name__ == "__main__":
    print("Starting script execution...")
    try:
        # Create all necessary directories first
        print("Creating necessary directories...")
        os.makedirs('results_nn/models', exist_ok=True)
        os.makedirs('results_nn/plots', exist_ok=True)
        print("Directories created successfully.")
        
        # Load data
        print("Attempting to load data...")
        signature_score, segment_score = load_data()
        print("Data loaded successfully!")
        print("Signature score shape: {}".format(signature_score.shape))
        print("Segment score shape: {}".format(segment_score.shape))

        # Example: Train model for the first signature
        if not signature_score.empty:
            first_signature_name = signature_score.index[0]
            print("\nTraining model for signature: {}".format(first_signature_name))

            # Prepare data for the first signature
            X_data = segment_score.T  # Transpose so samples are rows, features are columns
            y_data = signature_score.loc[first_signature_name].values

            print("X_data shape: {}".format(X_data.shape))
            print("y_data shape: {}".format(y_data.shape))

            # Ensure X_data and y_data have the same number of samples
            if X_data.shape[0] != len(y_data):
                raise ValueError("Mismatch in number of samples between X ({}) and y ({})".format(
                    X_data.shape[0], len(y_data)))

            feature_names = list(X_data.columns)

            # Option 1: Build a single prediction model
            print("\nStarting model training...")
            model_info = build_prediction_model(X_data.values, y_data, first_signature_name, feature_names)
            print("\nModel training completed for {}.".format(first_signature_name))
            print("  Test Correlation: {:.4f}".format(model_info['test_correlation']))
            print("  Test MSE: {:.4f}".format(model_info['test_mse']))
        else:
            print("No signature scores found to train.")

    except Exception as e:
        print("An error occurred: {}".format(str(e)))
        import traceback
        traceback.print_exc()

    print("\nScript finished.") 