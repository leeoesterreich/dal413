import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet, ElasticNetCV, MultiTaskElasticNetCV
from sklearn.metrics import roc_curve, auc, accuracy_score, roc_auc_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import os
from joblib import Parallel, delayed, dump, load
from sklearn.impute import SimpleImputer
import warnings
import pickle
from sklearn.linear_model import LogisticRegressionCV
from scipy.stats import pearsonr

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

def clean_data(X, y):
    """Clean data by removing NaN values and standardizing"""
    # Remove samples with NaN in y
    mask = ~np.isnan(y)
    X = X[mask]
    y = y[mask]
    
    # Standardize features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    
    return X, y

def train_model(X, y, alpha=1.0, l1_ratio=0.5):
    """Train elastic net model"""
    model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, random_state=42)
    model.fit(X, y)
    return model

def evaluate_model(model, X, y):
    """Evaluate model performance"""
    y_pred = model.predict(X)
    mse = mean_squared_error(y, y_pred)
    correlation = np.corrcoef(y, y_pred)[0, 1]
    
    # Calculate binary classification metrics
    threshold = np.percentile(y, 66.67)  # Top 1/3 threshold
    y_binary = (y > threshold).astype(int)
    y_pred_binary = (y_pred > threshold).astype(int)
    
    roc_auc = roc_auc_score(y_binary, y_pred)
    accuracy = accuracy_score(y_binary, y_pred_binary)
    
    return {
        'mse': mse,
        'correlation': correlation,
        'roc_auc': roc_auc,
        'accuracy': accuracy
    }

def plot_model_performance(y_true, y_pred, title="Model Performance"):
    """Plot actual vs predicted values and ROC curve"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Actual vs Predicted plot
    ax1.scatter(y_true, y_pred, alpha=0.5)
    ax1.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2)
    ax1.set_xlabel('Actual Values')
    ax1.set_ylabel('Predicted Values')
    ax1.set_title('Actual vs Predicted')
    
    # ROC curve
    threshold = np.percentile(y_true, 66.67)
    y_binary = (y_true > threshold).astype(int)
    fpr, tpr, _ = roc_curve(y_binary, y_pred)
    roc_auc = auc(fpr, tpr)
    
    ax2.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    ax2.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax2.set_xlabel('False Positive Rate')
    ax2.set_ylabel('True Positive Rate')
    ax2.set_title('ROC Curve')
    ax2.legend(loc="lower right")
    
    plt.tight_layout()
    return fig

def main():
    # Load data
    signature_score, segment_score = load_data()
    
    # Create results directory
    os.makedirs('results/models', exist_ok=True)
    os.makedirs('results/plots', exist_ok=True)
    
    # Store all metrics
    all_metrics = []
    
    # Train model for each signature
    for signature in signature_score.index:
        print(f"\nTraining model for signature: {signature}")
        
        # Prepare data
        X = segment_score.values.T
        y = signature_score.loc[signature].values
        
        # Clean data
        X_clean, y_clean = clean_data(X, y)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X_clean, y_clean, test_size=0.2, random_state=42)
        
        # Train model
        model = train_model(X_train, y_train)
        
        # Evaluate model
        train_metrics = evaluate_model(model, X_train, y_train)
        test_metrics = evaluate_model(model, X_test, y_test)
        
        # Save model and feature names
        model_info = {
            'model': model,
            'feature_names': list(segment_score.index)
        }
        safe_name = signature.replace('/', '_').replace('\\', '_')
        dump(model_info, f'results/models/{safe_name}_model.pkl')
        
        # Plot and save performance plots
        fig = plot_model_performance(y_test, model.predict(X_test),
                                   title=f"Model Performance - {signature}")
        fig.savefig(f'results/plots/{safe_name}_performance.png')
        plt.close(fig)
        
        # Store metrics
        metrics = {
            'signature': signature,
            'train_mse': train_metrics['mse'],
            'test_mse': test_metrics['mse'],
            'train_correlation': train_metrics['correlation'],
            'test_correlation': test_metrics['correlation'],
            'train_roc_auc': train_metrics['roc_auc'],
            'test_roc_auc': test_metrics['roc_auc'],
            'train_accuracy': train_metrics['accuracy'],
            'test_accuracy': test_metrics['accuracy']
        }
        all_metrics.append(metrics)
        
        print("Training metrics:", train_metrics)
        print("Testing metrics:", test_metrics)
    
    # Save all metrics
    metrics_df = pd.DataFrame(all_metrics)
    metrics_df.to_csv('results/all_model_metrics.csv', index=False)
    print("\nAll metrics saved to results/all_model_metrics.csv")

if __name__ == "__main__":
    main() 