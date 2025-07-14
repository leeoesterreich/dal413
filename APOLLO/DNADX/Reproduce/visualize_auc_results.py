import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet
from sklearn.metrics import roc_auc_score, mean_squared_error, roc_curve, auc
import os
from datetime import datetime

def load_data():
    """Load the signature scores and segment scores, align by exact sample names"""
    print("Loading data...")
    signature_score = pd.read_pickle("results/signature_score.pkl")
    segment_score = pd.read_pickle("results/segment_score.pkl")
    
    # Find common samples using exact names
    common_samples = set(signature_score.columns) & set(segment_score.columns)
    print(f"\nFound {len(common_samples)} common samples")
    
    # Align by common samples
    signature_score = signature_score[list(common_samples)]
    segment_score = segment_score[list(common_samples)]
    
    return signature_score, segment_score

def train_test_analysis(signature_score, segment_score, signature_name):
    """Build prediction model for a target signature using segment scores"""
    print(f"\nBuilding prediction model for {signature_name}...")
    
    # Prepare data
    X = segment_score.T.astype(float)
    y = signature_score.loc[signature_name].astype(float)
    
    # Remove columns in X that are all NaN
    X = X.dropna(axis=1, how='all')
    # Fill remaining NaNs with column medians
    X = X.fillna(X.median())
    
    # Keep only samples where y has finite values
    mask = np.isfinite(y)
    X = X[mask]
    y = y[mask]
    
    if len(X) < 10:  # Require at least 10 samples
        print(f"Not enough valid samples after filtering (need at least 10, got {len(X)})")
        return None
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # Train model
    model = ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42)
    model.fit(X_train, y_train)
    
    # Get predictions
    train_pred = model.predict(X_train)
    test_pred = model.predict(X_test)
    
    # Calculate metrics
    train_corr = np.corrcoef(y_train, train_pred)[0, 1]
    test_corr = np.corrcoef(y_test, test_pred)[0, 1]
    train_mse = mean_squared_error(y_train, train_pred)
    test_mse = mean_squared_error(y_test, test_pred)
    
    # For ROC/AUC, binarize using median as threshold
    y_train_bin = (y_train > np.median(y_train)).astype(int)
    y_test_bin = (y_test > np.median(y_test)).astype(int)
    
    fpr_train, tpr_train, _ = roc_curve(y_train_bin, train_pred)
    fpr_test, tpr_test, _ = roc_curve(y_test_bin, test_pred)
    train_auc = auc(fpr_train, tpr_train)
    test_auc = auc(fpr_test, tpr_test)
    
    return {
        'signature': signature_name,
        'train_auc': train_auc,
        'test_auc': test_auc,
        'train_corr': train_corr,
        'test_corr': test_corr,
        'train_mse': train_mse,
        'test_mse': test_mse
    }

def load_metrics():
    """Load model metrics from CSV file"""
    metrics_file = '/ix1/alee/LO_LAB/Personal/Daisong/Projects/APOLLO/DNADX/Reproduce/results/plots/model_metrics.csv'
    if not os.path.exists(metrics_file):
        raise FileNotFoundError(f"Metrics file not found: {metrics_file}")
    
    metrics_df = pd.read_csv(metrics_file)
    print(f"Loaded metrics for {len(metrics_df)} signatures")
    return metrics_df

def create_visualizations(metrics_df):
    """Create various visualizations of the model metrics"""
    # Create timestamp for unique filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create plots directory if it doesn't exist
    os.makedirs('results/plots', exist_ok=True)
    
    # 1. Create a figure with multiple subplots
    plt.figure(figsize=(20, 15))
    
    # 1.1 Distribution of metrics
    plt.subplot(2, 2, 1)
    sns.boxplot(data=metrics_df[['train_correlation', 'test_correlation']])
    sns.stripplot(data=metrics_df[['train_correlation', 'test_correlation']], 
                 color='black', alpha=0.3, size=4)
    plt.title('Distribution of Train/Test Correlations')
    plt.ylabel('Correlation')
    
    # 1.2 Top 20 signatures by test correlation
    top_20 = metrics_df.sort_values('test_correlation', ascending=False).head(20)
    plt.subplot(2, 2, 2)
    plt.scatter(range(len(top_20)), top_20['test_correlation'], alpha=0.6)
    plt.plot(range(len(top_20)), top_20['test_correlation'], 'r--', alpha=0.3)
    plt.xticks(range(len(top_20)), [s.split('_')[0] for s in top_20['signature']], rotation=45, ha='right')
    plt.title('Top 20 Signatures by Test Correlation')
    plt.ylabel('Correlation')
    
    # 1.3 Train vs Test correlation comparison
    plt.subplot(2, 2, 3)
    plt.scatter(metrics_df['train_correlation'], metrics_df['test_correlation'], alpha=0.5)
    plt.plot([0, 1], [0, 1], 'r--', alpha=0.3)  # Diagonal line
    plt.xlabel('Train Correlation')
    plt.ylabel('Test Correlation')
    plt.title('Train vs Test Correlation')
    
    # 1.4 Distribution of MSE
    plt.subplot(2, 2, 4)
    sns.boxplot(data=metrics_df[['train_mse', 'test_mse']])
    sns.stripplot(data=metrics_df[['train_mse', 'test_mse']], 
                 color='black', alpha=0.3, size=4)
    plt.title('Distribution of Train/Test MSE')
    plt.ylabel('MSE')
    
    plt.tight_layout()
    plt.savefig(f'results/plots/metrics_distribution_{timestamp}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Create AUC-specific visualizations
    plt.figure(figsize=(20, 10))
    
    # 2.1 Distribution of AUC scores
    plt.subplot(1, 2, 1)
    sns.boxplot(data=metrics_df[['train_auc', 'test_auc']])
    sns.stripplot(data=metrics_df[['train_auc', 'test_auc']], 
                 color='black', alpha=0.3, size=4)
    plt.title('Distribution of Train/Test AUC Scores')
    plt.ylabel('AUC Score')
    
    # 2.2 Top 20 signatures by AUC
    top_20_auc = metrics_df.sort_values('test_auc', ascending=False).head(20)
    plt.subplot(1, 2, 2)
    plt.scatter(range(len(top_20_auc)), top_20_auc['test_auc'], label='Test', alpha=0.6)
    plt.scatter(range(len(top_20_auc)), top_20_auc['train_auc'], label='Train', alpha=0.6)
    plt.xticks(range(len(top_20_auc)), [s.split('_')[0] for s in top_20_auc['signature']], rotation=45, ha='right')
    plt.title('Top 20 Signatures by AUC Score')
    plt.ylabel('AUC Score')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'results/plots/auc_distribution_{timestamp}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Create a summary table of top performers
    top_10_corr = metrics_df.sort_values('test_correlation', ascending=False).head(10)
    top_10_auc = metrics_df.sort_values('test_auc', ascending=False).head(10)
    
    # Save top performers to CSV with timestamp
    top_10_corr.to_csv(f'results/plots/top_10_correlation_{timestamp}.csv', index=False)
    top_10_auc.to_csv(f'results/plots/top_10_auc_{timestamp}.csv', index=False)
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print("=" * 80)
    print("\nCorrelation Metrics:")
    print(f"Mean Train Correlation: {metrics_df['train_correlation'].mean():.3f}")
    print(f"Mean Test Correlation: {metrics_df['test_correlation'].mean():.3f}")
    print(f"Max Test Correlation: {metrics_df['test_correlation'].max():.3f}")
    
    print("\nAUC Metrics:")
    print(f"Mean Train AUC: {metrics_df['train_auc'].mean():.3f}")
    print(f"Mean Test AUC: {metrics_df['test_auc'].mean():.3f}")
    print(f"Max Test AUC: {metrics_df['test_auc'].max():.3f}")
    
    print("\nMSE Metrics:")
    print(f"Mean Train MSE: {metrics_df['train_mse'].mean():.3f}")
    print(f"Mean Test MSE: {metrics_df['test_mse'].mean():.3f}")
    print(f"Min Test MSE: {metrics_df['test_mse'].min():.3f}")
    
    print("\nTop 5 Signatures by Test Correlation:")
    print("=" * 80)
    for i, (_, row) in enumerate(top_10_corr.head(5).iterrows(), 1):
        print(f"{i}. {row['signature']}")
        print(f"   Test Correlation: {row['test_correlation']:.3f}")
        print(f"   Test AUC: {row['test_auc']:.3f}")
        print(f"   Test MSE: {row['test_mse']:.3f}")
        print("-" * 80)

def main():
    try:
        # Load metrics
        metrics_df = load_metrics()
        
        # Create visualizations
        create_visualizations(metrics_df)
        
        print("\nVisualization complete! Check the results/plots directory for output files.")
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 