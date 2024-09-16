#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Title: Machine Learning Model Evaluation Script

Author: Erika Castaneda

Version: Python 3.11.5 

Description:
    This script trains and evaluates multiple machine learning models on a dataset 
    to compare their performance. It handles class imbalance using SMOTE, scales features 
    using StandardScaler, and performs hyperparameter tuning with GridSearchCV. The models 
    evaluated include Logistic Regression, Support Vector Machine (SVM), Random Forest,
    and K-Nearest Neighbors (KNN). The script generates and displays evaluation metrics, 
    confusion matrices, ROC curves, and summary results for each model.

Inputs:
    - `data.csv`: A CSV file containing the dataset with features and target variable.

Outputs:
    - Confusion matrices and ROC curves for each model, displayed as plots.
    - Classification reports for each model, printed to the console.
    - A summary DataFrame containing accuracy, AUC, precision, recall, and F1 
      scores for each model and feature set, printed to the console.

Usage:
    1. Ensure all required libraries are installed.
    2. Replace 'data.csv' with the path to your dataset.
    3. Run the script in a Python environment that supports plotting (e.g., JupyterLab).
    4. Review the output plots and printed metrics to evaluate model performance.

Notes:
    - Adjust the feature sets, hyperparameters, and file paths as needed based on your specific dataset and analysis requirements.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, precision_score, recall_score, f1_score, accuracy_score, roc_curve, auc, confusion_matrix
from imblearn.over_sampling import SMOTE
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

# Function to train and evaluate a model
def train_and_evaluate_model(model, X_train, y_train, X_test, y_test, model_name):
    """
    Trains and evaluates a machine learning model.
    
    Parameters:
    - model: The machine learning model to be trained.
    - X_train: Features for training.
    - y_train: Target variable for training.
    - X_test: Features for testing.
    - y_test: Target variable for testing.
    - model_name: Name of the model for display purposes.
    
    Returns:
    - DataFrame containing the evaluation metrics.
    """
    # Apply SMOTE to handle class imbalance
    smote = SMOTE(sampling_strategy='minority', random_state=42)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

    # Train the model
    model.fit(X_train_resampled, y_train_resampled)
    y_pred = model.predict(X_test)

    # Evaluate model performance
    accuracy = accuracy_score(y_test, y_pred)

    # Show and save confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(6, 4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False)
    plt.title(f'{model_name} Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.show()  # Display in JupyterLab

    # Generate and save ROC curve
    fpr, tpr, _ = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
    roc_auc = auc(fpr, tpr)
    plt.figure(figsize=(6, 4))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'{model_name} (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'{model_name} ROC Curve')
    plt.legend(loc='lower right')
    plt.show()

    # Display classification report
    class_report = classification_report(y_test, y_pred, target_names=['Class 0', 'Class 1'])
    print(f'{model_name} Classification Report:\n{class_report}')

    # Calculate macro-averaged precision, recall, and F1-score
    macro_precision = precision_score(y_test, y_pred, average='macro')
    macro_recall = recall_score(y_test, y_pred, average='macro')
    macro_f1 = f1_score(y_test, y_pred, average='macro')

    # Create a DataFrame with the results
    results_df = pd.DataFrame({
        'Model': [model_name],
        'Accuracy': [accuracy],
        'AUC': [roc_auc],
        'Precision': [macro_precision],
        'Recall': [macro_recall],
        'F1': [macro_f1]
    })

    return results_df

# Load and prepare data (replace 'xy_data' with actual DataFrame)
xy_data = pd.read_csv('your_data.csv')  # Adjust as needed

# Experiment with different feature sets and models
features_list = [
    ['coverage_dif'],
    ['Heterozigosity_dif'],
    ['GC'],
    ['GC_dif'],
    ['coverage_dif', 'GC'],
    ['coverage_dif', 'Heterozigosity_dif'],
    ['GC', 'Heterozigosity_dif'],
    ['GC', 'Heterozigosity_dif', 'coverage_dif'],
    ['GC_dif', 'Heterozigosity_dif', 'coverage_dif']
]

results = []

# Iterate over each feature set
for features in features_list:
    X = xy_data[features]
    y = xy_data['class']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Apply SMOTE
    smote = SMOTE(random_state=42)
    X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

    # Apply StandardScaler
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_resampled)
    X_test_scaled = scaler.transform(X_test)

    # Train and evaluate Logistic Regression
    param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100]}
    grid_search = GridSearchCV(LogisticRegression(), param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X_train_scaled, y_train)
    best_C = grid_search.best_params_['C']
    logreg_model = LogisticRegression(class_weight='balanced', C=best_C)
    results.append(train_and_evaluate_model(logreg_model, X_train_scaled, y_train, X_test_scaled, y_test, 'Logistic Regression'))

    # Train and evaluate SVM
    param_grid = {'C': [0.1, 1, 10], 'gamma': [0.01, 0.1, 1]}
    svm_model = SVC(probability=True)
    grid_search = GridSearchCV(svm_model, param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X_train_scaled, y_train)
    best_params = grid_search.best_params_
    svm_model = SVC(probability=True, kernel='rbf', class_weight='balanced', **best_params)
    results.append(train_and_evaluate_model(svm_model, X_train_scaled, y_train, X_test_scaled, y_test, 'Support Vector Machine'))

    # Train and evaluate Random Forest
    param_grid = {
        'n_estimators': [50, 100, 200],
        'max_depth': [None, 2, 10, 20],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4]
    }
    rf_model = RandomForestClassifier(random_state=42)
    grid_search = GridSearchCV(rf_model, param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X_train_scaled, y_train)
    best_params = grid_search.best_params_
    rf_model = RandomForestClassifier(class_weight='balanced', **best_params)
    results.append(train_and_evaluate_model(rf_model, X_train_scaled, y_train, X_test_scaled, y_test, 'Random Forest'))

    # Train and evaluate K-Nearest Neighbors
    param_grid = {'n_neighbors': [3, 5, 7, 9], 'p': [1, 2]}
    knn_model = KNeighborsClassifier()
    grid_search = GridSearchCV(knn_model, param_grid, cv=5, scoring='accuracy')
    grid_search.fit(X_train_scaled, y_train)
    best_params = grid_search.best_params_
    knn_model = KNeighborsClassifier(**best_params)
    results.append(train_and_evaluate_model(knn_model, X_train_scaled, y_train, X_test_scaled, y_test, 'K-Nearest Neighbors'))

# Display combined ROC curve for final comparison
X_final = xy_data[['coverage_dif', 'Heterozigosity_dif', 'GC']]
y_final = xy_data['class']
X_train_final, X_test_final, y_train_final, y_test_final = train_test_split(X_final, y_final, test_size=0.2, random_state=42)
scaler_final = StandardScaler()
X_train_scaled_final = scaler_final.fit_transform(X_train_final)
X_test_scaled_final = scaler_final.transform(X_test_final)
smote_final = SMOTE(sampling_strategy='minority', random_state=42)
X_train_resampled_final, y_train_resampled_final = smote_final.fit_resample(X_train_scaled_final, y_train_final)

# Train models on final feature set
rf_model = RandomForestClassifier(class_weight='balanced', n_estimators=100, max_depth=None, min_samples_split=2, min_samples_leaf=4)
svm_model = SVC(probability=True, kernel='rbf', class_weight='balanced', C=0.1, gamma=0.01)
knn_model = KNeighborsClassifier(n_neighbors=5, p=1)
logreg_model = LogisticRegression(class_weight='balanced', C=1)

# Fit models
rf_model.fit(X_train_resampled_final, y_train_resampled_final)
svm_model.fit(X_train_resampled_final, y_train_resampled_final)
knn_model.fit(X_train_resampled_final, y_train_resampled_final)
logreg_model.fit(X_train_resampled_final, y_train_resampled_final)

# Predict probabilities
y_score_rf = rf_model.predict_proba(X_test_scaled_final)[:, 1]
y_score_svm = svm_model.predict_proba(X_test_scaled_final)[:, 1]
y_score_knn = knn_model.predict_proba(X_test_scaled_final)[:, 1]
y_score_lr = logreg_model.predict_proba(X_test_scaled_final)[:, 1]

# ROC Curve for final comparison
fpr_rf, tpr_rf, _ = roc_curve(y_test_final, y_score_rf)
roc_auc_rf = auc(fpr_rf, tpr_rf)
fpr_svm, tpr_svm, _ = roc_curve(y_test_final, y_score_svm)
roc_auc_svm = auc(fpr_svm, tpr_svm)
fpr_knn, tpr_knn, _ = roc_curve(y_test_final, y_score_knn)
roc_auc_knn = auc(fpr_knn, tpr_knn)
fpr_lr, tpr_lr, _ = roc_curve(y_test_final, y_score_lr)
roc_auc_lr = auc(fpr_lr, tpr_lr)

plt.figure(figsize=(10, 6))
plt.plot(fpr_rf, tpr_rf, color='blue', lw=2, label=f'Random Forest (AUC = {roc_auc_rf:.2f})')
plt.plot(fpr_svm, tpr_svm, color='green', lw=2, label=f'SVM (AUC = {roc_auc_svm:.2f})')
plt.plot(fpr_knn, tpr_knn, color='orange', lw=2, label=f'KNN (AUC = {roc_auc_knn:.2f})')
plt.plot(fpr_lr, tpr_lr, color='red', lw=2, label=f'Logistic Regression (AUC = {roc_auc_lr:.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve Comparison')
plt.legend(loc="lower right")
plt.show()

# Combine results from all experiments
results_df = pd.concat(results, ignore_index=True)
print(results_df)
