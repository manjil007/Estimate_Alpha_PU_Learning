import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.metrics import matthews_corrcoef, accuracy_score, brier_score_loss, f1_score, average_precision_score
from collections import Counter
from sklearn.utils import shuffle
from sklearn.metrics import confusion_matrix
from scipy.sparse import csr_matrix
from PULSNAR import PULSNAR
import sys
import os
import argparse
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy import stats


def fetch_gene_labels(ifile=None):
    """
    fetch gene CUI and labels from the file. create a dictionary with CUI as key and label as value
    """
    # read the input file
    if ifile is not None:
        df = pd.read_csv(ifile)
    else:
        print("input file1 is missing")
        exit(-1)
    df.columns = df.columns.str.replace(' ', '')

    cui_label_dict = {v.strip().strip('"'): df['gene_label'][i].strip().strip('"') if df['gene_label'][i].strip().strip('"') != 'NULL' else
    v.strip().strip('"').split(' CUI')[0] for i, v in enumerate(df['gene_CUI'])}

    return cui_label_dict


def generate_ml_data(cui_name_dict, ifile=None):
    """
    read input file and generate ML data
    """
    # all genes
    if ifile is not None:
        df = pd.read_csv(ifile)
    else:
        print("input file1 is missing")
        exit(-1)
    df.columns = df.columns.str.replace(' ', '')

    # get all positives, unlabeled and features from the dataframe
    pos_genes_cui = [p.strip(' ').strip('"') for p in df['positive_genes_cui']]
    all_genes_cui = [p.strip(' ').strip('"') for p in df['unknown_genes_cui']]
    all_features = [p.strip(' ').strip('"') for p in df['features']]
    features_to_remove = set([p.strip(' ').strip('"') for p in df['features_to_remove']])

    # select features for all genes
    gene_features_dict = {}
    for i, u in enumerate(all_genes_cui):
        p = pos_genes_cui[i]
        f = all_features[i]
        if f not in features_to_remove:
            gene_features_dict.setdefault(p, set()).add(f)
            gene_features_dict.setdefault(u, set()).add(f)

    # determine unique genes and features
    all_gene_list = sorted(list(gene_features_dict.keys()))
    all_feature_list = np.sort(list(set().union(*gene_features_dict.values())))
    print("Unique genes and features count: ", len(all_gene_list), len(all_feature_list))

    # generate binary features for each gene
    all_gene_labels = []
    all_gene_names = []
    all_gene_features = np.zeros((len(all_gene_list), len(all_feature_list)), dtype=np.uint8)

    pos_genes_cui = set(pos_genes_cui)
    for j, g in enumerate(all_gene_list):
        all_gene_names.append(cui_name_dict[g])
        f = list(gene_features_dict[g])
        idx = np.isin(all_feature_list, f).nonzero()[0]
        all_gene_features[j][idx] = 1
        if g in pos_genes_cui:
            all_gene_labels.append(1)
        else:
            all_gene_labels.append(0)

    # remove genes without any features
    gene_feature_count = np.sum(all_gene_features, axis=1)
    # print(min(gene_feature_count[p_idx]), max(gene_feature_count[p_idx]), np.mean(gene_feature_count[p_idx]), np.median(gene_feature_count[p_idx]))
    ix = np.where(gene_feature_count >= round(np.median(gene_feature_count)))[0]  # to ensure gene has at least one feature
    return csr_matrix(all_gene_features[ix]), np.asarray(all_gene_labels)[ix], np.asarray(all_gene_names)[ix]


def ml_performance_before_pul(X, y, performance_metrics, params, rseed=0):
    """
    Run XGBoost and compute performance
    """
    disp_msg = "before PU"

    # ML model
    X, y = shuffle(X, y, random_state=rseed)
    r = Counter(y)[0] / Counter(y)[1]
    print("Train and test the model (neg/pos): ", r)
    params['scale_pos_weight'] = r
    params['random_state'] = rseed
    bst = XGBClassifier(**params)
    preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

    # performance metrics
    auc_val = roc_auc_score(y, preds[:, 1])
    mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
    acc_val = accuracy_score(y, np.round(preds[:, 1]))
    brier_val = brier_score_loss(y, preds[:, 1])
    f1_val = f1_score(y, np.round(preds[:, 1]))
    aps_val = average_precision_score(y, preds[:, 1])
    tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
    recall_val = tp / (tp + fn)
    ppv_val = tp / (tp + fp)

    # update dictionary
    performance_metrics['AUC'].append(auc_val)
    performance_metrics['MCC'].append(mcc_val)
    performance_metrics['Accuracy'].append(acc_val)
    performance_metrics['Brier_loss'].append(brier_val)
    performance_metrics['F1_score'].append(f1_val)
    performance_metrics['APS'].append(aps_val)
    performance_metrics['Sensitivity'].append(recall_val)
    performance_metrics['Precision'].append(ppv_val)

    # display results
    print(f"\nAUC-ROC {disp_msg}: ", auc_val)
    print(f"MCC {disp_msg}: ", mcc_val)
    print(f"Accuracy {disp_msg}: ", acc_val)
    print(f"Brier score loss {disp_msg}: ", brier_val)
    print(f"F1 score {disp_msg}: ", f1_val)
    print(f"APS {disp_msg}: ", aps_val)
    print(f"sensitivity (recall) {disp_msg}: ", recall_val)
    print(f"precision (PPV) {disp_msg}: ", ppv_val)

    return performance_metrics


def run_pu_classifier(X, y, genes, rseed=0):
    """
    Run PU method to estimate alpha
    """
    # update pulsnar_args.yaml if you want to override the default parameters
    if len(sys.argv) < 2:
        user_param_file = '/home/mpradhan/Intern_Research_Project/data/yaml_files/pulsnar_residues_alpha.yaml'
    else:
        user_param_file = sys.argv[1]

    # check if results folder exist. if not, create it
    if not os.path.exists("results"):
        os.makedirs("results")

    # instantiate PULSNARClassifier
    pls = PULSNAR.PULSNARClassifier(scar=False, csrdata=False, classifier='xgboost',
                                        n_clusters=0, max_clusters=10, covar_type='full', top50p_covars=False,
                                        bin_method='rice', bw_method='hist', lowerbw=0.015, upperbw=0.5, optim='local',
                                        calibration=False, calibration_data='PU', calibration_method='isotonic',
                                        calibration_n_bins=100, smooth_isotonic=False,
                                        classification_metrics=False,
                                        n_iterations=1, kfold=5, kflips=1,
                                        pulsnar_params_file=user_param_file)

    # get results
    X, y, genes = shuffle(X, y, genes, random_state=rseed)
    res = pls.pulsnar(X, y, tru_label=y, rec_list=genes)
    # print(res)
    print("\nEstimated alpha: {0}".format(res['estimated_alpha']))
    return res['estimated_alpha']


def ml_performance_after_pul(X, y, genes, performance_metrics, params, alpha, unlab_gene_preds, rseed=0, preds_file=None):
    """
    Run XGBoost and compute performance
    """
    alpha = 0.2
    disp_msg = "after PU"
    # find probable positives among unlabeled examples
    df = pd.read_csv(preds_file, sep="\t", header=0)
    pred_genes = df['rec_id'].to_numpy()
    pred_calibrated_prob = df['calibrated_prob'].to_numpy()
    pred_predicted_prob = df['predicted_prob'].to_numpy()

    # added predictions to dictionary
    for j, rec in enumerate(pred_genes):
        unlab_gene_preds.setdefault(rec, {})
        unlab_gene_preds[rec].setdefault('calibrated_prob', []).append(pred_calibrated_prob[j])
        unlab_gene_preds[rec].setdefault('predicted_prob', []).append(pred_predicted_prob[j])

    # select top alpha*U unlabeled examples as probable positives
    # idx = np.where(pred_calibrated_prob > 0.5)[0]
    idx = np.argsort(pred_calibrated_prob)[::-1][:int(alpha*len(pred_calibrated_prob))] # select top alpha*U unlabeled
    probable_pos_genes = pred_genes[idx]
    print("number of probable positives: ", len(probable_pos_genes))

    # flip labels for probable positives
    idx = np.isin(genes, probable_pos_genes).nonzero()[0]
    print("number of true positives: ", np.sum(y))
    y[idx] = 1
    print("number of true + probable positives: ", np.sum(y))

    # ML model
    X, y, genes = shuffle(X, y, genes, random_state=rseed)
    r = Counter(y)[0] / Counter(y)[1]
    print("Train and test the model (neg/pos): ", r)
    params['scale_pos_weight'] = r
    params['random_state'] = rseed
    bst = XGBClassifier(**params)
    preds = cross_val_predict(bst, X, y, cv=5, method='predict_proba')

    # performance metrics
    auc_val = roc_auc_score(y, preds[:, 1])
    mcc_val = matthews_corrcoef(y, np.round(preds[:, 1]))
    acc_val = accuracy_score(y, np.round(preds[:, 1]))
    brier_val = brier_score_loss(y, preds[:, 1])
    f1_val = f1_score(y, np.round(preds[:, 1]))
    aps_val = average_precision_score(y, preds[:, 1])
    tn, fp, fn, tp = confusion_matrix(y, np.round(preds[:, 1])).ravel()
    recall_val = tp / (tp + fn)
    ppv_val = tp / (tp + fp)

    # update dictionary
    performance_metrics['AUC'].append(auc_val)
    performance_metrics['MCC'].append(mcc_val)
    performance_metrics['Accuracy'].append(acc_val)
    performance_metrics['Brier_loss'].append(brier_val)
    performance_metrics['F1_score'].append(f1_val)
    performance_metrics['APS'].append(aps_val)
    performance_metrics['Sensitivity'].append(recall_val)
    performance_metrics['Precision'].append(ppv_val)

    # display results
    print(f"\nAUC-ROC {disp_msg}: ", auc_val)
    print(f"MCC {disp_msg}: ", mcc_val)
    print(f"Accuracy {disp_msg}: ", acc_val)
    print(f"Brier score loss {disp_msg}: ", brier_val)
    print(f"F1 score {disp_msg}: ", f1_val)
    print(f"APS {disp_msg}: ", aps_val)
    print(f"sensitivity (recall) {disp_msg}: ", recall_val)
    print(f"precision (PPV) {disp_msg}: ", ppv_val)

    return performance_metrics, unlab_gene_preds


def generate_bar_chart(before_metrics, after_metrics, metric_list, plt_file="metrics_comparison.png"):
    """
    generate barcharts using the performance metrics
    """
    # Plotting
    fig, ax = plt.subplots()

    # Set positions for the bars
    x = np.arange(len(metric_list))
    bar_width = 0.20

    # Plot bars with error bars
    ax.bar(x, calculate_mean_performance(before_metrics, metric_list), bar_width, label='XGBoost only',
           yerr=calculate_95ci(before_metrics, metric_list), capsize=2)
    ax.bar(x + bar_width, calculate_mean_performance(after_metrics, metric_list), bar_width, label='XGBoost+PULSCAR',
           yerr=calculate_95ci(after_metrics, metric_list), capsize=2)

    # Add labels, title, and legend
    ax.set_xlabel('ML model performance metrics', fontsize=12)
    ax.set_ylabel('Scores', fontsize=12)
    # ax.set_title('Comparison of metrics with/without PULSCAR', fontsize=13, fontweight='bold')
    ax.set_xticks(x + bar_width / 2)
    ax.set_xticklabels(metric_list, rotation=90)
    ax.legend(loc='lower left')
    # ax.legend()
    # ax.grid(visible=True, which='major', axis='y', color='0.25', linestyle='-')
    # ax.grid(visible=True, which='minor', axis='both', color='0.45', linestyle='--')
    ax.minorticks_on()
    ax.xaxis.set_tick_params(which='minor', bottom=False)

    # Save plot as a high-quality PNG file
    plt.savefig(plt_file, dpi=300, bbox_inches='tight')

    # Show plot
    plt.show()


def calculate_mean_performance(performance_metrics, metric_list):
    """
    calculate mean for each of the performance metric
    """
    mean_values = []
    for m in metric_list:
        mean_values.append(np.mean(performance_metrics[m]))
    return mean_values


def calculate_95ci(performance_metrics, metric_list):
    """
    calculate 95% CI for each of the performance metric
    """
    ci_values = []
    for m in metric_list:
        data = performance_metrics[m]
        ci_95p = stats.t.interval(0.95, len(data) - 1, loc=np.mean(data), scale=stats.sem(data))
        ci_values.append([np.mean(data) - ci_95p[0], ci_95p[1] - np.mean(data)])
    return np.transpose(ci_values)

def save_mean_predictions(gene_preds, opfile=None):
    """
    save the mean predictions for all unlabeled genes
    """
    with open(opfile, 'w') as fo:
        hdr = "gene\tmean_predicted_probs\tmean_calibrated_probs\n"
        fo.write(hdr)
        # write predictions for genes
        for gene, preds in gene_preds.items():
            line = str(gene) + "\t" + str(np.mean(preds['predicted_prob'])) + "\t" + str(np.mean(preds['calibrated_prob'])) + "\n"
            fo.write(line)

def main():
    """
    This code calls PU classifier to estimate the alpha. Then using true positives, probable positives, and probable negatives, it runs XGBoost
    model and calculates classification performance metrics. It also calculates classification performance metrics without running the PU classifier.
    In this code, after metrics are based on class 1 = labeled positives + imputed positives and class 0 = imputed negatives
    """
    # local variables to save values
    n_iterations = 1
    alphas = []
    unlab_gene_preds = {}
    before_metrics = {'AUC': [], 'MCC': [], 'Accuracy': [], 'Brier_loss': [], 'F1_score': [], 'APS': [], 'Sensitivity': [], 'Precision': []}
    after_metrics = {'AUC': [], 'MCC': [], 'Accuracy': [], 'Brier_loss': [], 'F1_score': [], 'APS': [], 'Sensitivity': [], 'Precision': []}
    metric_list = ['AUC', 'MCC', 'Accuracy', 'Brier_loss', 'F1_score', 'APS', 'Sensitivity', 'Precision']
    xgb_params = {'max_depth': 4, 'n_jobs': 16, 'eval_metric': 'logloss'}

    # get command line parameters
    parser = argparse.ArgumentParser()
    # parser.add_argument("-inpfile", help="provide input file (output of CQL)", required=True)
    # parser.add_argument("-gene_label_file", help="provide input file containing genes and their labels", required=True, default="IOData/gene_label_data.csv")
    parser.add_argument("-pltfile", help="provide filename to save bar plot", required=True, default="IOData/results/temp.png")
    parser.add_argument("-predfile", help="provide filename to save predictions", required=True, default="IOData/results/preds.tsv")
    args = parser.parse_args()

    # fetch labels for all genes
    # orig_cui_name_dict = fetch_gene_labels(ifile=args.gene_label_file)

    # generate positive and unknown data
    # orig_X, orig_y, orig_genes = generate_ml_data(orig_cui_name_dict, ifile=args.inpfile)

    orig_X = pd.read_csv("/home/mpradhan/Intern_Research_Project/data/X_balanced_validation.csv").to_numpy()
    orig_y = pd.read_csv("/home/mpradhan/Intern_Research_Project/data/y_balanced_true_validation.csv").to_numpy().flatten()
    orig_genes = np.arange(len(orig_y))

    print("Data shape, positive, unlabeled: ", orig_X.shape, len(orig_y), np.sum(orig_y), len(orig_genes) - np.sum(orig_y))

    for itr in range(n_iterations):
        print("\nRunning models for iteration: {0}".format(itr + 1))

        # Classification performance before applying PUL
        X, y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        before_metrics = ml_performance_before_pul(X, y, before_metrics, xgb_params, rseed=itr)  # all data

        # run PU method
        X, y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        est_alpha = run_pu_classifier(X, y, genes, rseed=itr)
        alphas.append(est_alpha)

        # classification performance after applying PUL
        X, y, genes = deepcopy(orig_X), deepcopy(orig_y), deepcopy(orig_genes)
        after_metrics, unlab_gene_preds = ml_performance_after_pul(X, y, genes, after_metrics, xgb_params, est_alpha, unlab_gene_preds, rseed=itr,
                                                                   preds_file="predictions.tsv")

    print("\nbefore metrics: ", before_metrics)
    print("\nafter metrics: ", after_metrics)
    print("\nalphas: ", alphas)

    # generate barchart using performance metrics
    print("\nmean alpha and 95% CI: ", np.mean(alphas), stats.t.interval(0.95, len(alphas) - 1, loc=np.mean(alphas), scale=stats.sem(alphas)))
    print("\nbefore Recall and 95% CI: ", np.mean(before_metrics['Sensitivity']), stats.t.interval(0.95, len(before_metrics['Sensitivity']) - 1, loc=np.mean(before_metrics['Sensitivity']), scale=stats.sem(before_metrics['Sensitivity'])))
    print("\nafter Recall and 95% CI: ", np.mean(after_metrics['Sensitivity']), stats.t.interval(0.95, len(after_metrics['Sensitivity']) - 1, loc=np.mean(after_metrics['Sensitivity']), scale=stats.sem(after_metrics['Sensitivity'])))

    generate_bar_chart(before_metrics, after_metrics, metric_list, plt_file=args.pltfile)

    # save mean predictions
    save_mean_predictions(unlab_gene_preds, opfile=args.predfile)

if __name__ == "__main__":
    main()