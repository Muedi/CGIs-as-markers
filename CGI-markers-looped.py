"""
There will be a corresponding R script and environment
The R script will read the idat or GEO-series files and return methylation values
This script will make use of these. The meth values will be used for training a
RF classifier and employing different feature selection to the data beforehand
"""
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve
from numpy import random as rd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.feature_selection import f_classif, SelectFromModel, SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from scipy.stats import randint as sp_randint

from sklearn.model_selection import RandomizedSearchCV
from sklearn import metrics
from tabulate import tabulate


def ROC_Curve_new(y_predicted_test, Y_test, auc):
    # plots a crude ROC-Curve
    false_positive, true_positive, _ = roc_curve(Y_test.ravel(), y_predicted_test.ravel())
    plot = plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(false_positive, true_positive, color='darkorange', label='Random Forest')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve (area = %0.2f)' % auc)
    plt.legend(loc='best')
    plt.show()


def Print_Metrics(Y_test, y_predict_test):
    # Define a method to print the model performance
    print('\nModel performance on the test data set:')

    logloss_test = metrics.log_loss(Y_test, y_predict_test)
    accuracy_test = metrics.accuracy_score(Y_test, y_predict_test)
    F1_test = metrics.f1_score(Y_test, y_predict_test)
    precision_test = metrics.precision_score(Y_test, y_predict_test, average='binary')
    auc_test = metrics.roc_auc_score(Y_test, y_predict_test)
    r2_test = metrics.r2_score(Y_test, y_predict_test)

    header = ["Metric", "Test"]
    table  = [
               ["logloss",   logloss_test],
               ["accuracy",  accuracy_test],
               ["precision", precision_test],
               ["F1",        F1_test],
               ["r2",        r2_test],
               ["AUC",       auc_test]
             ]

    print(tabulate(table, header, tablefmt="fancy_grid"))


def Plot_predictor_importance(best_model, feature_columns):
    # Method to plot the predictor importance
    feature_importances = best_model.feature_importances_
    sorted_idx = np.argsort(feature_importances)
    y_pos  = np.arange(sorted_idx.shape[0]) + .5
    fig, ax = plt.subplots()
    ax.barh(y_pos,
            feature_importances[sorted_idx],
            align='center',
            color='green',
            ecolor='black',
            height=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(feature_columns)
    ax.invert_yaxis()
    ax.set_xlabel('Relative Importance')
    ax.set_title('Predictor Importance')
    plt.show()
    return fig


def Print_confusion_matrix(cm, auc, heading):
    # Define a method to print the Confusion Matrix and the performance metrics
    print('\n', heading)
    print(cm)
    true_negative  = cm[0,0]
    true_positive  = cm[1,1]
    false_negative = cm[1,0]
    false_positive = cm[0,1]
    total = true_negative + true_positive + false_negative + false_positive
    accuracy = (true_positive + true_negative)/total
    precision = (true_positive)/(true_positive + false_positive)
    recall = (true_positive)/(true_positive + false_negative)
    misclassification_rate = (false_positive + false_negative)/total
    F1 = (2*true_positive)/(2*true_positive + false_positive + false_negative)
    print('accuracy.................%7.4f' % accuracy)
    print('precision................%7.4f' % precision)
    print('recall...................%7.4f' % recall)
    print('F1.......................%7.4f' % F1)
    print('auc......................%7.4f' % auc)


def get_islands_Rnbeads(index_file, meth_site_file):
    """
    opens the betas from an index and a value csv file.
    betas were extracted with the RnBeads package for R
    the files are returned in the same order, so that the islands can be taken
    as index for the meth values
    """
    island_anno = pd.read_csv(index_file)
    beta_means = pd.read_csv(meth_site_file)
    beta_means["UCSC_CpG_Islands_Name"] = (island_anno["Chromosome"].astype("str") + ":" + island_anno["Start"].astype("str") + "-" + island_anno["End"].astype("str")).values
    beta_means = beta_means[beta_means != 'NA'].dropna()
    beta_means.set_index(['UCSC_CpG_Islands_Name'], inplace = True)
    return beta_means.astype("float")


# Anova feature selection
def anova_k_select(X_train, Y_train, k=10):
    """
    univariant anova for feature selection
    returns df with selected features and corresponding data in the same form as the input
    """
    select = SelectKBest(f_classif, k)
    select.fit(X_train, Y_train)
    X_train_selected = select.transform(X_train)
    index = X_train.T[select.get_support()].T.index
    columns = X_train.T[select.get_support()].index
    ret_def = pd.DataFrame(X_train_selected, index = index, columns = columns)
    return ret_def


# lasso deature selection
def lasso_selection(X_train, Y_train, features=5, threshold=0.1):
    """ lasso with CV is used to get a linear combination witb l1 penalty,
    meaning, that non-important features will get smaller coefficients, promoting sparsety in the set
    one can pick the remaining features with big coefficients above a certain threshold
    """
    random_state=rd.RandomState(0)
    lasso_clf = LogisticRegression(random_state=random_state, penalty="l1")
    select = SelectFromModel(lasso_clf, threshold=threshold)
    select.fit(X_train, Y_train)
    X_train_selected = select.transform(X_train)
    n_features = X_train_selected.shape[1]

    # Reset the threshold untill the number of features equals features-value passed to the function.
    # fitting the metatransformer.
    while n_features > features:
        select.threshold += 0.005
        # select.threshold += 0.1
        X_train_selected = select.transform(X_train)
        n_features = X_train_selected.shape[1]

    # return dataframe in form of input X_train
    index = X_train.T[select.get_support()].T.index
    columns = X_train.T[select.get_support()].index
    ret_def = pd.DataFrame(X_train_selected, index = index, columns = columns)
    return ret_def


# MI selection
def mutual_info_k_select(X_train, Y_train, k=10):
    """
    Estimate mutual information for a discrete target variables
    returns df with selected features and corresponding data in the same form as the input
    """
    select = SelectKBest(mutual_info_classif, k)
    select.fit(X_train, Y_train)
    X_train_selected = select.transform(X_train)
    index = X_train.T[select.get_support()].T.index
    columns = X_train.T[select.get_support()].index
    ret_def = pd.DataFrame(X_train_selected, index = index, columns = columns)
    return ret_def


# Select if Sites or islands should be used.
sites = True
if sites:
# loads the normalized data from RNBeads sites
    meth_site_file = "Output_RnBeads/complete_HCC_sites_meth.csv"
    data_sites = pd.read_csv(meth_site_file)
    data_sites.dropna(inplace=True)
    data_sites = data_sites.T
    data_sites["Class"] = ["HCC" if ("HCC" in x or "carcinoma" in x) else ("blood") for x in data_sites.index]
    data_sites["Class_code"] = [0 if x == "blood" else 1 for x in data_sites["Class"]]

    data=data_sites

else:
    # loads the normalized data from RNBeads  CpG-islands
    index_file = "Output_RnBeads/complete_HCC_islands_anno.csv"
    meth_file = "Output_RnBeads/complete_HCC_islands_meth.csv"
    data = get_islands_Rnbeads(index_file, meth_file)
    data = data.T
    data["Class"] = ["HCC" if ("HCC" in x or "carcinoma" in x) else ("blood") for x in data.index]
    data["Class_code"] = [0 if x == "blood" else 1 for x in data["Class"]]


data.drop(["PBMC_25", "PBMC_26"], inplace=True)

# init random state
seed = rd.randint(1,400)
random_state = rd.RandomState(seed)
# random_state = rd.RandomState(0)

# data processing
X = data.drop(['Class'], axis=1).iloc[:, :-1]
Y = data.drop(['Class'], axis=1).iloc[:, -1]
# will be used for collecting selected features:
all_dem_features = pd.DataFrame([0] * len(X.columns), index = X.columns)
# array that contains the names of three used feature selections
selectors = ["ANOVA", "MI", "LASSO"]
# change to true, if the complete set should be used and a train/test
# split of 0.6/0.4 should be applied to the dataset:
complete_set = False

for selector in selectors:
    mis_pred_samples = pd.DataFrame([0]*len(data.index), columns=["mis-prediction"], index=data.index )
    sum_mse = []
    sum_accuracY_test = []
    sum_F1_test = []
    sum_precision_test = []
    sum_recall_test = []
    sum_auc_test = []
    sum_r2_test = []
    # Number of desired iterations
    for i in range(10):
        if complete_set:
            seed = rd.randint(1,400)
            random_state = rd.RandomState(seed)
            X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.4, random_state=random_state)
        else:
            X_train = X[[True if ("PBMC" in x or "carcinoma" in x or "HCC_" in x) else False for x in X.index]]
            X_train = X_train.sample(frac=1)  # shuffle dataframe
            Y_train = Y.loc[X_train.index]
            # x_test with only cfDNA and cfHCC data
            X_test = X[[True if ("cf" in x) else False for x in X.index]]
            Y_test = Y.loc[X_test.index]
            # decrease X to the samples in training/test does nothing,
            # if the training and test sets have all samples
            # if you drop the whole blood samples, these are removed from X to
            # prevent
            index = list(X_train.index.values) + list(X_test.index.values)
            X = X.loc[index]
            Y = Y.loc[X.index]
        # parameters to work with durich random search
        param_grid = {"n_estimators": range(100, 500, 2),
                      "max_depth": range(4, 50, 2),
                      "min_samples_leaf": range(2, 100, 2),
                      "min_samples_split": sp_randint(2, 10),
                      "bootstrap": [True, False],
                      "criterion": ["gini", "entropy"]}
        rf_classifier = RandomForestClassifier(class_weight = 'balanced', random_state=random_state)
        # number of iterations for the search
        n_iter_search = 500
        estimator = RandomizedSearchCV(rf_classifier,
                                       param_distributions = param_grid,
                                       n_iter = n_iter_search,
                                       scoring = 'roc_auc',
                                       verbose = 0,
                                       n_jobs = 1,
                                       random_state=random_state)
        #feature selections:
        if selector == "ANOVA":
            X_train = anova_k_select(X_train, Y_train, 5)
        if selector == "MI":
            X_train = mutual_info_k_select(X_train, Y_train, 5)
        if selector == "LASSO":
            X_train = lasso_selection(X_train, Y_train, 5)

        # xtest is filtered for the selected features
        X_test = X_test[X_train.columns]
        fit = estimator.fit(X_train, Y_train)

        print("\n______________________________________________________________________\n")
        print("Feature selection method: %s Cross-validation fold: %i" % (selector, i))
        print("\n______________________________________________________________________\n")
        best_model = estimator.best_estimator_
        best_model.random_state = random_state
        # prints best model, this could be more than one, if the scorefunctions get the same value
        print('\nbest_model:\n', best_model)

        # print features and importances, plot them with function
        feature_importances = pd.DataFrame(best_model.feature_importances_,index = X_train.columns,columns=['importance']).sort_values('importance', ascending=True)
        # add selected features to the all_dem_features dataframe
        all_dem_features.loc[feature_importances.index] += 1
        # print("\n Selected Features:")
        # print(feature_importances)
        # fig = Plot_predictor_importance(best_model, feature_importances.index)

        # predicts the training set for evaluation
        y_predicted_train = best_model.predict(X_train)
        cm = confusion_matrix(Y_train, y_predicted_train)
        auc = roc_auc_score(np.array(Y_train), y_predicted_train)
        Print_confusion_matrix(cm, auc, 'Confusion matrics of the training dataset')

        # actual prediction of the test set
        y_predicted_test = best_model.predict(X_test)
        ntotal = len(Y_test)
        correct = Y_test == y_predicted_test
        numCorrect = sum(correct)
        percent = round( (100.0*numCorrect)/ntotal, 6)
        print("\nCorrect classifications on test data: {0:d}/{1:d} {2:8.3f}%".format(numCorrect, ntotal, percent))
        prediction_score = 100.0*best_model.score(X_test, Y_test)
        print('Random Forest Prediction Score on test data: %8.3f' % prediction_score)

        cm = confusion_matrix(Y_test, y_predicted_test)
        auc = roc_auc_score(Y_test, y_predicted_test)
        Print_confusion_matrix(cm, auc, 'Confusion matrics of the test dataset')

        plot = ROC_Curve_new(y_predicted_test, Y_test.values, auc)
        ## plot.savefig()

        Print_Metrics(Y_test, y_predicted_test)
        output_df = pd.DataFrame([Y_test.values, y_predicted_test], index=["Class", "Prediction"], columns=Y_test.index ).T
        print(output_df[output_df["Class"] != output_df["Prediction"]])

        # sum of all metrics is represented as an array
        sum_mse.append(metrics.mean_squared_error(Y_test, y_predicted_test))
        sum_accuracY_test.append(metrics.accuracy_score(Y_test, y_predicted_test))
        sum_F1_test.append(metrics.f1_score(Y_test, y_predicted_test))
        sum_precision_test.append(metrics.precision_score(Y_test, y_predicted_test, average='binary'))
        sum_recall_test.append(metrics.recall_score(Y_test, y_predicted_test, average='binary'))
        sum_auc_test.append(metrics.roc_auc_score(Y_test, y_predicted_test))
        sum_r2_test.append(metrics.r2_score(Y_test, y_predicted_test))

        mis_pred_samples.loc[output_df[output_df["Class"] != output_df["Prediction"]].index] += 1


    # print average of metrics
    print("\nSet: ", set, "\t Selector: ", selector, "\n")
    header = ["Metric", "Mean", "Standard Deviation"]
    table  = [
        ["accuracy", np.mean(sum_accuracY_test), np.std(sum_accuracY_test)],
        ["precision", np.mean(sum_precision_test), np.std(sum_precision_test)],
        ["F1",        np.mean(sum_F1_test), np.std(sum_F1_test)],
        ["r2",        np.mean(sum_r2_test), np.std(sum_r2_test)],
        ["AUC",       np.mean(sum_auc_test), np.std(sum_auc_test)],
        ["Recall",    np.mean(sum_recall_test), np.std(sum_recall_test)]
        ]

    print(tabulate(table, header, tablefmt="fancy_grid"))

    print("\n mispredicted Samples:")
    print(mis_pred_samples[mis_pred_samples > 0].dropna())

    print("Selected features:")
    print(all_dem_features[all_dem_features > 0].dropna())

    print("Selected features greater or equal to 5:")
    print(all_dem_features[all_dem_features >= 5].dropna())
