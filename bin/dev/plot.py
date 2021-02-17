import argparse
import matplotlib.pyplot as plt
import numpy as np



def main():
    args = parseCmd()
    label = [   '$\it{Arabidopsis~thaliana}$\nspecies excluded',
                '$\it{Arabidopsis~thaliana}$\nfamily excluded',
                '$\it{Arabidopsis~thaliana}$\norder excluded',
                '$\it{Caenorhabditis~elegans}$\nspecies excluded',
                '$\it{Caenorhabditis~elegans}$\nfamily excluded',
                '$\it{Caenorhabditis~elegans}$\norder excluded',
                '$\it{Danio~rerio}$\norder excluded',
                '$\it{Drosophila~melanogaster}$\nspecies excluded',
                '$\it{Drosophila~melanogaster}$\nfamily excluded',
                '$\it{Drosophila~melanogaster}$\norder excluded',
                '$\it{Medicago~truncatula}$\norder excluded']
    trans_sn = {'braker1' : [47.83, 47.83, 47.83, 48.57, 48.57, 48.57, 27.92, 44.33, 44.33, 44.33, 40.89],
            'braker2' :     [54.69, 50.62, 48.86, 49.15, 35.18, 35.58, 24.91, 50.88, 42.81, 38.89, 45.25],
            'test1' :       [54.56, 51.27, 50.17, 50.98, 45.01, 44.69, 31.74, 50.49, 44.49, 42.50, 47.70]}

    trans_sp = {'braker1' : [57.94, 57.94, 57.94, 58.44, 58.44, 58.44, 24.08, 57.78, 57.78, 57.78, 37.46],
            'braker2' :     [68.39, 66.58, 64.19, 63.02, 53.68, 54.90, 19.56, 67.90, 60.89, 57.99, 41.86],
            'test1' :       [78.37, 78.00, 76.54, 74.52, 70.97, 70.87, 33.88, 78.45, 72.74, 71.10, 54.22]}

    cds_sn= {'braker1' :    [81.70, 81.70, 81.70, 85.21, 85.21, 85.21, 80.52, 77.95, 77.95, 77.95, 78.44],
            'braker2' :     [83.52, 81.77, 81.05, 84.56, 74.77, 75.51, 75.39, 80.16, 74.75, 71.55, 78.92],
            'test1' :       [83.34, 82.4, 82.02, 84.78, 81.98, 81.98, 77.2, 80.1, 76.93, 76.02, 80.03],
            'test2' :       [83.23, 82.37, 82.0, 84.25, 81.57, 81.66, 77.12, 79.64, 76.78, 75.71, 79.72]}

    cds_sp= {'braker1' :    [81.49, 81.49, 81.49, 86.68, 86.68, 86.68, 71.59, 80.20, 80.20, 80.20, 65.6],
            'braker2' :     [86.42, 87.12, 86.33, 90.30, 88.14, 88.56, 69.17, 87.12, 84.74, 82.88, 71.24],
            'test1' :       [87.96, 88.39, 87.98, 91.95, 90.9, 90.9, 70.45, 87.91, 86.19, 85.42, 71.71],
            'test2' :       [90.02, 90.6, 90.34, 92.54, 91.51, 91.51, 77.46, 90.97, 88.98, 88.04, 77.38]}



    tests = {'braker1' : 'tab:orange', 'braker2' : 'tab:red', 'test1' : 'tab:green'}
    x = np.arange(len(label))
    fig, ((ax1, ax2)) = plt.subplots(2, 1)
    fig.subplots_adjust(hspace=.05)
    fig.suptitle('Transcript Accuracy', fontsize=40)
    wid = 0.24
    i = -1.5 * wid
    ax1.set_ylabel('Sensitivity', fontsize=30)
    ax1.set_ylim(0,100)
    plt.xticks(x, label, fontsize=20)
    ax1.set_xticks([],[])
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    for t in tests.keys():
        ax1.bar(x +i, trans_sn[t], color=tests[t], width=wid)
        i += wid
    i = -1.5 * wid
    ax2.set_ylabel('Specificity', fontsize=30)
    ax2.set_ylim(0,100)
    for t in tests.keys():
        ax2.bar(x +i, trans_sp[t], color=tests[t], width=wid)
        i += wid
    fig.legend(["BRAKER1", "BRAKER2", "Decision_Rule"], prop={'size': 30})
    '''
    i = -1.5 * wid
    for t in tests.keys():
        ax3.bar(x +i, cds_sn[t], color=tests[t], width=wid)
        i += wid
    i = -1.5 * wid
    for t in tests.keys():
        ax4.bar(x +i, cds_sp[t], color=tests[t], width=wid)
        i += wid
    i = -1.5 * wid
    '''
    plt.show()
def parseCmd():
    """Parse command line arguments

    Returns:
        dictionary: Dictionary with arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--test_dirs', type=str,
        help='')
    parser.add_argument('--test_label', type=str,
        help='')
    return parser.parse_args()

if __name__ == '__main__':
    main()
