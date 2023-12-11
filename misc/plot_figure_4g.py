import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
import numpy as np
import argparse
import os

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-in', dest = 'in_path', required=True, help='Path to csv file containing normalized single-cell data.')
    parser.add_argument('-out', dest = 'out_path', required=True, help='Path to output folder.')
    args = parser.parse_args()
    return args


def split(word: str):
    return [char for char in word]

def main (args):

    in_path = args.in_path
    out_path = args.out_path

    df = pd.read_csv(in_path, delimiter=';', decimal=',')
    col = [i for i in df.columns if 'MXP_c' in i or 'batch' in i or 'cluster' in i]

    crtl_list = ['BM1.2', 'BM1.3', 'BM2.3']
    inf_list = ['BM2.1', 'BM2.2', 'BM4.1', 'BM1.1', 'BM3.1']

    df = df[col]
    df = df[df['cluster']=='MO/M']
    df['inf'] = ['crtl' if df.iloc[i]['batch'] in crtl_list else 'NB' for i in range(0, df.shape[0])]
    df = df.drop(['batchlegend', 'batch', 'cluster'], axis=1)

    col_list = list(df.columns)
    col_list.remove('inf')

    df.inf = df.inf.astype("category")

    for marker in col_list:

        medians = df.groupby(['inf'])[marker].median()
        medians = medians.reindex(index=['crtl', 'NB'], copy=True)
        vertical_offset = df[marker].median() * 0.05 

        ax = sns.boxplot(x='inf', y=marker, data=df, linewidth=0.4, fliersize=1, width=0.5, order=['crtl', 'NB'])

        plot = add_stat_annotation(ax, data=df, x='inf', y=marker, box_pairs=[('crtl', 'NB')], test='Mann-Whitney', text_format='star', loc='outside', verbose=2)

        for xtick in ax.get_xticks():
            ax.text(xtick,medians[xtick] + vertical_offset,medians[xtick], horizontalalignment='center',size='x-small',color='black',weight='semibold')

        plt.tight_layout()
        plt.savefig(os.path.join(out_path, marker +'_expr_crtl_vs_NB.eps'), format='eps')
        plt.savefig(os.path.join(out_path, marker +'_expr_crtl_vs_NB.pdf'), format='pdf')
        plt.close()


if __name__ == "__main__":
        args = parse()
        main(args)