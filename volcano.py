import argparse
import csv
import io
import math
import os

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from PIL import Image

COLOR = {'UP': 'red', 'DOWN': 'blue', 'NOT SIG': 'grey'}


def volcano_plot(
    filepath: str,
    src_data: list[str],
    genes_label: list[str],
    title: str,
    thr_x: float,
    thr_y: float,
    fmt_save: str,
) -> str | Image.Image:
    '''Plot volcano plot of gene analysis results.

    Args:
        filepath (str): Path to analysis result TXT file from Mageck.
        src_data (list[str]): Data source to visualize. Options are "Negative" and "Positive".
        genes_label (list[str]): List of the genes to label.
        title (str): Title of the figure.
        thr_x (float): Threshold of log2 fold change for interesting gene.
        thr_y (float): Threshold of -log10 P-value for interesting gene.
        fmt_save (str): Format of figure to save. If not given, figure will not be saved.

    Returns:
        str | PIL.Image.Image: If save_fmt is not given, a PIL.Image.Image object of the figure is returned. Otherwise, the path to the temporary output figure file.
    '''
    # Data preparation
    with open(filepath, 'r') as f:
        header, *data = list(csv.reader(f, delimiter='\t'))

    idx_x_neg = header.index('neg|lfc')
    idx_y_neg = header.index('neg|p-value')
    idx_x_pos = header.index('pos|lfc')
    idx_y_pos = header.index('pos|p-value')
    idx_id = header.index('id')

    up = []
    up_id = []
    down = []
    down_id = []
    not_sig = []
    not_sig_id = []

    for sample in data:
        id = sample[idx_id]

        if 'Negative' in src_data:
            x_neg = float(sample[idx_x_neg])
            y_neg = -1 * math.log10(float(sample[idx_y_neg]))

            if y_neg > thr_y and x_neg < -thr_x:
                down.append([x_neg, y_neg])

                if (id in genes_label) or ('<all_significant>' in genes_label):
                    down_id.append(id)
                else:
                    down_id.append(None)
            else:
                not_sig.append([x_neg, y_neg])

                if id in genes_label:
                    not_sig_id.append(id)
                else:
                    not_sig_id.append(None)

        if 'Positive' in src_data:
            x_pos = float(sample[idx_x_pos])
            y_pos = -1 * math.log10(float(sample[idx_y_pos]))

            if y_pos > thr_y and x_pos > thr_x:
                up.append([x_pos, y_pos])

                if (id in genes_label) or ('<all_significant>' in genes_label):
                    up_id.append(id)
                else:
                    down_id.append(None)
            else:
                not_sig.append([x_pos, y_pos])

                if id in genes_label:
                    not_sig_id.append(id)
                else:
                    not_sig_id.append(None)

    up = np.array(up)
    down = np.array(down)
    not_sig = np.array(not_sig)

    fig, ax = plt.subplots()
    ax.axhline(thr_y, c='k', ls='--')
    ax.axvline(thr_x, c='k', ls='--')
    ax.axvline(-thr_x, c='k', ls='--')

    text = []

    if len(up):
        ax.scatter(up[:, 0], up[:, 1], c='red', alpha=0.5, label='UP')

        for i in range(len(up)):
            if up_id:
                text.append(
                    ax.text(x=up[i, 0], y=up[i, 1], s=up_id[i], fontsize=8)
                )

    if len(down):
        ax.scatter(down[:, 0], down[:, 1], c='blue', alpha=0.5, label='DOWN')

        for i in range(len(down)):
            if down_id:
                text.append(
                    ax.text(
                        x=down[i, 0], y=down[i, 1], s=down_id[i], fontsize=8
                    )
                )

    if len(not_sig):
        ax.scatter(
            not_sig[:, 0], not_sig[:, 1], c='grey', alpha=0.5, label='NOT SIG'
        )

        for i in range(len(not_sig)):
            if not_sig_id:
                text.append(
                    ax.text(
                        x=not_sig[i, 0],
                        y=not_sig[i, 1],
                        s=not_sig_id[i],
                        fontsize=8,
                    )
                )

    if text:
        adjust_text(text, arrowprops=dict(arrowstyle='-', color='k'))

    if title:
        fig.title(title)

    ax.set_xlabel("$log_{2}$ fold change", size=15)
    ax.set_ylabel("-$log_{10}$ P value", size=15)

    if fmt_save:
        fname_save = (
            os.path.basename(filepath).rsplit('.', 1)[0] + f'.{fmt_save}'
        )
        path_save = os.path.join(os.path.dirname(filepath), fname_save)
        fig.savefig(path_save)

        return path_save
    else:
        buf = io.BytesIO()
        fig.savefig(buf)
        buf.seek(0)

        return Image.open(buf)
