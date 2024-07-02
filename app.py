import csv
import math
import os

import gradio as gr

from volcano import volcano_plot


def upload_file(
    filepath: str,
) -> list[
    gr.UploadButton,
    gr.CheckboxGroup,
    gr.Dropdown,
    gr.Textbox,
    gr.Number,
    gr.Number,
    gr.Button,
]:
    '''Function for uploading file.

    Args:
        filepath (str): Temporary path to the uploaded file.

    Returns:
        list[gradio.UploadButton | gradio.CheckboxGroup | gradio.Dropdown | gradio.Textbox | gradio.Number | gradio.Button]: List of the gradio modules to visualize.
    '''
    fname = os.path.basename(filepath)

    with open(filepath, 'r') as f:
        header, *data = list(csv.reader(f, delimiter='\t'))

    if 'id' in header:
        gr.Error('"id" header does not exist.')

    if 'neg|lfc' in header:
        gr.Error('"neg|lfc" header does not exist.')

    if 'neg|p-value' in header:
        gr.Error('"neg|p-value" header does not exist.')

    if 'pos|lfc' in header:
        gr.Error('"pos|lfc" header does not exist.')

    if 'pos|p-value' in header:
        gr.Error('"pos|p-value" header does not exist.')

    idx_id = header.index('id')
    genes = sorted([sample[idx_id] for sample in data])
    genes.insert(0, '<all_significant>')

    return [
        gr.UploadButton(label=f'Uploaded {fname}'),
        gr.CheckboxGroup(visible=True),
        gr.Dropdown(choices=genes, value='<all_significant>', visible=True),
        gr.Textbox(visible=True),
        gr.Number(value=1.0, visible=True),
        gr.Number(
            value=round(-1 * math.log10(0.05), 2),
            visible=True,
        ),
        gr.Dropdown(visible=True),
        gr.Button(visible=True),
    ]


def draw_plot(
    filepath: str,
    src_data: list[str],
    genes_label: list[str],
    title: str,
    thr_x: str,
    thr_y: str,
    fmt_save: str,
) -> list[gr.Image, gr.Dropdown, gr.DownloadButton]:
    '''Draw the volcano plot.

    Args:
        filepath (str): Path to analysis result TXT file from Mageck.
        src_data (list[str]): Data source to visualize. Options are "Negative" and "Positive".
        genes_label (list[str]): List of the genes to label.
        title (str): Title of the figure.
        thr_x (str): Threshold of log2 fold change for interesting gene.
        thr_y (str): Threshold of -log10 P-value for interesting gene.
        fmt_save (str): Format of figure to save. If not given, figure will not be saved.

    Returns:
        list[gr.Image, gr.DownloadButton]: List of the gradio modules to update.
    '''
    img = volcano_plot(
        filepath=filepath,
        src_data=src_data,
        genes_label=genes_label,
        title=title,
        thr_x=thr_x,
        thr_y=thr_y,
        fmt_save=None,
    )
    path_save = volcano_plot(
        filepath=filepath,
        src_data=src_data,
        genes_label=genes_label,
        title=title,
        thr_x=thr_x,
        thr_y=thr_y,
        fmt_save=fmt_save,
    )
    fname = os.path.basename(path_save)

    return [
        gr.Image(value=img, visible=True),
        gr.DownloadButton(
            label=f'Download {fname}', value=path_save, visible=True
        ),
    ]


with gr.Blocks() as app:
    with gr.Tab('Volcano Plot'):
        with gr.Row():
            with gr.Column():
                ulbtn_data = gr.UploadButton(
                    label='Upload', file_count='single'
                )
                ckb_src_data = gr.CheckboxGroup(
                    ['Negative', 'Positive'],
                    label='Data source to visualize',
                    interactive=True,
                    visible=False,
                )
                dpdn_gene_label = gr.Dropdown(
                    label='Genes to label',
                    multiselect=True,
                    interactive=True,
                    visible=False,
                )
                txt_title = gr.Textbox(
                    label='Title of the figure (optional)',
                    interactive=True,
                    visible=False,
                )
                nmb_thr_x = gr.Number(
                    label='Threshold of log2 fold change for interesting gene',
                    interactive=True,
                    visible=False,
                )
                nmb_thr_y = gr.Number(
                    label='Threshold of -log10 P-value for interesting gene',
                    interactive=True,
                    visible=False,
                )
                dpdn_fmt = gr.Dropdown(
                    choices=[('PNG', 'png'), ('PDF', 'pdf'), ('SVG', 'svg')],
                    value='svg',
                    label='The format of figure to download',
                    interactive=True,
                    visible=False,
                )
                btn_draw = gr.Button(
                    value='Draw Plot', interactive=True, visible=False
                )

            with gr.Column():
                img_fig = gr.Image(visible=False)
                dlbtn_fig = gr.DownloadButton(visible=False)

        ulbtn_data.upload(
            upload_file,
            ulbtn_data,
            [
                ulbtn_data,
                ckb_src_data,
                dpdn_gene_label,
                txt_title,
                nmb_thr_x,
                nmb_thr_y,
                dpdn_fmt,
                btn_draw,
            ],
        )
        btn_draw.click(
            draw_plot,
            [
                ulbtn_data,
                ckb_src_data,
                dpdn_gene_label,
                txt_title,
                nmb_thr_x,
                nmb_thr_y,
                dpdn_fmt,
            ],
            [img_fig, dlbtn_fig],
        )


if __name__ == '__main__':
    app.launch(server_port=8000)
