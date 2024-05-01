#!/usr/bin/env python

import argparse

import plotly.graph_objects as go
import numpy as np
import pandas as pd

from plotly.subplots import make_subplots

from RGTools.BedTable import BedTable3, BedTable6, BedTable6Plus

def init_fig(num_annotation_tracks: int) -> go.Figure:
    fig = make_subplots(rows=num_annotation_tracks+1, 
                        cols=1, 
                        shared_xaxes=True, 
                        vertical_spacing=0.01, 
                        row_heights=[1] + [0.1]*num_annotation_tracks,
                        )

    return fig

def make_gwas_scatter(gwas_summary_bed_table: BedTable6, 
                    plot_chr: str, 
                    plot_start_loc: int, 
                    plot_end_loc: int, 
                    marker_setting: dict = None, 
                    ) -> go.Scatter:
    genome_loc = gwas_summary_bed_table.region_subset(plot_chr, plot_start_loc, plot_end_loc).get_start_locs()
    neg_log10_p = -np.log10(gwas_summary_bed_table.region_subset(plot_chr, plot_start_loc, plot_end_loc).get_region_scores())

    if not marker_setting:
        marker_setting = dict(
            size=5,
            color="black",
        )
    
    return go.Scatter(x=genome_loc, 
                      y=neg_log10_p, 
                      mode="markers", 
                      marker=marker_setting, 
                      name="GWAS summary", 
                      )

def add_annotation_track(fig: go.Figure,
                         annotation_bed_path: str, 
                         annotation_bed_type: str,
                         region_chr: str,
                         region_start: int,
                         region_end: int,
                         row: int, 
                         col: int, 
                         ) -> None:

    match annotation_bed_type:
        case "bed3":
            annotation_bed_table = BedTable3()
            annotation_bed_table.load_from_file(annotation_bed_path)
        case "bed6":
            annotation_bed_table = BedTable6()
            annotation_bed_table.load_from_file(annotation_bed_path)
        case "bed6Plus_with_gene_name":
            annotation_bed_table = BedTable6Plus(extra_column_names=["gene_name"], 
                                                 extra_column_dtype=["str"],)
            annotation_bed_table.load_from_file(annotation_bed_path)
        case _:
            raise ValueError(f"Unrecognized annotation file type: {annotation_bed_type}")

    annotation_bed_table_plot_subset = annotation_bed_table.region_subset(region_chr, 
                                                                          region_start, 
                                                                          region_end, 
                                                                          )

    for region in annotation_bed_table_plot_subset.iter_regions():
        match annotation_bed_type:
            case "bed3":
                region_name = ""
            case "bed6":
                region_name = region["name"]
            case "bed6Plus_with_gene_name":
                region_name = region["gene_name"]

        fig.add_trace(go.Scatter(x=[region["start"], region["end"]],
                                 y=[0, 0], 
                                 mode="lines",
                                 line=dict(width=10, color="blue"),
                                 name=region_name,
                                 ),
                                 row=row,
                                 col=col,
                                 )

def set_parser(parser:argparse.ArgumentParser):
    parser.add_argument("--gwas_summary_bed_path", 
                        help="path to GWAS summary bed file (bed6 format).", 
                        )
    
    parser.add_argument("--annotation_bed_path", 
                        help="path to annotation bed file (bed3 format).", 
                        action="append",
                        )

    parser.add_argument("--annotation_name", 
                        help="Name of annotations. ", 
                        action="append",
                        )
    
    parser.add_argument("--annotation_file_type",
                        help="Type of annotation file. ", 
                        action="append",
                        choices=["bed3", "bed6", "bed6Plus_with_gene_name"],
                        )

    parser.add_argument("--plot_chr", 
                        help="Chromosome to plot.", 
                        )
    
    parser.add_argument("--plot_start_loc", 
                        help="Start location of the plot.", 
                        type=int, 
                        )
    
    parser.add_argument("--plot_end_loc", 
                        help="End location of the plot.", 
                        type=int, 
                        )

def main(args):
    num_annotation_tracks = len(args.annotation_bed_path)

    fig = init_fig(num_annotation_tracks)

    # Plot GWAS manhattan plot
    gwas_summary_bed_table = BedTable6()
    gwas_summary_bed_table.load_from_file(args.gwas_summary_bed_path)
    
    gwas_scatter = make_gwas_scatter(gwas_summary_bed_table, args.plot_chr, args.plot_start_loc, args.plot_end_loc)

    fig.add_trace(gwas_scatter, 
                  row=1, 
                  col=1, 
                  )

    # Plot annotation tracks
    for anno_ind in range(len(args.annotation_bed_path)):
        annotation_bed_path = args.annotation_bed_path[anno_ind]
        annotation_bed_type = args.annotation_file_type[anno_ind]

        add_annotation_track(fig, 
                             annotation_bed_path, 
                             annotation_bed_type,
                             region_chr=args.plot_chr,
                             region_start=args.plot_start_loc,
                             region_end=args.plot_end_loc,
                             row=anno_ind+2, 
                             col=1, 
                             )

    # update layout for Manhattan plot
    fig.update_layout(xaxis=dict(range=[args.plot_start_loc, args.plot_end_loc]))

    # update global layout for the plot
    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
        showlegend=False,  # No legend
        title_text='GWAS Visualization',  # Title
        title_x=0.5,  # Center the title
    )

    # update layout for annotation tracks
    for annot_ind, annot_name in enumerate(args.annotation_names):
        fig.add_annotation(
            text=annot_name, 
            x=args.plot_start_loc,  # Position at the very beginning of the x-axis
            y=0,  # Centered vertically
            xref='paper',  # Use 'paper' to position relative to the entire figure
            yref='paper',
            showarrow=False,  # Do not show an arrow pointing to the text
            align='left',
            xanchor='right',  # Anchor text at the right to keep it outside the plot area
            yanchor='middle', 
            borderpad=10, 
            row=annot_ind+2, 
            col=1, 
        )
        fig.update_layout({
            "yaxis{:d}".format(annot_ind+2):dict(
                ticks='',
                showticklabels=False, 
                range=[-1, 1],
            ), 
        })

    fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="GWAS Visualizer")
    set_parser(parser)

    args = parser.parse_args()

    main(args=args)
