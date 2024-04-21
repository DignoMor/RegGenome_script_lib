#!/usr/bin/env python

import argparse

import plotly.graph_objects as go
import numpy as np
import pandas as pd

from plotly.subplots import make_subplots

from RGTools.BedTable import BedTable3, BedTable6

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

def set_parser(parser:argparse.ArgumentParser):
    parser.add_argument("--gwas_summary_bed_path", 
                        help="path to GWAS summary bed file (bed6 format).", 
                        )
    
    parser.add_argument("--annotation_bed_path", 
                        help="path to annotation bed file (bed3 format).", 
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
    gwas_summary_bed_table = BedTable6()
    gwas_summary_bed_table.load_from_file(args.gwas_summary_bed_path)

    annotation_bed_table = BedTable3()
    annotation_bed_table.load_from_file(args.annotation_bed_path)

    fig = make_subplots(rows=2, cols=1, 
                        row_heights=[0.9, 0.1],
                        shared_xaxes=True,
                        )

    gwas_scatter = make_gwas_scatter(gwas_summary_bed_table, args.plot_chr, args.plot_start_loc, args.plot_end_loc)

    fig.add_trace(gwas_scatter, 
                row=1, 
                col=1, 
                )

    fig.update_layout(xaxis=dict(range=[args.plot_start_loc, args.plot_end_loc]))

    annotation_bed_table_plot_subset = annotation_bed_table.region_subset(args.plot_chr, args.plot_start_loc, args.plot_end_loc)

    for region in annotation_bed_table_plot_subset.iter_regions():
        fig.add_trace(go.Scatter(x=[region["start"], region["end"]],
                                y=[1,1], 
                                mode="lines",
                                line=dict(width=10, color="red"),
                                ),
                                row=2,
                                col=1,
                                )

    fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
        showlegend=False,  # No legend
        xaxis1=dict(
            showline=False,  # Remove axis line
            showgrid=False,  # Remove grid
            zeroline=False,  # Remove zero line
            ticks='',  # Remove ticks
            showticklabels=False  # Remove tick labels
        ),
        yaxis2=dict(
            showline=False,
            showgrid=False,
            zeroline=False,
            ticks='',
            showticklabels=False, 
            range=[0, 2],
        ), 
    )
    
    fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="GWAS Visualizer")
    set_parser(parser)

    args = parser.parse_args()

    main(args=args)
