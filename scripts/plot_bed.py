#!/usr/bin/env python

# Plot figures from bed files

import argparse

import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from RGTools.BedTable import BedTable3, BedTable6

class PlotBed:
    @staticmethod
    def set_parser(parser):
        subparser = parser.add_subparsers(dest="subcommand")

        plot_size_distribution_parser = subparser.add_parser("plot_size_distribution",
                                                             help="Plot size distribution of bed file.",
                                                             )
        
        PlotBed.set_parser_plot_size_distribution(plot_size_distribution_parser)

    @staticmethod
    def set_parser_plot_size_distribution(parser):
        PlotBed.set_parser_bed_io(parser)
        
        parser.add_argument("--opath", "-O",
                            help="Output path for figure.",
                            required=True,
                            dest="opath",
                            )
        
        parser.add_argument("--sample", "-S",
                            help="Sample name.",
                            required=True,
                            dest="sample",
                            )
        
        parser.add_argument("--nbins", "-N",
                            help="Number of bins for histogram.",
                            type=int,
                            default=20,
                            dest="nbins",
                            )
        
        parser.add_argument("--size_max", "-M",
                            help="Maximum size for x-axis.",
                            type=int,
                            default=1000,
                            dest="size_max",
                            )
        
    @staticmethod
    def set_parser_bed_io(parser):
        parser.add_argument("--inpath", "-I", 
                            help="Input path for bed file.", 
                            required=True, 
                            dest="inpath", 
                            )
        
        parser.add_argument("--input_file_type", "-T",
                            help="Input file type.",
                            required=True,
                            choices=PlotBed.get_filetype2class_dict().keys(),
                            dest="input_file_type",
                            )
    
    @staticmethod
    def get_filetype2class_dict():
        return {"bed3": BedTable3, 
                "bed6": BedTable6, 
                }
    
    @staticmethod
    def plot_size_distribution(bt, x_lim=(0, 1000), nbins=20, 
                               title="Size distribution"):
        '''
        Plot size distribution of a BedTable object.
        
        Keyword arguments:
        - bt: BedTable object for regions
        - x_lim: Tuple of two integers for x-axis limit
        - title: Title of the plot

        Return plt.Figure 
        '''
        sizes = []

        for region in bt.iter_regions():
            size = region["end"] - region["start"]
            sizes.append(size)

        sizes = np.array(sizes)

        density = gaussian_kde(sizes)
        x_vals = np.linspace(x_lim[0], x_lim[1], 100)
        cumulative_density = np.cumsum(density(x_vals)) * (x_vals[1] - x_vals[0]) 

        fig, axes = plt.subplots(1, 2, figsize=(10,6))

        axes[0].hist(sizes, 
                     bins=nbins, 
                     density=True,
                     alpha=0.3, 
                     color="gray", 
                     )
        axes[0].plot(x_vals, density(x_vals), color='k', label='Density')
        axes[0].set_xlim(x_lim)
        axes[0].set_title("Size distribution")

        axes[1].hist(sizes,
                     bins=nbins, 
                     cumulative=True,
                     density=True,
                     alpha=0.3, 
                     color="gray", 
                     )
        axes[1].plot(x_vals, cumulative_density, color='k', label='Density')
        axes[1].set_xlim(x_lim)
        axes[1].set_title("Cumulative size distribution")

        fig.suptitle(title)

        fig.tight_layout()
    
        return fig

    @staticmethod
    def plot_size_distribution_main(args):
        bt = PlotBed.get_filetype2class_dict()[args.input_file_type]()
        bt.load_from_file(args.inpath)

        fig = PlotBed.plot_size_distribution(bt, 
                                             x_lim=(0, args.size_max), 
                                             nbins=args.nbins, 
                                             title="Size distribution of {}".format(args.sample),
                                             )

        fig.savefig(args.opath)

    @staticmethod
    def main(args):
        if args.subcommand == "plot_size_distribution":
            PlotBed.plot_size_distribution_main(args)
        else:
            raise ValueError("Invalid subcommand.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot figures from bed files")
    PlotBed.set_parser(parser)
    args = parser.parse_args()
    PlotBed.main(args)
