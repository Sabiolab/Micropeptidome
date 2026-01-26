#!/usr/bin/python3

import os
import sys
import argparse
from src.pipeline import Pipeline


class Annotator:
    def __init__(self):
        self.args = self.__get_args()

    def __get_args(self):
        self.main_parser = argparse.ArgumentParser(
            description="Annotator",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        self.mode_parser = self.main_parser.add_argument_group("Mode input options")
        self.mode_parser.add_argument(
            "mode",
            metavar="Mode",
            help="Mode to run the pipeline for.\nList of Modes: smorf_types"
        )

        # Parse the first set of arguments to get the mode
        self.args = self.main_parser.parse_args(sys.argv[1:2])
        self.mode = self.args.mode

        # Check for supported modes
        supported_modes = ["smorf_types"]
        if self.mode not in supported_modes:
            self.main_parser.error(f"Unsupported mode '{self.mode}'. Supported modes: {', '.join(supported_modes)}")

        # Create a new parser for mode-specific arguments
        self.parser = argparse.ArgumentParser(
            description=f"Run pipeline in {self.mode} mode",
            prog=" ".join(sys.argv[0:2]),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        self.general_args = self.parser.add_argument_group("General Parameters")
        self.general_args.add_argument("mode", metavar=self.mode)
        self.general_args.add_argument("--outdir", "-o", help="Inform the output directory", default="Annotator_output")
        self.general_args.add_argument("--threads", "-p", help="Number of threads to be used.", type=int, default=1)

        # Add mode-specific arguments
        self.__configure_mode()

        # Parse the full set of arguments
        args = self.parser.parse_args()

        # Ensure output directory exists
        os.makedirs(args.outdir, exist_ok=True)

        return args

    def __configure_mode(self):
        if self.mode == 'smorf_types':
            self.__set_annotator_mode()

    def __set_annotator_mode(self):
        self.modeArguments = self.parser.add_argument_group("Training mode options")
        self.modeArguments.add_argument("--smorf_gtf", help="Provide the smORF GTF file", required=True)
        self.modeArguments.add_argument("--ensembl_gtf", help="Provide the ENSEMBL GTF file", required=True)
        self.modeArguments.add_argument("--intersect_output", help="Provide the intersect output file", default=f"{self.general_args.get_default('outdir')}/intersect.gtf")
        self.modeArguments.add_argument("--non_intersect_output", help="Provide the non-intersect output file", default=f"{self.general_args.get_default('outdir')}/nonintersect.gtf")
        self.modeArguments.add_argument("--output_file", help="Provide the output file", default=f"{self.general_args.get_default('outdir')}/smORF_annotation.txt")


    def execute(self):
        if self.mode == 'smorf_types':
            pipeline = Pipeline(args=self.args)
            pipeline.annotate()

if __name__ == '__main__':
    print("""
 ▗▄▖ ▗▖  ▗▖▗▖  ▗▖ ▗▄▖▗▄▄▄▖▗▄▖▗▄▄▄▖▗▄▖ ▗▄▄▖ 
▐▌ ▐▌▐▛▚▖▐▌▐▛▚▖▐▌▐▌ ▐▌ █ ▐▌ ▐▌ █ ▐▌ ▐▌▐▌ ▐▌
▐▛▀▜▌▐▌ ▝▜▌▐▌ ▝▜▌▐▌ ▐▌ █ ▐▛▀▜▌ █ ▐▌ ▐▌▐▛▀▚▖
▐▌ ▐▌▐▌  ▐▌▐▌  ▐▌▝▚▄▞▘ █ ▐▌ ▐▌ █ ▝▚▄▞▘▐▌ ▐▌
                                           
                                           
                                           """)

    print("""Annotator is a platform that annotates smORF types.""")
    Annotator = Annotator()
    Annotator.execute()
