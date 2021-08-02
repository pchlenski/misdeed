#!python3
"""
MiSDEED command line interface.
"""

import argparse
import logging
import os
import os.path
import sys
from uuid import uuid4

from src import __version__
from src.OmicsGenerator import OmicsGenerator

class options:
	def __init__(self):
# 		parser = argparse.ArgumentParser(
# 			description=f"MiSDEED (v{__version__}): microbiome synthetic data engine for experimental design",
# 			usage="""misdeed [options] output_directory
# mode:	single
# 	multiple
# 	casecontrol
# 	infer
# """
# 		)

# 		# help dialogue
# 		if len(sys.argv[1:]) < 1:
# 			parser.print_help()
# 			sys.exit(2)

# 		# interpret command
# 		# parser.add_argument("command", type=str, help="What command to run.")
# 		args = parser.parse_args(sys.argv[1:2])
# 		getattr(self, args.command)()
# 		# TODO: more helpful error message when command is wrong

	# def single(self):
		parser = argparse.ArgumentParser(
			usage = "misdeed [options] node_sizes output_directory"
		)

		# Generator params
		parser.add_argument(
			"node_sizes",
			help="Number of dimensions for each node, separated by commas"
		)
		parser.add_argument(
			"output_path",
			help="Directory in which to store synthetic data"
		)
		parser.add_argument(
			"--node_names",
			default=None,
			help="Names for each node, separated by commas. If set, must match node_sizes in length."
		)
		parser.add_argument(
			"--time_points", "-t",
			type=int,
			default=100,
			help="Number of time points to generate."
		)
		parser.add_argument(
			"--discard_first",
			type=int,
			default=0,
			help="Number of time points to generate and discard at initialization."
		)

		# Simulation params
		parser.add_argument(
			"--n_samples", "-n",
			type=int,
			default=1,
			help="Number of samples to generate."
		)
		parser.add_argument(
			"--extinct_fraction",
			type=float,
			default=0,
			help="Fraction of clades to set to zero when generating samples."
		)
		parser.add_argument(
			"--noise_variance", "-v",
			type=float,
			default=1e-2,
			help="Variance of biological noise added to system."
		)
		parser.add_argument(
			"--n_reads", "-r",
			type=int,
			default=1e5,
			help="Number of reads drawn per sample."
		)
		parser.add_argument(
			"--time_step",
			type=float,
			default=1e-2,
			help="Size of time increment between samples."
		)
		parser.add_argument(
			"--downsample", "-d",
			type=int,
			default=1,
			help="Downsampling coefficient (keep every d-th sample)."
		)

		# Case-control params
		parser.add_argument(
			"--case_fraction", "-c",
			type=float,
			default=0,
			help="Fraction of samples to generate as cases (default = 0)"
		)
		parser.add_argument(
			"--intervention_node",
			type=str,
			default=None,
			help="Node to which intervention is applied (default = first node)"
		)
		parser.add_argument(
			"--effect_size",
			type=float,
			default=0.1
		)

		# if len(sys.argv[2:]) < 2:
		if len(sys.argv[1:]) < 2:
			parser.print_help()
			sys.exit(2)
		# TODO: add arg for matrix

		# Process args
		# args = parser.parse_args(sys.argv[2:])
		args = parser.parse_args(sys.argv[1:])
		
		node_sizes = [int(x) for x in args.node_sizes.split(",")]

		if args.node_names == None:
			node_names = []
		else:
			node_names = args.node_names.split(",")

		# Path checking
		# if os.path.isdir(output_path)

		# Set up generator
		print("==========================SETUP========================")
		print()
		gen = OmicsGenerator(
			nodes=node_names,
			node_sizes=node_sizes,
			time_points=args.time_points,
			discard_first=args.discard_first,
			init_full=True
		)
		print(gen)

		print("========================SIMULATION=====================")
		print()
		# Run generator
		if args.n_samples == 1:
			print(f"Generating single sample...")
			x,y,z = gen.generate(
				noise_var=args.noise_variance,
				n_reads=args.n_reads,
				dt=args.time_step,
				downsample=args.downsample
			)
		elif args.case_fraction == 0:
			print(f"Generating {args.n_samples} samples...")
			x,y,z = gen.generate_multiple(
				args.n_samples,
				extinct_fraction=args.extinct_fraction,
				noise_var=args.noise_variance,
				n_reads=args.n_reads,
				dt=args.time_step,
				downsample=args.downsample
			)
		else:
			print(f"Generating {args.n_samples} samples (case-control)...")

			if args.intervention_node == None:
				args.intervention_node = gen.nodes[0].name 

			x_control,y_control,z_control,x_case,y_case,z_case = gen.case_control(
				args.n_samples,
				case_frac=args.case_fraction,
				node_name=args.intervention_node,
				effect_size=args.effect_size,
				extinct_fraction=args.extinct_fraction,
				noise_var=args.noise_variance,
				n_reads=args.n_reads,
				dt=args.time_step,
				downsample=args.downsample
			)

		# Save outputs
		print("=========================SAVING========================")
		print()
		path = f"{args.output_path}/{uuid4()}"
		os.mkdir(path)
		if args.case_fraction == 0:
			gen.save(x, f"{path}/X")
			gen.save(y, f"{path}/Y")
			gen.save(z, f"{path}/Z")
		else:
			os.mkdir(f"{path}/control")
			os.mkdir(f"{path}/case")
			gen.save(x_control, f"{path}/control/X")
			gen.save(x_case, f"{path}/case/X")
			gen.save(y_control, f"{path}/control/Y")
			gen.save(y_case, f"{path}/case/Y")
			gen.save(z_control, f"{path}/control/Z")
			gen.save(z_case, f"{path}/case/Z")

options()