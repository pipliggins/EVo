import argparse

from dgs import run_evo


def main(argv=None):
    my_parser = argparse.ArgumentParser(
        prog="dgs", description="Run EVo: a thermodynamic magma degassing model"
    )

    # Add the arguments
    my_parser.add_argument("chem", metavar="chem.yaml", help="the magma chemistry file")

    my_parser.add_argument(
        "env", metavar="env.yaml", help="the run environment settings file"
    )

    my_parser.add_argument(
        "--output-options",
        help="use selected output options from output.yaml file",
    )

    my_parser.add_argument(
        "-o", "--output", help="the folder location to write the results to"
    )

    # Parse in files
    args = my_parser.parse_args(argv)

    f_chem = args.chem  # set chemical compositions file
    f_env = args.env  # set environment file

    if args.output_options:
        f_out = args.output_options  # set output file as an optional input
        run_evo(f_chem, f_env, f_out, folder=args.output)
    else:
        f_out = None
        run_evo(f_chem, f_env, f_out, folder=args.output)


if __name__ == "__main__":
    main()
