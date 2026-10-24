
#!/usr/bin/env python
# coding: utf-8

import os
import re
import tempfile
import argparse
from time import time
from pathlib import Path
import sys


print(f"Welcome to the LaTeX Runner!\nPython Version: {sys.version}\n")


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def get_cli_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Compile a LaTeX document using an isolated temporary "
            "working directory and copy the generated PDF to the "
            "specified output directory."
        )
    )

    parser.add_argument(
        "filename",
        nargs="?",
        type=str,
        help="Path to the input LaTeX (.tex) file."
    )

    parser.add_argument(
        "--utils",
        nargs="*",
        default="utils",
        help=(
            "Additional files or directories to copy into the temporary "
            "compilation directory, such as LaTeX utilities, style files, "
            "data files, or figures."
        )
    )

    parser.add_argument(
        "--compiler",
        default="xelatex",
        help=(
            "LaTeX compiler to use. Examples: xelatex, pdflatex, "
            "or lualatex. Default: xelatex."
        )
    )

    parser.add_argument(
        "--outdir",
        default=".",
        help=(
            "Directory where the generated PDF will be copied. "
            "Default: current directory."
        )
    )

    parser.add_argument(
        "--stdin",
        action="store_true",
        help=(
            "Enable interactive input and confirmation prompts. "
            "By default, the runner operates non-interactively."
        )
    )

    parser.add_argument(
        "--logs",
        action="store_true",
        help=(
            "Keep LaTeX compilation log and auxiliary files when "
            "supported by the compilation workflow."
        )
    )

    # parser.add_argument(
    #     "--defaults",
    #     action="store_true",
    #     help=(
    #         "Use default configuration values where applicable."
    #     )
    # )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help=(
            "Display detailed LaTeX compiler output instead of "
            "suppressing it to the runtime log."
        )
    )

    parser.add_argument(
        "--times",
        type=int,
        default=2,
        help=(
            "Number of LaTeX compilation passes. At least two passes "
            "are required. Default: 2."
        )
    )

    return parser


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def file_exists(p: Path) -> bool:
    if not p.exists():
        print(f"{p} does NOT EXIST!")
        return False

    return True


def expand_glob(pattern: str):
    """Support simple * glob via regex fallback"""

    if "*" not in pattern:
        return Path(pattern)

    path_obj = Path(pattern)

    base = (
        path_obj.parent
        if path_obj.parent != Path()
        else Path(".")
    )

    regex_pattern = re.compile(
        pattern.replace("*", ".*")
    )

    matches = [
        base / f
        for f in os.listdir(base)
        if regex_pattern.match(f)
    ]

    return matches if matches else pattern


def flatten(lst):
    out = []

    for x in lst:
        if isinstance(x, list):
            out.extend(flatten(x))
        else:
            out.append(x)

    return out


# ----------------------------------------------------------------------
# Interactive input
# ----------------------------------------------------------------------
def collect_interactive(args, defaults):

    if not args.stdin:
        return defaults.copy()

    cfg = defaults.copy()

    def prompt(text, key, transform=lambda x: x):

        value = input(text).strip()

        if value.lower() == "q":
            return None

        if value.lower() == "d":
            return True

        elif len(value):
            cfg[key] = transform(value)

        return cfg

    # Sources
    result = prompt(
        "Enter source paths (comma separated): ",
        "SRCS",
        lambda srcs: flatten(
            [
                expand_glob(s.strip())
                for s in srcs.split(",")
            ]
        )
        if "," in srcs
        else expand_glob(srcs)
    )

    if not isinstance(result, dict):
        return cfg if result else result

    # Output directory
    result = prompt(
        "Output directory: ",
        "OUTDIR"
    )

    if not isinstance(result, dict):
        return cfg if result else result

    # Verbose
    result = prompt(
        "Print all compilation details? [y/n]: ",
        "VERBOSE",
        lambda x: x.lower().startswith("y")
    )

    if not isinstance(result, dict):
        return cfg if result else result

    # Keep logs
    result = prompt(
        "Keep logs? [y/n]: ",
        "LOGS",
        lambda x: x.lower().startswith("y")
    )

    if not isinstance(result, dict):
        return cfg if result else result

    return cfg


# ----------------------------------------------------------------------
# Resolve filename
# ----------------------------------------------------------------------
def resolve_filename(name: str) -> Path:

    if not name:
        raise ValueError(
            "A LaTeX filename is required."
        )

    p = Path(
        os.path.expandvars(
            os.path.expanduser(name)
        )
    )

    if p.suffix.lower() != ".tex":
        p = p.with_suffix(".tex")

    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found"
        )

    return p


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():

    parser = get_cli_arguments()
    args = parser.parse_args()

    if args.times < 1:
        parser.error("--times must be at least 1")

    host_dir = Path.cwd()

    # --------------------------------------------------------------
    # Resolve filename
    # --------------------------------------------------------------
    try:
        filename = resolve_filename(args.filename)

    except (ValueError, FileNotFoundError) as e:
        parser.error(str(e))

    defaults = {
        "SRCS": args.utils,
        "OUTDIR": args.outdir,
        "VERBOSE": args.verbose,
        "LOGS": args.logs,
    }

    user_cfg = collect_interactive(
        args,
        defaults
    )

    if user_cfg is None:
        return

    SRCS = user_cfg["SRCS"]
    OUTDIR = Path(user_cfg["OUTDIR"])
    VERBOSE = user_cfg["VERBOSE"]
    LOGS = user_cfg["LOGS"]

    # ------------------------------------------------------------------
    # Locate latex helper package
    # ------------------------------------------------------------------
    latex_dir = None

    search_root = Path.cwd()

    target = Path(
        "BI-HUB/MATH-IT/PKGS"
    )

    while search_root != search_root.parent:

        candidate = search_root / target

        if candidate.exists():
            latex_dir = candidate
            break

        search_root = search_root.parent

    # Web execution fallback:
    # Search relative to this script if the normal cwd-based search fails.
    if latex_dir is None:

        runner_dir = Path(__file__).resolve().parent

        candidates = [
            runner_dir / "PKGS",
            runner_dir / "BI-HUB/MATH-IT/PKGS",
        ]

        for candidate in candidates:

            if candidate.exists():
                latex_dir = candidate
                break

    if not latex_dir:

        print(
            "LaTeX helper package not found."
        )

        return

    sys.path.append(str(latex_dir))

    from createLatex import (
        copy_files, 
        run_latex, 
        run_checks
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    summary = [
        "\n" + "=" * 33,
        "\tINPUT SUMMARY",
        "-" * 33,
        f"File: {os.path.relpath(filename)}",
        f"Source files: {SRCS}",
        f"output directory: {os.path.relpath(OUTDIR)}",
        f"Host directory: {host_dir}",
        f"LaTeX compiler: {args.compiler}",
        f"Prints out ALL compilation details? {VERBOSE}",
        f"keep logs files? {LOGS}",
        (
            "Location of Python code used by LaTeX Runner: "
            f"{os.path.relpath(latex_dir)}"
        ),
        "-" * 33 + "\n"
    ]

    for i in range(3, len(summary) - 1):
        summary[i] = f" {i - 2}. {summary[i]}"

    print("\n".join(summary))

    # ------------------------------------------------------------------
    # Interactive pause
    #
    # Only pause when --stdin is explicitly supplied.
    # This allows Flask to run the script noninteractively.
    # ------------------------------------------------------------------
    if args.stdin:
        input("Press ENTER to continue...")

    # ------------------------------------------------------------------
    # Compilation
    # ------------------------------------------------------------------
    with tempfile.TemporaryDirectory(prefix="TeXRun-") as tmp_dir:

        work_dir = Path(tmp_dir)

        print(f"Working dir: {work_dir}")

        t0 = time()
        filename = copy_files(str(filename), SRCS=SRCS, work_dir=str(work_dir))
        print(f"Setup done in {time() - t0:.2f}s")

        try:

            os.chdir(work_dir)

            print(f"Compiling {filename}...")

            t0 = time()

            run_latex(filename, compiler=args.compiler, verbose=VERBOSE, times=args.times)

            print(f"Done in {time() - t0:.2f}s")

        except Exception as e:

            print(f"Compilation error: {e}")

            log = work_dir / "logs/latex.log"

            if log.exists():

                print(log.read_text()[-500:])

            return

        finally:

            os.chdir(host_dir)

        run_checks(
            filename, 
            work_dir=tmp_dir, 
            output_dir=str(OUTDIR), 
            KeepFiles=LOGS, 
            verbose=VERBOSE
        )


if __name__ == "__main__":
    main()