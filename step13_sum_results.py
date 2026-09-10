import argparse
from collections import defaultdict
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

from FileRead import writecol


POPULATIONS = ("Milky Way", "M31")
POPULATION_RE = re.compile(
    r"^(Milky Way|M31)\s*:\s*.*\bN_exp=([+\-0-9.eE]+)"
)
FILTER_COLORS = {
    "F150W": "tab:blue",
    "F277W": "tab:red",
    "r": "black",
    "H": "tab:purple",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Sum star-weighted expected event counts and plot PBH constraints "
            "for the Milky Way and M31 lens populations."
        )
    )
    parser.add_argument(
        "scale_by",
        type=float,
        help="Additional multiplicative scale applied after star weighting.",
    )
    parser.add_argument("input_files", nargs="+", help="Monte Carlo result files.")
    parser.add_argument(
        "--output",
        default="sum_results.pdf",
        help=(
            "Two-population output path; _MW and _M31 are added for the "
            "single-population figures (default: sum_results.pdf)."
        ),
    )
    return parser.parse_args()


def read_counts(input_files, scale_by):
    # counts[population][(log10_mass, filter)] = star-weighted expected events
    counts = {population: defaultdict(float) for population in POPULATIONS}

    for input_file in input_files:
        current_number_of_stars = None
        current_filter = None
        current_log10_mass = None

        with input_file.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith("median_log_R"):
                    current_number_of_stars = None
                    current_filter = None
                    current_log10_mass = None
                    continue

                if line.startswith("number_of_stars"):
                    current_number_of_stars = float(line.split(maxsplit=1)[1])
                    if current_number_of_stars < 0.0:
                        raise ValueError(
                            f"Negative number_of_stars in {input_file} "
                            f"on line {line_number}"
                        )
                    continue

                if line.startswith("filt_name"):
                    current_filter = line.split(maxsplit=1)[1]
                    continue

                if line.startswith("log10_mass"):
                    current_log10_mass = float(line.split(maxsplit=1)[1])
                    continue

                match = POPULATION_RE.match(line)
                if match is None:
                    continue

                if (
                    current_number_of_stars is None
                    or current_filter is None
                    or current_log10_mass is None
                ):
                    raise ValueError(
                        f"Incomplete metadata before population result in "
                        f"{input_file} on line {line_number}"
                    )

                population = match.group(1)
                n_exp_per_star = float(match.group(2))
                key = (current_log10_mass, current_filter)
                counts[population][key] += (
                    n_exp_per_star * current_number_of_stars * scale_by
                )

    for population in POPULATIONS:
        if not counts[population]:
            raise ValueError(f"No {population} results found in the input files")

    return counts


def filter_color(filter_name, filter_index):
    return FILTER_COLORS.get(filter_name, f"C{filter_index % 10}")


def population_series(counts, population):
    filters = sorted({filter_name for mass, filter_name in counts[population]})
    series = {}

    for filter_name in filters:
        masses = np.array(
            sorted(
                mass
                for mass, this_filter in counts[population]
                if this_filter == filter_name
            )
        )
        n_exp = np.array(
            [counts[population][(mass, filter_name)] for mass in masses]
        )
        series[filter_name] = (10.0**masses, n_exp)

    return series


def add_five_percent_headroom(axis, zero_bottom=False):
    default_ymin, default_ymax = axis.get_ylim()
    if zero_bottom:
        default_ymin = 0.0
    axis.set_ylim(default_ymin, 1.05 * default_ymax)


def plot_results(counts, populations, output_file):
    n_columns = len(populations)
    fig, axes = plt.subplots(
        2,
        n_columns,
        figsize=(6, 7) if n_columns == 1 else (10, 7),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )

    for column, population in enumerate(populations):
        event_axis = axes[0, column]
        constraint_axis = axes[1, column]
        series = population_series(counts, population)

        for filter_index, (filter_name, (masses, n_exp)) in enumerate(
            series.items()
        ):
            positive = n_exp > 0.0
            color = filter_color(filter_name, filter_index)
            event_axis.plot(
                masses,
                n_exp,
                marker=".",
                linestyle="none",
                color=color,
                label=filter_name,
            )
            constraint_axis.plot(
                masses[positive],
                3.0 / n_exp[positive],
                marker=".",
                linestyle="none",
                color=color,
                label=filter_name,
            )

        event_axis.set_title(population)
        event_axis.set_xscale("log")
        event_axis.set_xlim(1.0e-11, 1.0e-5)
        add_five_percent_headroom(event_axis, zero_bottom=True)
        event_axis.legend()

        constraint_axis.set_xscale("log")
        constraint_axis.set_yscale("log")
        constraint_axis.set_xlim(1.0e-11, 1.0e-5)
        add_five_percent_headroom(constraint_axis)
        constraint_axis.set_xlabel(r"$M_{\mathrm{PBH}}/M_{\odot}$")

    axes[0, 0].set_ylabel("Expected Number of Detectable Lensing Events")
    axes[1, 0].set_ylabel(r"Constraint $3/N_{\mathrm{exp}}$")

    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def single_population_output_path(output_file, suffix):
    extension = output_file.suffix or ".pdf"
    stem = output_file.stem if output_file.suffix else output_file.name
    return output_file.with_name(f"{stem}_{suffix}{extension}")


def write_f150w_constraints(counts, population, output_file):
    series = population_series(counts, population)
    if "F150W" not in series:
        return

    masses, n_exp = series["F150W"]
    positive = n_exp > 0.0
    writecol(
        str(output_file),
        [masses[positive], 3.0 / n_exp[positive]],
    )
    print(f"Wrote {output_file}")


def main():
    args = parse_args()
    input_files = [Path(path) for path in args.input_files]
    output_file = Path(args.output)
    counts = read_counts(input_files, args.scale_by)

    output_files = (
        (("Milky Way",), single_population_output_path(output_file, "MW")),
        (("M31",), single_population_output_path(output_file, "M31")),
        (POPULATIONS, output_file),
    )
    for populations, this_output_file in output_files:
        plot_results(counts, populations, this_output_file)
        print(f"Wrote {this_output_file}")

    scale_label = str(args.scale_by)
    write_f150w_constraints(
        counts,
        "Milky Way",
        Path(f"constraint_sum_F150W_{scale_label}x.txt"),
    )
    write_f150w_constraints(
        counts,
        "M31",
        Path(f"constraint_sum_F150W_M31_{scale_label}x.txt"),
    )


if __name__ == "__main__":
    main()
