import argparse
from collections import defaultdict
from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


POPULATIONS = ("Milky Way", "M31")
POPULATION_RE = re.compile(
    r"^(Milky Way|M31)\s*:\s*.*\bN_exp=([+\-0-9.eE]+)"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot expected Milky Way and M31 lensing events, stacked by "
            "stellar-radius bin."
        )
    )
    parser.add_argument("target_filter", help="Filter to plot, for example F150W.")
    parser.add_argument(
        "input_file",
        nargs="?",
        default="all_monte_carlo.txt",
        help="Monte Carlo results file (default: all_monte_carlo.txt).",
    )
    parser.add_argument(
        "--output",
        default="hist_by_rad.pdf",
        help=(
            "Two-panel output path; _MW and _M31 are added for the single-panel "
            "figures (default: hist_by_rad.pdf)."
        ),
    )
    return parser.parse_args()


def read_expected_events(input_file, target_filter):
    # values[population][(log10_mass, median_log_R)] = summed N_exp
    values = {population: defaultdict(float) for population in POPULATIONS}

    current_median_log_R = None
    current_filt_name = None
    current_log10_mass = None

    with input_file.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("median_log_R"):
                current_median_log_R = float(line.split(maxsplit=1)[1])
                current_filt_name = None
                current_log10_mass = None
                continue

            if line.startswith("filt_name"):
                current_filt_name = line.split(maxsplit=1)[1]
                continue

            if line.startswith("log10_mass"):
                current_log10_mass = float(line.split(maxsplit=1)[1])
                continue

            match = POPULATION_RE.match(line)
            if match is None:
                continue

            if (
                current_median_log_R is None
                or current_filt_name is None
                or current_log10_mass is None
            ):
                raise ValueError(
                    "Incomplete metadata before population result on line "
                    f"{line_number}"
                )

            if current_filt_name == target_filter:
                population = match.group(1)
                n_exp = float(match.group(2))
                key = (current_log10_mass, current_median_log_R)
                values[population][key] += n_exp

    if not any(values[population] for population in POPULATIONS):
        raise ValueError(f"No entries found for filt_name = {target_filter}")

    for population in POPULATIONS:
        if not values[population]:
            raise ValueError(
                f"No {population} entries found for filt_name = {target_filter}"
            )

    return values


def plot_expected_events(values, target_filter, populations, output_file):
    all_keys = [
        key
        for population in populations
        for key in values[population]
    ]
    mass_values = np.array(sorted({mass for mass, radius in all_keys}))
    median_log_R_values = np.array(
        sorted({radius for mass, radius in all_keys})
    )

    if len(mass_values) > 1:
        bar_width = 0.9 * np.min(np.diff(mass_values))
    else:
        bar_width = 0.09

    cmap = plt.colormaps["gist_rainbow"]
    if len(median_log_R_values) == 1:
        radius_padding = 0.5
    else:
        radius_padding = 0.0
    norm = mpl.colors.Normalize(
        vmin=median_log_R_values.min() - radius_padding,
        vmax=median_log_R_values.max() + radius_padding,
    )

    n_panels = len(populations)
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=(6, 4) if n_panels == 1 else (10, 4),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )
    axes = axes[0]

    for axis, population in zip(axes, populations):
        bottom = np.zeros(len(mass_values))

        for median_log_R in median_log_R_values:
            n_exp_values = np.array(
                [
                    values[population].get((log10_mass, median_log_R), 0.0)
                    for log10_mass in mass_values
                ]
            )
            axis.bar(
                mass_values,
                n_exp_values,
                width=bar_width,
                bottom=bottom,
                color=cmap(norm(median_log_R)),
                edgecolor="none",
            )
            bottom += n_exp_values

        axis.set_title(population)
        axis.set_xlabel(r"$\log_{10}(M_{\mathrm{PBH}}/M_{\odot})$")
        axis.set_ylim(bottom=0.0)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

    axes[0].set_ylabel("Expected Number of Detectable Lensing Events")
    fig.suptitle(target_filter)

    scalar_mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(
        scalar_mappable,
        ax=axes,
        pad=0.02,
        fraction=0.04,
    )
    colorbar.set_label(r"Binned $\log_{10}(R/R_{\odot})$")

    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def single_population_output_path(output_file, suffix):
    extension = output_file.suffix or ".pdf"
    stem = output_file.stem if output_file.suffix else output_file.name
    return output_file.with_name(f"{stem}_{suffix}{extension}")


def main():
    args = parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output)

    values = read_expected_events(input_file, args.target_filter)
    output_files = (
        (
            ("Milky Way",),
            single_population_output_path(output_file, "MW"),
        ),
        (
            ("M31",),
            single_population_output_path(output_file, "M31"),
        ),
        (POPULATIONS, output_file),
    )

    for populations, this_output_file in output_files:
        plot_expected_events(
            values,
            args.target_filter,
            populations,
            this_output_file,
        )
        print(f"Wrote {this_output_file}")


if __name__ == "__main__":
    main()
