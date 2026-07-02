#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ROOT = Path("project_primary_FASTA")
F2DIR = ROOT / "12_suppfig2_read_support"
REFDIR = ROOT / "10_chr14_mapback/01_refs"
FLAGSTAT = F2DIR / "02_tables/raw_snapshot/flagstat_summary.tsv"
OUTDIR = F2DIR / "04_charts"
TABLEDIR = F2DIR / "02_tables/figure_data"

OUTDIR.mkdir(parents=True, exist_ok=True)
TABLEDIR.mkdir(parents=True, exist_ok=True)


MODELS = [
    {
        "row": 0,
        "col": 0,
        "test": "test01_direct_noGap",
        "row_label": "Direct",
        "column_label": "No inserted gap",
        "adopted": False,
    },
    {
        "row": 0,
        "col": 1,
        "test": "test05_direct_gap10k",
        "row_label": "Direct",
        "column_label": "10-kb N-spacer model",
        "adopted": False,
    },
    {
        "row": 1,
        "col": 0,
        "test": "test02_B_noGap",
        "row_label": "+ ptg000010l",
        "column_label": "No inserted gap",
        "adopted": False,
    },
    {
        "row": 1,
        "col": 1,
        "test": "test06_B_gap10k",
        "row_label": "+ ptg000010l",
        "column_label": "10-kb N-spacer model",
        "adopted": False,
    },
    {
        "row": 2,
        "col": 0,
        "test": "test03_G_noGap",
        "row_label": "+ ptg000020l",
        "column_label": "No inserted gap",
        "adopted": True,
    },
    {
        "row": 2,
        "col": 1,
        "test": "test07_G_gap10k",
        "row_label": "+ ptg000020l",
        "column_label": "10-kb N-spacer model",
        "adopted": False,
    },
    {
        "row": 3,
        "col": 0,
        "test": "test04_BG_noGap",
        "row_label": "+ ptg000010l and ptg000020l",
        "column_label": "No inserted gap",
        "adopted": False,
    },
    {
        "row": 3,
        "col": 1,
        "test": "test08_BG_gap10k",
        "row_label": "+ ptg000010l and ptg000020l",
        "column_label": "10-kb N-spacer model",
        "adopted": False,
    },
]


COMPONENT_STYLE = {
    "Chr14_A1": {
        "display": "ptg000025l",
        "color": "#0072B2",
        "width": 1.00,
    },
    "Chr14_A2": {
        "display": "ptg000027l",
        "color": "#56B4E9",
        "width": 1.00,
    },
    "Chr14_B": {
        "display": "ptg000010l",
        "color": "#A6A6A6",
        "width": 0.82,
    },
    "Chr14_E": {
        "display": "ptg000015l\nsegment",
        "color": "#009E73",
        "width": 1.10,
    },
    "Chr14_G": {
        "display": "ptg000020l",
        "color": "#E69F00",
        "width": 0.82,
    },
    "Chr14_H1": {
        "display": "ptg000003l",
        "color": "#CC79A7",
        "width": 1.00,
    },
    "Chr14_H3": {
        "display": "ptg000004l\nsegment",
        "color": "#D55E00",
        "width": 1.10,
    },
    "GAP": {
        "display": "10 kb\nN",
        "color": "#FFFFFF",
        "width": 0.46,
    },
}


def read_manifest(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    for row in rows:
        row["start_1based"] = int(row["start_1based"])
        row["end_1based"] = int(row["end_1based"])

    return rows


def parse_flagstat(path):
    """
    Raw file columns:
    test, total, mapped, secondary, primary_mapped, primary, supplementary

    Validate that:
    total = secondary + primary + supplementary
    """
    data = {}

    with path.open() as handle:
        for line_no, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue

            fields = raw.rstrip("\n").split("\t")

            if len(fields) != 7:
                raise ValueError(
                    f"Expected 7 columns at line {line_no}, "
                    f"found {len(fields)}"
                )

            test = fields[0]
            values = list(map(int, fields[1:]))

            total, mapped, secondary, primary_mapped, primary, supplementary = values

            if total != secondary + primary + supplementary:
                raise ValueError(
                    f"Flagstat decomposition failed for {test}: "
                    f"{total} != {secondary} + {primary} + {supplementary}"
                )

            data[test] = {
                "total": total,
                "mapped": mapped,
                "secondary": secondary,
                "primary_mapped": primary_mapped,
                "primary": primary,
                "supplementary": supplementary,
            }

    return data


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def text_color(hex_color):
    dark_colors = {
        "#0072B2",
        "#009E73",
        "#CC79A7",
        "#D55E00",
    }
    return "white" if hex_color in dark_colors else "black"


def draw_model(ax, manifest, adopted, metrics=None):
    ax.set_axis_off()

    left_margin = 0.12
    right_margin = 0.12
    box_y = 0.43
    box_h = 0.30
    gap_between = 0.035

    widths = [
        COMPONENT_STYLE[row["component"]]["width"]
        for row in manifest
    ]

    total_width_units = (
        sum(widths)
        + gap_between * max(len(widths) - 1, 0)
    )

    available_width = 1.0 - left_margin - right_margin
    scale = available_width / total_width_units

    x = left_margin

    for index, row in enumerate(manifest):
        component = row["component"]
        style = COMPONENT_STYLE[component]
        width = style["width"] * scale

        if component == "GAP":
            rectangle = Rectangle(
                (x, box_y),
                width,
                box_h,
                facecolor="white",
                edgecolor="#666666",
                linewidth=0.8,
                hatch="////",
                zorder=2,
            )
        else:
            rectangle = Rectangle(
                (x, box_y),
                width,
                box_h,
                facecolor=style["color"],
                edgecolor="white",
                linewidth=0.8,
                zorder=2,
            )

        ax.add_patch(rectangle)

        ax.text(
            x + width / 2,
            box_y + box_h / 2,
            style["display"],
            ha="center",
            va="center",
            fontsize=6.1 if component != "GAP" else 5.7,
            color=text_color(style["color"]),
            linespacing=0.92,
            zorder=3,
        )

        x += width

        if index < len(manifest) - 1:
            ax.plot(
                [x, x + gap_between * scale],
                [box_y + box_h / 2] * 2,
                color="#555555",
                linewidth=0.7,
                zorder=1,
            )
            x += gap_between * scale

    if adopted:
        ax.add_patch(
            Rectangle(
                (0.015, 0.08),
                0.97,
                0.84,
                fill=False,
                edgecolor="#0072B2",
                linewidth=2.0,
                zorder=5,
            )
        )

        ax.text(
            0.98,
            0.91,
            "adopted current-best model",
            ha="right",
            va="top",
            fontsize=6.8,
            fontweight="bold",
            color="#0072B2",
        )

    if metrics is not None:
        ax.text(
            0.5,
            0.17,
            (
                f"primary mapped: {metrics['primary_mapped']:,}   |   "
                f"supplementary: {metrics['supplementary']:,}"
            ),
            ha="center",
            va="center",
            fontsize=5.8,
            color="#555555",
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


flagstat = parse_flagstat(FLAGSTAT)

model_rows = []
manifest_by_test = {}

for model in MODELS:
    test = model["test"]
    manifest_path = REFDIR / f"{test}.manifest.tsv"

    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)

    manifest = read_manifest(manifest_path)
    manifest_by_test[test] = manifest

    model_rows.append({
        "test": test,
        "row_label": model["row_label"],
        "column_label": model["column_label"],
        "adopted": model["adopted"],
        "component_order": ",".join(
            row["component"] for row in manifest
        ),
        "candidate_length": manifest[-1]["end_1based"],
        "primary_mapped": flagstat[test]["primary_mapped"],
        "supplementary": flagstat[test]["supplementary"],
    })

write_tsv(
    TABLEDIR / "SuppFig2A_model_definitions.tsv",
    model_rows,
    [
        "test",
        "row_label",
        "column_label",
        "adopted",
        "component_order",
        "candidate_length",
        "primary_mapped",
        "supplementary",
    ],
)


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def create_figure(show_metrics, suffix):
    fig, axes = plt.subplots(
        4,
        2,
        figsize=(11.8, 5.5 if show_metrics else 5.0),
    )

    column_titles = [
        "No inserted gap",
        "10-kb N-spacer model",
    ]

    for col, title in enumerate(column_titles):
        axes[0, col].set_title(
            title,
            fontsize=10,
            fontweight="bold",
            pad=8,
        )

    row_labels = [
        "Direct",
        "+ ptg000010l",
        "+ ptg000020l",
        "+ ptg000010l and ptg000020l",
    ]

    for row, label in enumerate(row_labels):
        axes[row, 0].text(
            -0.02,
            0.50,
            label,
            transform=axes[row, 0].transAxes,
            ha="right",
            va="center",
            fontsize=8,
            fontweight="bold",
            clip_on=False,
        )

    for model in MODELS:
        metrics = (
            flagstat[model["test"]]
            if show_metrics
            else None
        )

        draw_model(
            axes[model["row"], model["col"]],
            manifest_by_test[model["test"]],
            model["adopted"],
            metrics=metrics,
        )

    fig.text(
        0.5,
        0.025,
        (
            "Schematics are not to scale. "
            "Source-contig-derived segments are shown using shortened IDs."
        ),
        ha="center",
        va="bottom",
        fontsize=7,
        color="#555555",
    )

    fig.subplots_adjust(
        left=0.19,
        right=0.985,
        top=0.90,
        bottom=0.10,
        hspace=0.24,
        wspace=0.10,
    )

    for extension, kwargs in [
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": 600}),
    ]:
        fig.savefig(
            OUTDIR / f"SuppFig2A_alternative_models_{suffix}.{extension}",
            bbox_inches="tight",
            **kwargs,
        )

    plt.close(fig)


create_figure(
    show_metrics=False,
    suffix="clean",
)

create_figure(
    show_metrics=True,
    suffix="with_flagstat",
)

print("Generated:")
print(OUTDIR / "SuppFig2A_alternative_models_clean.pdf")
print(OUTDIR / "SuppFig2A_alternative_models_clean.svg")
print(OUTDIR / "SuppFig2A_alternative_models_clean.png")
print(OUTDIR / "SuppFig2A_alternative_models_with_flagstat.pdf")
print(OUTDIR / "SuppFig2A_alternative_models_with_flagstat.svg")
print(OUTDIR / "SuppFig2A_alternative_models_with_flagstat.png")
print(TABLEDIR / "SuppFig2A_model_definitions.tsv")
