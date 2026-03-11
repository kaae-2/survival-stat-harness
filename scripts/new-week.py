#!/usr/bin/env python3

import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 2 or not sys.argv[1].strip():
        raise SystemExit("Usage: python scripts/new-week.py <week-number>")

    week = sys.argv[1].strip()
    week_dir = Path(f"r-survival-stats-week{week}")
    week_dir.mkdir(parents=True, exist_ok=True)

    (week_dir / f"practical-{week}.R").write_text(
        """# Practical exercises for this week
# TODO: add your practical tasks here.
""",
        encoding="utf-8",
    )

    (week_dir / f"prep-{week}.R").write_text(
        """# Prep exercises for this week
# TODO: add your prep tasks here.
""",
        encoding="utf-8",
    )

    (week_dir / f"r-survival-stats-week{week}.Rproj").write_text(
        """Version: 1.0

RestoreWorkspace: Default
SaveWorkspace: Default
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX
""",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
