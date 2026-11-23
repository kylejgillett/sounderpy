#!/usr/bin/env python3
"""
Write docs/source/reference/references.rst, the Parameters & References page, from sounderpy_references.py:
a table of every parameter sounding_params(...).calc() returns with its source, and the full reference list.

Run from any directory (paths are found from this script's own location):

    python docs/scripts/make_references_page.py
    python docs/scripts/make_references_page.py --check   # exit 1 if the page is out of date
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sounderpy_references import GROUPS, PARAMETERS, REFS, SHORT, TAGS  # noqa: E402

PAGE = Path(__file__).resolve().parents[1] / "source" / "reference" / "references.rst"

DICTS = (("general", "General values (``general``)",
          "Surface values and values at every level of the profile."),
         ("thermo", "Thermodynamic values (``thermo``)",
          "``sb_``, ``ml_``, and ``mu_`` are the surface-based, mixed-layer (lowest 100 hPa), and most-unstable parcels."),
         ("kinem", "Kinematic values (``kinem``)",
          "Layers are in m AGL (``0_to_500`` is 0-500 m); ``eil`` is the effective inflow layer. Wind vectors are "
          "(u, v) in kt."))

# the winter dictionary's values, by topic (all in ``winter``)
WINTER_TOPICS = (("winter_ptype", "Precipitation type",
                  "From the precipitation source down to the ground. Phases: ``SN`` snow, ``RA`` rain, ``DZ`` drizzle, "
                  "``PL`` ice pellets, ``FZ`` freezing."),
                 ("winter_snow", "Clouds and snow growth", "RH over ice is used below 0 °C."),
                 ("winter_hazards", "Surface and winter hazards",
                  "The freezing fog and drizzle potentials are ingredient checks, not forecasts: ``HIGH`` when every "
                  "check is met, ``MOD`` when one is short, ``LOW`` otherwise."),
                 ("winter_winds", "Winds", "Wind vectors are (u, v) in kt."),
                 ("winter_profiles", "Profiles", "For drawing your own panels."))


def rst_text(text):
    """A reference's text in rst: DOIs with < or > become explicit links (rst can't auto-link them)."""
    def link(match):
        url = match.group(0)
        if "<" not in url and ">" not in url:
            return url
        shown = url.replace("<", r"\<").replace(">", r"\>")
        return f"`{shown} <{url.replace('<', '%3C').replace('>', '%3E')}>`__"
    return re.sub(r"https?://\S+?(?=,?\s|$)", link, text)


def heading(text, char):
    return [text, char * len(text), ""]


def table(name, note):
    """A parameter table: every PARAMETERS row for `name`, after a one-line note."""
    lines = [note, "", ".. list-table::", "   :header-rows: 1", "   :widths: 26 44 12 18", "   :class: param-table", "",
             "   * - Key(s)", "     - Parameter", "     - Units", "     - Source"]
    for dict_name, keys, parameter, units, refs in PARAMETERS:
        if dict_name != name:
            continue
        source = ", ".join(f":ref:`{SHORT[k]} <ref-{k}>`" for k in refs) or "SounderPy"
        lines += [f"   * - ``{keys}``", f"     - {parameter}", f"     - {units}", f"     - {source}"]
    return lines + [""]


def build():
    lines = [".. _references:", "", "=" * 30, "📖 Parameters & References", "=" * 30, "",
             "The parameters SounderPy computes, where each one comes from, and the full list of references behind",
             "SounderPy's calculations, data sources, and plots. The same references are listed in the header at the",
             "top of each SounderPy module.", "",
             "*" * 63, ""]

    lines += heading("Citing SounderPy", "=")
    lines += ["If SounderPy contributes to your work, please cite its Journal of Open Source Software article:", "",
              f"    {rst_text(REFS['sounderpy_joss'])}", "",
              "and the references below for the parameters and data you use.", "", "*" * 63, ""]

    lines += heading("Parameters", "=")
    lines += [".. _params-severe:", ""] + heading("Severe & general (``sounding_params``)", "-")
    lines += ["``spy.sounding_params(clean_data).calc()`` returns these values in four dictionaries: ``general``,",
              "``thermo``, ``kinem``, and ``intrp`` (the profile interpolated to 100 m levels). See",
              "the :doc:`helper tools </reference/helpertools>` page.", ""]
    for name, title, note in DICTS:
        lines += heading(title, "^") + table(name, note)

    lines += [".. _params-winter:", ""] + heading("Winter (``winter_params``)", "-")
    lines += ["``spy.winter_params(clean_data).calc()`` returns everything ``sounding_params`` does, plus a ``winter``",
              "dictionary with these values. The :doc:`winter soundings tutorial </tutorials/winter_soundings>` explains",
              "how they're used on the winter sounding.", ""]
    for name, title, note in WINTER_TOPICS:
        lines += heading(title, "^") + table(name, note)
    lines += ["*" * 63, ""]

    lines += heading("References", "=")
    for group, keys in GROUPS:
        label = "ref-winter" if group == "Winter parameters" else "ref-" + re.sub(r"\W+", "-", group.lower()).strip("-")
        lines += [f".. _{label}:", ""] + heading(group, "-")
        for key in keys:
            lines += [f".. _ref-{key}:", "", f"{rst_text(REFS[key])} *({TAGS[key]})*", ""]
    return "\n".join(lines).rstrip() + "\n"


def main():
    parser = argparse.ArgumentParser(description="Write the Parameters & References docs page.")
    parser.add_argument("--check", action="store_true", help="exit 1 if the page is out of date")
    args = parser.parse_args()
    page = build()
    if args.check:
        current = PAGE.read_text(encoding="utf-8") if PAGE.exists() else ""
        print("up to date" if current == page else f"{PAGE} is out of date")
        return 0 if current == page else 1
    PAGE.write_text(page, encoding="utf-8")
    print(f"wrote {PAGE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
