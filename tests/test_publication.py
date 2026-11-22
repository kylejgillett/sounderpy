"""End-to-end publication figures: build_publication_sounding's panels and layouts, widths, options, and saved size.

Every test builds a figure offline from synthetic profiles (the map tests read Natural Earth and MetPy shapes, cached
after the first run), so the whole file is marked ``slow``.
"""

import matplotlib.pyplot as plt
import numpy.ma as ma
import pytest
from metpy.plots import SkewT

from sounderpy.plot import publication
from sounderpy.plot.publication import PUBLICATION_WIDTHS, build_publication_sounding

from helpers import make_clean_profile, make_second_profile


pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def two_profiles():
    return [make_clean_profile(), make_second_profile()]


def texts(ax):
    return [text.get_text() for text in ax.texts]


@pytest.mark.parametrize("panels", [
    ("skewt",),
    ("skewt", "hodograph"),
    ("skewt", "hodograph", "parameters"),
    ("parameters", "hodograph", "skewt"),
    ("map", "skewt", "parameters"),
    ("skewt", "hodograph", "map", "parameters"),
])
def test_panels_are_placed_by_importance(two_profiles, panels):
    fig, axes = build_publication_sounding(two_profiles, panels=panels, dpi=50)
    assert list(axes) == [panel for panel in publication.PUBLICATION_PANELS if panel in panels]
    assert isinstance(axes["skewt"], SkewT)
    assert fig.get_figwidth() == pytest.approx(PUBLICATION_WIDTHS["double"])
    fig.canvas.draw()
    boxes = {panel: (ax.ax if panel == "skewt" else ax).get_position() for panel, ax in axes.items()}
    inches = {panel: (box.width * fig.get_figwidth(), box.height * fig.get_figheight()) for panel, box in boxes.items()}
    areas = [inches[panel][0] * inches[panel][1] for panel in ("skewt", "hodograph", "map") if panel in inches]
    assert areas == sorted(areas, reverse=True)                                     # more important, bigger
    layout = publication.PUBLICATION_LAYOUTS[tuple(axes)]
    for panel in ("hodograph", "map"):                                              # filling their spaces, any shape
        if panel in inches:
            assert inches[panel][0] == pytest.approx(layout[panel][2], rel=0.03)
            assert inches[panel][1] == pytest.approx(layout[panel][3], rel=0.03)
    for a, b in [(a, b) for a in boxes.values() for b in boxes.values() if a is not b]:   # none overlap
        assert a.x1 <= b.x0 + 1e-6 or b.x1 <= a.x0 + 1e-6 or a.y1 <= b.y0 + 1e-6 or b.y1 <= a.y0 + 1e-6


def test_default_puts_the_table_under_the_skewt_and_hodograph(two_profiles):
    fig, axes = build_publication_sounding(two_profiles, dpi=50)
    skew, hodograph, table = (axes["skewt"].ax.get_position(), axes["hodograph"].get_position(),
                              axes["parameters"].get_position())
    assert skew.x1 < hodograph.x0 and skew.y1 == pytest.approx(hodograph.y1)    # side by side, tops lined up
    assert min(skew.y0, hodograph.y0) > table.y1                                # the table underneath
    assert table.x0 <= skew.x0 and table.x1 >= hodograph.x1 - 0.01             # across the figure


def test_three_panels_put_a_big_skewt_beside_two_stacked(two_profiles):
    fig, axes = build_publication_sounding(two_profiles, panels=("skewt", "hodograph", "map"), dpi=50)
    skew, hodograph, map_ = (axes["skewt"].ax.get_position(), axes["hodograph"].get_position(),
                             axes["map"].get_position())
    assert skew.x1 < hodograph.x0 and hodograph.x0 == pytest.approx(map_.x0)    # skew-t left, the others right
    assert hodograph.y0 > map_.y1                                               # stacked
    assert skew.height > 1.5 * hodograph.height                                 # the skew-t is the big one


def test_four_panels_put_the_table_in_a_strip_across_the_bottom(two_profiles):
    fig, axes = build_publication_sounding(two_profiles, panels=publication.PUBLICATION_PANELS, dpi=50)
    skew, hodograph, map_, table = (axes["skewt"].ax.get_position(), axes["hodograph"].get_position(),
                                    axes["map"].get_position(), axes["parameters"].get_position())
    assert skew.x1 < hodograph.x0 and hodograph.x0 == pytest.approx(map_.x0)    # skew-t left, the others right
    assert hodograph.y0 > map_.y1                                               # the hodograph over the map
    assert min(skew.y0, map_.y0) > table.y1                                     # the table underneath
    assert table.x0 <= skew.x0 and table.x1 >= map_.x1 - 0.01                  # across the figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    names = [text for text in axes["parameters"].texts if text.get_text() in ("MLCAPE", "STP")]
    boxes = sorted((text.get_window_extent(renderer) for text in names), key=lambda box: box.x0)
    assert len(boxes) == 2 and boxes[0].y0 == pytest.approx(boxes[1].y0)      # the names in one row
    texts = sorted((text.get_window_extent(renderer) for text in axes["parameters"].texts if text.get_text().strip()
                    and round(text.get_window_extent(renderer).y0) == round(boxes[0].y0)), key=lambda box: box.x0)
    assert all(left.x1 < right.x0 for left, right in zip(texts, texts[1:]))    # the columns don't overlap


@pytest.mark.parametrize(("options", "message"), [
    ({"panels": ("hodograph", "map")}, "'panels' must name the skew-t"),
    ({"panels": ("skewt", "radar")}, "'panels' must name the skew-t"),
    ({"panels": ("skewt", "parameters"), "parameters": None}, "the 'parameters' panel needs at least one parameter"),
    ({"parameters": tuple(publication.PUBLICATION_PARAMETERS)}, "too wide for the table"),
])
def test_panel_option_errors(two_profiles, options, message):
    with pytest.raises(ValueError, match=message):
        build_publication_sounding(two_profiles, dpi=50, **options)


def test_width_in_inches_and_bad_settings(two_profiles):
    fig, _ = build_publication_sounding(two_profiles, width=5, dpi=50)
    assert fig.get_figwidth() == pytest.approx(5)
    for options, message in (({"width": "wide"}, "'width' must be one of"),
                             ({"width": -1}, "'width' must be one of"), ({"wind_units": "mph"}, "'wind_units' must be"),
                             ({"colors": ["red"]}, "2 profiles need at least that many colors"),
                             ({"ylim": (500, 1000)}, "'ylim' must be")):
        with pytest.raises(ValueError, match=message):
            build_publication_sounding(two_profiles, dpi=50, **options)


def test_one_profile_dict(two_profiles):
    fig, axes = build_publication_sounding(two_profiles[0], dpi=50)
    assert len(axes["skewt"].ax.get_legend_handles_labels()[0]) == 0  # no labels, no legend
    assert axes["skewt"].ax.get_legend() is None and not fig.legends


def test_labels_make_a_legend(two_profiles):
    fig, _ = build_publication_sounding(two_profiles, labels=["RAOB", "HRRR F03"], dpi=50)
    assert [text.get_text() for text in fig.legends[0].get_texts()] == ["RAOB", "HRRR F03"]


def test_panel_labels(two_profiles):
    _, axes = build_publication_sounding(two_profiles, dpi=50)
    assert "(a)" in texts(axes["skewt"].ax) and "(b)" in texts(axes["hodograph"])
    _, axes = build_publication_sounding(two_profiles, panel_labels=None, dpi=50)
    assert "(a)" not in texts(axes["skewt"].ax)


def test_parameter_table(two_profiles):
    _, axes = build_publication_sounding(two_profiles, parameters=("mlcape", "stp"), dpi=50)
    cells = texts(axes["parameters"])
    assert "MLCAPE" in cells and "STP" in cells and r"J kg$^{-1}$" in cells
    # each value is right of its row's name, and the second profile's value right of the first's
    by_row = {}
    for text in axes["parameters"].texts:
        by_row.setdefault(round(text.get_position()[1], 3), []).append(text)
    for row in by_row.values():
        xs = [text.get_position()[0] for text in row]
        assert xs == sorted(xs)


def test_missing_parameter_reads_na():
    params = {"thermo": {"mlcape": ma.masked}}
    assert publication._format_publication_parameter(params, "mlcape")[2] == publication.MISSING_VALUE


def test_zoomed_skewt_fits_its_temperature_range(two_profiles):
    _, axes = build_publication_sounding(two_profiles, ylim=(1000, 600), panels=("skewt",), dpi=50)
    assert axes["skewt"].ax.get_xlim() != pytest.approx((-40, 40))
    assert axes["skewt"].ax.get_ylim() == pytest.approx((1000, 600))
    # a range given in the style is kept
    _, axes = build_publication_sounding(two_profiles, ylim=(1000, 600), style={"temperature_range": (-10, 40)},
                                         panels=("skewt",), dpi=50)
    assert axes["skewt"].ax.get_xlim()[0] == pytest.approx(-10, abs=0.5)


def test_watermark_is_off_unless_asked(two_profiles):
    fig, _ = build_publication_sounding(two_profiles, dpi=50)
    assert not fig.texts
    fig, _ = build_publication_sounding(two_profiles, dpi=50, watermark={"title": "Gillett et al. (2026)", "subtitle": ""})
    assert [text.get_text() for text in fig.texts] == ["Gillett et al. (2026)"]


def test_map_needs_a_location(two_profiles):
    profile = make_clean_profile()
    del profile["site_info"]["site-latlon"]
    with pytest.raises(ValueError, match="the 'map' panel needs each profile's location"):
        build_publication_sounding(profile, panels=("skewt", "map"), dpi=50)


def test_save(two_profiles, tmp_path):
    path = tmp_path / "figure.pdf"
    build_publication_sounding(two_profiles, save=True, filename=str(path), dpi=50)
    assert path.read_bytes()[:4] == b"%PDF"


@pytest.mark.parametrize(("width", "inches"), [("single", 3.2), ("double", 6.5), ("ams_two_column", 5.5),
                                               ("agu_single", 3.35), ("agu_double", 6.7), (4.0, 4.0)])
def test_journal_widths(two_profiles, width, inches):
    fig, _ = build_publication_sounding(two_profiles, width=width, panels=("skewt", "hodograph"), dpi=50)
    assert fig.get_figwidth() == pytest.approx(inches)


def test_small_panels_warn(two_profiles):
    import warnings
    from sounderpy.utils import SounderPyWarning

    # every panel at the single-column width is too small to read
    with pytest.warns(SounderPyWarning, match="A panel is only"):
        build_publication_sounding(two_profiles, width="single", panels=("skewt", "hodograph", "map"), dpi=50)
    with warnings.catch_warnings():
        warnings.simplefilter("error", SounderPyWarning)
        build_publication_sounding(two_profiles, width="single", panels=("skewt",), dpi=50)
        build_publication_sounding(two_profiles, dpi=50)


def hodograph_markers(ax):
    """Centers (display pixels) of the white height-marker circles on a hodograph."""
    ax.figure.canvas.draw()
    return [ax.transData.transform(line.get_xydata()[0]) for line in ax.lines if line.get_gid() == "height_marker"]


def test_height_labels_are_s_then_km(two_profiles):
    _, axes = build_publication_sounding(two_profiles[0], dpi=100)
    assert {"S", "1", "3", "6", "12"} <= set(texts(axes["hodograph"]))
    _, axes = build_publication_sounding(two_profiles[0], dpi=100, style={"surface_label": "SFC"}, hodo_heights=(0, 500, 9000))
    assert {"SFC", ".5", "9"} <= set(texts(axes["hodograph"]))


def test_height_markers_never_overlap(two_profiles):
    # three profiles that differ little, so many markers fall in the same place
    import copy
    from metpy.units import units
    third = copy.deepcopy(two_profiles[0])
    third["u"] = third["u"] + 1 * units.kts
    fig, axes = build_publication_sounding(two_profiles + [third], dpi=100)
    centers = hodograph_markers(axes["hodograph"])
    size = publication.PUBLICATION_STYLE["height_marker_size"] * fig.dpi / 72
    assert len(centers) == 15
    for i, a in enumerate(centers):
        for b in centers[i + 1:]:
            assert ((a - b) ** 2).sum() ** 0.5 >= size


# ---------------------------------------------------------------------------------------------------------------
# MAP: radar and the user's own fields (map tests read Natural Earth / MetPy shapes, cached after the first run)
# ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize(("options", "message"), [({"panels": ("skewt", "map"), "map_radar": "nexrad"}, "'map_radar' must be"),
                                                  ({"map_radar": "mosaic"}, "'map_radar' needs the 'map' panel")])
def test_radar_option_errors(two_profiles, options, message):
    with pytest.raises(ValueError, match=message):
        build_publication_sounding(two_profiles, dpi=50, **options)


@pytest.mark.parametrize(("radar", "found"), [("mosaic", True), ("single", True), ("mosaic", False)])
def test_radar_caption_and_colorbar(monkeypatch, two_profiles, radar, found):
    import numpy as np
    requested = []

    def fake_radar(ax, profile, *args, **kwargs):
        # what SounderPy's map radar does: a reflectivity mesh, and a boxed label on the map
        requested.append(profile["site_info"]["site-latlon"])
        if found:
            ax.pcolormesh(np.linspace(-99, -95, 5), np.linspace(46, 50, 5), np.full((4, 4), 40.0), vmin=-32, vmax=95)
        ax.text(0.5, 0.03, "KTLX: 2024-05-06 | 21:04" if found else "RADAR UNAVAILABLE", transform=ax.transAxes)
        return found

    monkeypatch.setattr(publication, "_plot_map_mosaic_radar", fake_radar)
    monkeypatch.setattr(publication, "_plot_map_single_site_radar", fake_radar)
    fig, axes = build_publication_sounding(two_profiles, panels=("skewt", "hodograph", "map"), map_radar=radar, dpi=50)
    on_map = [text.get_text() for text in axes["map"].texts]
    assert not any("|" in text for text in on_map)   # the map helper's boxed label is gone
    if found:
        # a colorbar, with the scan time as the second line of its label
        (colorbar,) = axes["map"].child_axes
        assert colorbar.get_xlabel() == "reflectivity (dBZ)\nKTLX 2024-05-06 21:04"
        assert "KTLX 2024-05-06 21:04" not in on_map
    else:
        # no radar: no colorbar and no label on the map (the map radar's own SounderPyWarning says why)
        assert not axes["map"].child_axes
        assert "RADAR UNAVAILABLE" not in on_map
    assert len(requested) == 1


def test_radar_scan_time_above_the_map_without_a_colorbar(monkeypatch, two_profiles):
    import numpy as np

    def fake_radar(ax, profile, *args, **kwargs):
        ax.pcolormesh(np.linspace(-99, -95, 5), np.linspace(46, 50, 5), np.full((4, 4), 40.0), vmin=-32, vmax=95)
        ax.text(0.5, 0.03, "KTLX: 2024-05-06 | 21:04", transform=ax.transAxes)
        return True

    monkeypatch.setattr(publication, "_plot_map_single_site_radar", fake_radar)
    _, axes = build_publication_sounding(two_profiles, panels=("skewt", "hodograph", "map"), map_radar="single",
                                         style={"radar_colorbar": False}, dpi=50)
    caption = [text for text in axes["map"].texts if text.get_text() == "KTLX 2024-05-06 21:04"]
    assert len(caption) == 1 and caption[0].get_position()[1] > 1 and not axes["map"].child_axes


def test_user_contours_sit_under_the_map_lines(two_profiles):
    import cartopy.crs as ccrs
    import numpy as np
    from cartopy.mpl.feature_artist import FeatureArtist

    fig, axes = build_publication_sounding(two_profiles[0], panels=("skewt", "map"), dpi=50)
    lons, lats = np.meshgrid(np.linspace(-100, -94, 30), np.linspace(45, 51, 30))
    filled = axes["map"].contourf(lons, lats, lons + lats, transform=ccrs.PlateCarree())
    lines = [artist for artist in axes["map"].get_children() if isinstance(artist, FeatureArtist)]
    assert lines and all(artist.get_zorder() > filled.get_zorder() for artist in lines)   # no state fill over the field
    # the map is centered on the profile
    west, east, south, north = axes["map"].get_extent(ccrs.PlateCarree())
    lat, lon = two_profiles[0]["site_info"]["site-latlon"]
    assert (west + east) / 2 == pytest.approx(lon, abs=0.01) and (south + north) / 2 == pytest.approx(lat, abs=0.01)
    fig.canvas.draw()


def test_returned_figure_keeps_type42_fonts_when_saved_later(two_profiles, tmp_path):
    fig, _ = build_publication_sounding(two_profiles, dpi=50)
    for save in (fig.savefig, plt.savefig):
        path = tmp_path / f"{save.__qualname__}.pdf"
        plt.figure(fig.number)
        save(path)
        pdf = path.read_bytes()
        assert b"/Type3" not in pdf
        assert b"/CIDFontType2" in pdf or b"/TrueType" in pdf


@pytest.mark.parametrize(("options", "message"), [
    ({"parcel": "xx"}, "'parcel' must be"),
    ({"labels": ["one"]}, "2 profiles need at least that many labels"),
    ({"storm_motion_profile": 5}, "'storm_motion_profile' must be a profile's number"),
    ({"parameter_profiles": [3]}, "'parameter_profiles' must list profile numbers"),
    ({"parameters": ("mlcape", "nope")}, "Unknown parameter 'nope'"),
    ({"parameters": (("X", "thermo", "mlcape"),)}, "A parameter tuple needs"),
    ({"parameters": (("X", "thermo", "not_a_key", 1, "", 0),)}, "sounding_params has no 'thermo' value named 'not_a_key'"),
])
def test_bad_options_say_what_is_wrong(two_profiles, options, message):
    with pytest.raises(ValueError, match=message):
        build_publication_sounding(two_profiles, dpi=50, **options)


def test_cin_is_listed_as_negative(two_profiles):
    from sounderpy.plot.publication import _format_publication_parameter
    _, _, value = _format_publication_parameter({"thermo": {"mlcin": -46.4}}, "mlcin")
    assert value == "-46"


@pytest.mark.parametrize(("panels", "labels"), [(("skewt", "hodograph", "map"), ["00Z ABR RAOB", "HRRR F03 Ashley, ND"]),
                                                (("skewt", "hodograph"), ["RAOB", "HRRR F03"]),
                                                (publication.PUBLICATION_PANELS, ["RAOB", "HRRR F03"])])
def test_legend_is_one_row_under_every_panel(two_profiles, panels, labels):
    fig, axes = build_publication_sounding(two_profiles, labels=labels, panels=panels, dpi=50)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend = fig.legends[0]
    extent = legend.get_window_extent(renderer)
    assert extent.y1 <= min(ax.get_tightbbox(renderer).y0 for ax in fig.axes) + 1        # under every panel
    assert extent.x0 + extent.x1 == pytest.approx(fig.bbox.width, abs=2)                  # centered
    assert extent.y0 >= 0 and extent.x0 >= 0 and extent.x1 <= fig.bbox.width + 1        # inside the figure
    texts = legend.get_texts()
    assert [text.get_text() for text in texts] == labels
    assert len({round(text.get_window_extent(renderer).y0) for text in texts}) == 1   # one row


@pytest.mark.parametrize("options", [
    {},
    {"panels": ("skewt", "parameters"), "parameters": tuple(publication.PUBLICATION_PARAMETERS)[:8],
     "labels": ["a long profile name", "HRRR F03"]},
    {"panels": publication.PUBLICATION_PANELS, "labels": ["00Z ABR RAOB", "HRRR F03"]},
    {"panels": ("skewt", "hodograph", "map"), "labels": ["00Z ABR RAOB", "HRRR F03"]},
    {"width": "single", "panels": ("skewt",), "labels": ["RAOB", "HRRR"]},
    {"panels": ("skewt",), "ylim": (1000, 500)},
])
def test_everything_fits_inside_a_figure_of_the_journal_width(two_profiles, tmp_path, options):
    from PIL import Image
    fig, axes = build_publication_sounding(two_profiles, dpi=50, save=True, filename=str(tmp_path / "fig.png"), **options)
    width, height = fig.get_size_inches()
    assert width == pytest.approx(PUBLICATION_WIDTHS[options.get("width", "double")])
    drawn = fig.get_tightbbox(fig.canvas.get_renderer())
    assert drawn.x0 >= 0 and drawn.y0 >= 0 and drawn.x1 <= width + 0.01 and drawn.y1 <= height + 0.01
    if "parameters" in axes:   # every row of the table is inside the figure
        renderer = fig.canvas.get_renderer()
        assert min(text.get_window_extent(renderer).y0 for text in axes["parameters"].texts) >= 0
    saved = Image.open(tmp_path / "fig.png").size   # saved at exactly the figure's size
    assert saved[0] == pytest.approx(width * 50, abs=1) and saved[1] == pytest.approx(height * 50, abs=1)


def test_storm_relative_vectors_are_drawn_over_the_traces(two_profiles):
    from matplotlib.colors import to_rgba
    fig, axes = build_publication_sounding(two_profiles, dpi=100, sr_vector_heights=(1000, 6000))
    hodograph = axes["hodograph"]
    green = to_rgba(publication.PUBLICATION_STYLE["sr_vector_color"])
    vectors = [text for text in hodograph.texts if getattr(text, "arrow_patch", None) is not None
               and to_rgba(text.arrow_patch.get_edgecolor()) == green]
    traces = [line for line in hodograph.lines if line.get_gid() != "height_marker"
              and line.get_linewidth() == publication.PUBLICATION_STYLE["hodograph_linewidth"]]   # not marker leader lines
    trace_zorder = max(line.get_zorder() for line in traces)
    assert vectors and all(vector.get_zorder() > trace_zorder for vector in vectors)


@pytest.mark.parametrize("panels", [("skewt", "hodograph"), ("skewt", "hodograph", "map"), publication.PUBLICATION_PANELS])
@pytest.mark.parametrize("profiles", [2, 4])
def test_panel_labels_sit_clear_of_the_barbs_and_every_panel(two_profiles, panels, profiles):
    fig, axes = build_publication_sounding((two_profiles * 2)[:profiles], panels=panels, dpi=100)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    barbs_right = axes["skewt"].ax.get_tightbbox(renderer).x1
    labels = [text.get_window_extent(renderer) for ax in fig.axes for text in ax.texts if text.get_gid() == "panel label"]
    assert len(labels) == len(panels)
    panels_boxes = [(ax.ax if name == "skewt" else ax).get_window_extent(renderer) for name, ax in axes.items()]
    skew_box = panels_boxes[0]
    for label in labels:
        if label.x0 > skew_box.x1:                                       # beside the skew-t: right of its wind barbs
            assert label.x0 > barbs_right
        for box in panels_boxes:                                         # over no panel
            assert label.x1 <= box.x0 + 1 or label.x0 >= box.x1 - 1 or label.y0 >= box.y1 - 1 or label.y1 <= box.y0 + 1
