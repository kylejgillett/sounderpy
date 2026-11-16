"""Turn a docs build into redirect pages for the old docs address.

The SounderPy docs moved from kylejgillett.github.io/sounderpy/ to sounderpy.app/docs/. This takes a finished
Sphinx build and replaces every page with a small page that sends visitors (and search engines) to the same page
at the new address, keeping any #section. Everything else (images, CSS, notebooks) is removed. A 404.html catches
anything else under the old address.

    python docs/make_redirects.py docs/build/html
"""

import html
import sys
from pathlib import Path

NEW_BASE = "https://sounderpy.app/docs/"

PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SounderPy docs have moved</title>
<link rel="canonical" href="{url}">
<meta name="robots" content="noindex">
<meta http-equiv="refresh" content="0; url={url}">
<script>location.replace({url_js} + location.hash);</script>
</head>
<body style="font-family: sans-serif; background: #302d29; color: #f6ebd2; padding: 2rem;">
<p>The SounderPy docs have moved to <a style="color: #ddc386;" href="{url}">{url}</a>.</p>
</body>
</html>
"""

NOT_FOUND = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SounderPy docs have moved</title>
<meta name="robots" content="noindex">
<script>
  // /sounderpy/<path> on the old site is /docs/<path> on the new one
  var path = location.pathname.replace(/^\\/sounderpy\\/?/, "");
  location.replace("{base}" + path + location.search + location.hash);
</script>
</head>
<body style="font-family: sans-serif; background: #302d29; color: #f6ebd2; padding: 2rem;">
<p>The SounderPy docs have moved to <a style="color: #ddc386;" href="{base}">{base}</a>.</p>
</body>
</html>
"""


def main(build_dir):
    root = Path(build_dir)
    pages = [p for p in root.rglob("*.html")]
    keep = set()
    for page in pages:
        rel = page.relative_to(root).as_posix()
        url = NEW_BASE + ("" if rel == "index.html" else rel)
        page.write_text(PAGE.format(url=html.escape(url, quote=True), url_js=repr(url)))
        keep.add(page)
    (root / "404.html").write_text(NOT_FOUND.format(base=NEW_BASE))
    keep.add(root / "404.html")
    (root / ".nojekyll").write_text("")
    keep.add(root / ".nojekyll")
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_file() and path not in keep:
            path.unlink()
        elif path.is_dir() and not any(path.iterdir()):
            path.rmdir()
    print(f"{len(pages)} redirect pages written to {root}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "docs/build/html")
