#!/usr/bin/env python3
"""
build_standalone.py -- regenerate silent_mutation_standalone.html from the core.

The standalone file is a snapshot of the silent_mutation/ core: it embeds the
core (gzipped + base64) and runs it in the browser with Pyodide, so no Python
install is needed to use the tool. After you edit the core, run this to refresh
the snapshot:

    python build_standalone.py

It reads silent_mutation/webtool/index.html (the server UI), swaps the two
fetch('/api/...') calls for in-browser Pyodide calls, embeds the core + codon
table, and writes silent_mutation_standalone.html at the repo root.

DeepPrime is intentionally excluded (genet is Python-server only): deepprime_runner
and server.py are not embedded, and the DeepPrime checkbox reports "unavailable".
"""

from __future__ import annotations

import base64
import io
import pathlib
import tarfile

ROOT = pathlib.Path(__file__).resolve().parent
PKG = ROOT / "silent_mutation"
PYODIDE_VER = "314.0.2"
CDN = f"https://cdn.jsdelivr.net/pyodide/v{PYODIDE_VER}/full/"

# Core files to embed. deepprime_runner.py (genet) and server.py (Flask) are
# excluded -- they can't run under Pyodide and aren't needed for the standalone
# analyze/verify paths.
EMBED = [
    "silent_mutation/__init__.py",
    "silent_mutation/core/__init__.py",
    "silent_mutation/core/types.py",
    "silent_mutation/core/codon_utils.py",
    "silent_mutation/core/pam_finder.py",
    "silent_mutation/core/silent_finder.py",
    "silent_mutation/core/pegrna_builder.py",
    "silent_mutation/core/verify.py",
    "silent_mutation/io/__init__.py",
    "silent_mutation/io/sequence_loader.py",
    "silent_mutation/io/genome_loader.py",
    "silent_mutation/webtool/__init__.py",
    "silent_mutation/webtool/api.py",
    "data/reference/codon_table.csv",
]


def build_pkg_b64() -> str:
    """Reproducible gzipped tar of the embed set, base64-encoded."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz", compresslevel=9) as tar:
        for rel in EMBED:
            src = ROOT / rel
            if not src.exists():
                raise SystemExit(f"missing file to embed: {rel}")
            info = tarfile.TarInfo(name="./" + rel)
            data = src.read_bytes()
            info.size = len(data)
            info.mtime = 0
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    return base64.b64encode(buf.getvalue()).decode("ascii")


def transform_html(src_html: str, b64: str) -> str:
    s = src_html

    # Pyodide runtime <script> before </head>
    assert s.count("</head>") == 1
    s = s.replace("</head>", f'  <script src="{CDN}pyodide.js"></script>\n</head>')

    # embedded package (base64 text block) + loading banner right after <body>
    assert s.count("<body>") == 1
    banner = (
        "<body>\n"
        '<script type="text/plain" id="pkgb64">' + b64 + "</script>\n"
        '<div id="pyloading" style="position:fixed;top:0;left:0;right:0;z-index:9999;'
        "background:#0b7285;color:#fff;padding:8px 14px;font:14px system-ui,sans-serif;"
        'text-align:center">\u25CF Loading analysis engine\u2026 first load fetches the '
        "Python runtime (~10\u201320s, needs internet). Subsequent opens are cached &amp; fast.</div>"
    )
    s = s.replace("<body>", banner)

    # bootstrap + runPy, right after the $ helper
    anchor = "const $=id=>document.getElementById(id);"
    boot = anchor + "\n" + r"""
/* self-contained engine: run the real Python core in-browser via Pyodide */
const PY_INDEX_URL = "%CDN%";
let _pyodide = null;
const pyReady = (async () => {
  const banner = $('pyloading');
  try {
    _pyodide = await loadPyodide({ indexURL: PY_INDEX_URL });
    const b64 = $('pkgb64').textContent.trim();
    _pyodide.FS.writeFile('/tmp/pkg.tar.gz', Uint8Array.from(atob(b64), c => c.charCodeAt(0)));
    await _pyodide.runPythonAsync(`
import tarfile, sys, os, json as _json
os.makedirs('/app', exist_ok=True)
with tarfile.open('/tmp/pkg.tar.gz', 'r:gz') as _t:
    _t.extractall('/app')
if '/app' not in sys.path:
    sys.path.insert(0, '/app')
from silent_mutation.webtool import api as _api
`);
    if (banner) banner.style.display = 'none';
  } catch (e) {
    if (banner) banner.innerHTML = 'Could not load the analysis engine: ' + e +
      ' \u2014 check your internet connection and reload the page.';
    throw e;
  }
})();

async function runPy(fn, body) {
  await pyReady;
  _pyodide.globals.set('_body_json', JSON.stringify(body));
  const out = _pyodide.runPython(`_json.dumps(_api.${fn}(_json.loads(_body_json)))`);
  return JSON.parse(out);
}
""".replace("%CDN%", CDN)
    assert s.count(anchor) == 1
    s = s.replace(anchor, boot)

    # swap the two backend fetches for in-browser calls
    analyze_fetch = (
        "    const r=await fetch('/api/analyze',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});\n"
        "    data=await r.json();"
    )
    verify_fetch = (
        "    const r=await fetch('/api/verify',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});\n"
        "    data=await r.json();"
    )
    assert s.count(analyze_fetch) == 1, "analyze fetch block not found"
    assert s.count(verify_fetch) == 1, "verify fetch block not found"
    s = s.replace(analyze_fetch, "    data=await runPy('run_analyze', body);")
    s = s.replace(verify_fetch, "    data=await runPy('run_verify', body);")
    assert "fetch('/api/" not in s, "a backend fetch survived"
    return s


def main() -> None:
    src = (PKG / "webtool" / "index.html").read_text(encoding="utf-8")
    b64 = build_pkg_b64()
    out = transform_html(src, b64)
    dest = ROOT / "silent_mutation_standalone.html"
    dest.write_text(out, encoding="utf-8")
    print(f"wrote {dest}  ({len(out):,} bytes, embedded package {len(b64):,} b64 chars)")


if __name__ == "__main__":
    main()
