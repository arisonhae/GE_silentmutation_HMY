#!/usr/bin/env python3
"""
One-command launcher for the Silent Mutation web tool (server build).

    python run.py

Starts the Flask server on http://127.0.0.1:8502 and opens your browser.
The standalone analyze + verify path needs only Flask (see requirements.txt).
DeepPrime ranking additionally needs the genet stack (requirements-deepprime.txt);
without it the tool still works and DeepPrime just reports "unavailable".

For a zero-install option, just open silent_mutation_standalone.html in a browser.
"""
import os
import sys
import threading
import time
import webbrowser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PORT = int(os.environ.get("PORT", "8502"))


def _open_browser():
    time.sleep(1.3)
    webbrowser.open("http://127.0.0.1:%d/" % PORT)


if __name__ == "__main__":
    from silent_mutation.webtool.server import app
    threading.Thread(target=_open_browser, daemon=True).start()
    print("Serving the Silent Mutation web tool at http://127.0.0.1:%d/  (Ctrl+C to stop)" % PORT)
    app.run(host="127.0.0.1", port=PORT, debug=False)
