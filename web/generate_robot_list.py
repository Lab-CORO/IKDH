#!/usr/bin/env python3
"""
generate_robot_list.py — Regenerate the ROBOTS catalog array embedded in
web/ikdh.html from every file in robots/*.yaml.

Each entry's display name is taken verbatim from that file's `name:` field
(matching what RoboDK itself calls the robot), paired with its filename.

Usage (from the repository root):
    python3 web/generate_robot_list.py
"""

import glob
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROBOTS_DIR = os.path.join(ROOT, "robots")
HTML_PATH = os.path.join(ROOT, "web", "ikdh.html")

START_MARKER = "const ROBOTS = ["
END_MARKER = "];"


def _js_string(s):
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def load_entries():
    entries = []
    for path in sorted(glob.glob(os.path.join(ROBOTS_DIR, "*.yaml"))):
        with open(path) as f:
            for line in f:
                if line.startswith("name:"):
                    name = line.split(":", 1)[1].strip()
                    entries.append((name, os.path.basename(path)))
                    break
    return entries


def build_block(entries):
    lines = [START_MARKER]
    for name, yaml in entries:
        lines.append(f"  {{ name: {_js_string(name)}, yaml: {_js_string(yaml)} }},")
    lines.append(END_MARKER)
    return "\n".join(lines)


def main():
    entries = load_entries()
    block = build_block(entries)

    with open(HTML_PATH) as f:
        html = f.read()

    start = html.index(START_MARKER)
    end = html.index(END_MARKER, start) + len(END_MARKER)
    new_html = html[:start] + block + html[end:]

    with open(HTML_PATH, "w") as f:
        f.write(new_html)

    print(f"{len(entries)} robots written to {HTML_PATH}")


if __name__ == "__main__":
    main()
