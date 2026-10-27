"""Track-changes view of the student's vignettes.

Compares each vignette as Lukas Bolsinger wrote it (git commit 2de18bc,
"Add the student's four vignettes, unchanged") with the current version in
vignettes/, word by word, and writes one HTML page in which deleted text is
struck through in red and inserted text is underlined in green, like tracked
changes in Word.

Run from the package root:  python validation/09_vignette_track_changes.py
Output: validation/output/vignette_track_changes.html
"""

import difflib
import html
import re
import subprocess
from pathlib import Path

ORIGINAL_COMMIT = "2de18bc"
VIGNETTES = [
    ("introduction.Rmd", "introduction_p_ggeffects.Rmd"),
    ("parameterized_models.Rmd", "introduction_parameterized_model.Rmd"),
    ("heart_failure_egfr.Rmd", "p_ggeffects_demonstration.Rmd"),
    ("heart_failure_renal_index.Rmd", "p_ggeffects_heart_failure.Rmd"),
]
OUTPUT = Path("validation/output/vignette_track_changes.html")


def original(path):
    return subprocess.run(["git", "show", f"{ORIGINAL_COMMIT}:vignettes/{path}"],
                          capture_output=True, text=True, encoding="utf-8", check=True).stdout


def tokens(text):
    # words, runs of spaces, and single line breaks/punctuation
    return re.findall(r"\n|[ \t]+|\w+|[^\w\s]", text)


def redline(old, new):
    a, b = tokens(old), tokens(new)
    out, counts = [], {"deleted": 0, "inserted": 0}
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if op == "equal":
            out.append(html.escape("".join(a[i1:i2])))
            continue
        if op in ("delete", "replace"):
            out.append(f"<del>{html.escape(''.join(a[i1:i2]))}</del>")
            counts["deleted"] += sum(1 for t in a[i1:i2] if t.strip())
        if op in ("insert", "replace"):
            out.append(f"<ins>{html.escape(''.join(b[j1:j2]))}</ins>")
            counts["inserted"] += sum(1 for t in b[j1:j2] if t.strip())
    return "".join(out), counts


def main():
    sections, toc = [], []
    for current, student_name in VIGNETTES:
        old = original(current)
        new = Path("vignettes", current).read_text(encoding="utf-8")
        body, counts = redline(old, new)
        anchor = current.replace(".", "-")
        toc.append(f'<li><a href="#{anchor}">{current}</a> (student: {student_name}) &mdash; '
                   f'{counts["deleted"]} words deleted, {counts["inserted"]} inserted</li>')
        sections.append(f'<section id="{anchor}"><h2>{current}</h2>'
                        f'<p class="meta">Student\'s original: vignettes/{student_name} in p-ggeffects</p>'
                        f'<pre>{body}</pre></section>')

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Vignette track changes</title>
<style>
  :root {{ --bg: #ffffff; --ink: #1f2328; --muted: #59636e; --del: #b42318; --del-bg: #fdecea;
          --ins: #116329; --ins-bg: #e6f4ea; --rule: #d0d7de; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg: #0d1117; --ink: #e6edf3; --muted: #9198a1; --del: #ff7b72; --del-bg: #3d1214;
            --ins: #56d364; --ins-bg: #12331c; --rule: #30363d; }} }}
  body {{ background: var(--bg); color: var(--ink); font: 15px/1.5 system-ui, sans-serif;
         max-width: 980px; margin: 0 auto; padding: 24px 16px; }}
  h1 {{ font-size: 1.5rem; }} h2 {{ margin-top: 2.5rem; border-bottom: 1px solid var(--rule); }}
  .meta, .legend {{ color: var(--muted); }}
  pre {{ white-space: pre-wrap; word-wrap: break-word; font: 13px/1.55 ui-monospace, monospace; }}
  del {{ color: var(--del); background: var(--del-bg); text-decoration: line-through; }}
  ins {{ color: var(--ins); background: var(--ins-bg); text-decoration: underline; }}
</style></head><body>
<h1>Vignette track changes: p-ggeffects &rarr; termeffects</h1>
<p class="legend">The student's text (git commit {ORIGINAL_COMMIT}) compared with the current
termeffects vignettes. <del>Struck through</del>: removed. <ins>Underlined</ins>: added.
Everything else is unchanged.</p>
<ul>{''.join(toc)}</ul>
{''.join(sections)}
</body></html>
""", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
