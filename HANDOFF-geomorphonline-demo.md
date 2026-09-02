# Handoff: putting a corestone weathering demo on GeomorphOnline

Written 2026-09-02 by the session that built the first such exercise, the GRLP
gravel-river long profile, now live at

<https://geomorphonline.github.io/exercises/gravel-river-long-profile/>

Everything below was learned by building that one and deploying it. Follow it and
you should reach a working demo without rediscovering any of it. Read the
**Traps** section before writing code; most of it is not guessable.

---

## 0. Can corestone run in the browser at all?

**Yes. This is verified, not assumed** (2026-09-02):

- `pip wheel . --no-deps` produces `corestone-0.1.0.dev0-py3-none-any.whl`, a
  pure-Python wheel. Pyodide can install it.
- `artesian check numpy scipy matplotlib` reports all three as platform-only on
  PyPI **but bundled by Pyodide**, which is the good outcome.

Re-run `artesian check <your imports>` if you add a dependency. Anything with
compiled extensions that Pyodide does not bundle is a hard blocker: it needs a
wheel built for Emscripten, and for most scientific packages none exists.

The second constraint is speed. A model stepping in milliseconds animates
smoothly at 30 fps. If a weathering step takes seconds, do not animate it:
recompute on slider change instead, which is a perfectly good demo and avoids
the whole animation loop.

---

## 1. The pipeline

[artesian](https://github.com/MNiMORPH/artesian) compiles a Python model plus a
[Panel](https://panel.holoviz.org) interface into a self-contained WebAssembly
page that runs entirely in the reader's browser via Pyodide. No server, no
install, no accounts.

```sh
pip install "artesian @ git+https://github.com/MNiMORPH/artesian@main"
```

Not on PyPI yet (tracked as artesian#3), hence the git ref.

Build and view locally:

```sh
artesian build path/to/weathering_app.py -o _build -p . -r numpy -r scipy --serve
```

- `-p .` wheels **your** source tree and ships it, so the demo matches the code
  you are looking at.
- `-r` names things resolved in the browser: Pyodide-bundled packages, or
  anything with a pure-Python wheel.
- `--serve` is how to view it. Opening the page over `file://` trips the
  browser's cross-origin rules for web workers and fails opaquely.

---

## 2. Writing the app

Start from `examples/hillslope.py` in the artesian repository. It is a complete,
deliberately small demo and it is the file to copy. The GRLP one
(`~/models/GRLP/interactive_demo/grlp_panel.py`) is the realistic worked example.

Shape of it:

```python
import panel as pn
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from artesian.live import animator, reset_button, responsive

import corestone

# `sim`, NEVER `state`: panel exports pn.state and shadowing it fails silently.
sim = {"model": make_model(), "t": 0.}

def step():
    sim["model"].advance(dt)
    source.data = {"x": ..., "y": ...}

fig = figure(height=380, width=680, ...)
responsive(fig)                       # see section 3

pn.Column(
    pn.pane.Markdown("...", sizing_mode="stretch_width"),
    pn.Row(animator(step), reset_button(do_reset)),
    slider_a, slider_b,
    fig,
    sizing_mode="stretch_width",
).servable(title="...")
```

`animator()` gives a play/pause toggle wired to a periodic callback;
`reset_button()` a reset. That is all artesian offers for the app layer, on
purpose. It does not abstract your model or your plots, and you should build the
figure yourself.

---

## 3. Sizing: scale the app, do not merely stretch it

This took several iterations on GRLP and is the least obvious part.

**Responsive is not enough.** If elements merely fill the width, each keeps its
own intrinsic size: the plot grows without limit while 16 px text and ~18 px
slider handles stay put. On a wide screen the controls become small and fiddly
beside the model. What is wanted is *scaling* -- everything enlarging together,
the way zooming a PDF does.

So the app declares a **design width** and the embedding page scales the whole
thing above it:

```python
DESIGN_WIDTH = 900     # lay the app out to look right at this width
SLIDER_WIDTH = 520     # sliders stop here rather than spanning the pane

slider = pn.widgets.FloatSlider(..., sizing_mode="stretch_width",
                                max_width=SLIDER_WIDTH)
```

Bounds are in *layout* pixels, so scaling still enlarges them: a bounded slider
is not a small one on a large display.

For the figure, use `artesian.live.responsive(fig)`, which sets
`sizing_mode="scale_width"` and pins the aspect ratio. **Do not use
`stretch_width` for a plot.** It fills the width but pins the height, so the
aspect ratio drifts with the window: a 680x380 plot reads 1.79 as drawn, 2.89 in
a 1100 px column, 4.21 in a 1600 px one. Where meaning lives in a slope -- a
weathering front's depth against time, a fracture-spacing profile -- that is the
same data looking three times gentler to a reader with a wider monitor.

---

## 4. Deploying to GeomorphOnline

Repository: `GeomorphOnline/GeomorphOnline.github.io`, Jekyll with the
minimal-mistakes theme, GitHub Pages **legacy build from `master`, root path**.

### Layout

```
exercises/apps/          compiled apps + THE SHARED WHEELS
_pages/exercises/<name>.md   the human-facing page
```

### Put your app in `exercises/apps/`. This is not optional.

`panel` and `bokeh` are 28.9 MB and 6.1 MB, and they are shared by every app in
that directory. Your app adds roughly 30 KB. A separate directory per exercise
would add another 35 MB each, against GitHub Pages' **1 GB site limit** with the
repository already at ~160 MB.

```sh
artesian build weathering_app.py -o exercises/apps -p /home/awickert/models/corestone \
    -r numpy -r scipy
```

artesian only replaces superseded versions of its own distributions, so this
will not disturb `grlp_panel.*` or the shared wheels. (It used to delete
everything; that was fixed precisely for this.)

### The page

Copy `_pages/exercises/gravel-river-long-profile.md` and change the content. Keep
all of the front matter:

```yaml
---
title: "Exercise: ..."
layout: single
permalink: /exercises/<name>/
author_profile: false
sitemap: false          # keeps it out of sitemap.xml
classes: wide full-bleed
---
```

- `wide` drops the theme's right-hand table-of-contents padding.
- `full-bleed` is a custom class in `assets/css/main.scss` that reclaims the
  left sidebar space. Both are needed; neither alone does the job.
- **`sitemap: false` is obscurity, not privacy.** The repository and the site are
  public. Anyone with the URL can read the page. Do not put solutions there.

Copy the `<iframe>` **and the `<script>` beneath it** verbatim, changing only the
element `id` and the `src`. That script does two jobs: it sizes the frame to its
content (no fixed height can work, since the plot's height follows the reader's
window) and it applies the scaling from section 3.

Inline maths uses **single dollars**: `$D$`, not `$$D$$`. The site's MathJax
config sets `inlineMath: [['$','$'], ['\\(','\\)']]`, so `$$...$$` is display
maths and will break a symbol onto its own centred line.

---

## 5. Provenance: build from a commit, never a dirty tree

**The mistake I made, so you do not.** The first deployment was built by pointing
`artesian build` at `~/models/GRLP` while it sat on a feature branch with
uncommitted edits: 797 lines different from `master`, 51 of them in no commit
anywhere. Students would have been running code that could not be reconstructed
from git.

Nothing revealed it. A modified working tree still reports its version, so the
wheel had the same filename as a clean build -- same name, different bytes.

So, before building:

1. `git -C <source> status` and `git branch --show-current`. Build from a tag or
   a pushed commit, not from whatever is checked out.
2. Record what you built from. `exercises/apps/README.md` holds a provenance
   table -- add a row for your app: model version/commit, application source
   commit, build date, browser requirements.
3. If the app source lives in the model repository, push it before building, so
   the commit you record is fetchable by someone else.

For GRLP the model came from the `v2.1.0` tag, which is also what
`pip install grlp` gives a student who wants to run the same thing at home. That
is a good default if corestone has a release; otherwise a commit on `master`.

---

## 6. Traps, all of which cost me time

1. **`panel convert` can fail and still exit 0.** It prints "Failed to
   convert ... does not publish any Panel contents" and returns success.
   artesian now checks the page exists and raises. If you invoke `panel convert`
   directly, check its output yourself.
2. **Your model must be importable in the *building* environment.**
   `panel convert` executes the app to discover what it serves. Shipping the
   wheel for the browser is not sufficient; `pip install -e .` locally too.
3. **Never name the dict holding your model `state`.** `pn.state` exists; the
   shadowing is silent and confusing. Use `sim`.
4. **Panel's layout fills whatever height the frame gives it.** Measuring
   `body.scrollHeight` on a tall frame reads back the height you last set. The
   deployed script collapses the frame to `0px` before measuring. Adding padding
   to the measurement makes it grow a little on every observer callback.
5. **A zoomed frame reports its own width in local coordinates.** Measure the
   available width from `frame.parentElement`, never from the frame, or the
   factor oscillates.
6. **A zoomed frame needs an explicit pixel width**, not `width: 100%`: a
   percentage resolves against the parent and is *then* scaled by the zoom,
   overflowing the page.
7. **GitHub Pages `builds/latest` returns the previous completed build** while a
   new one is queued. Match the commit SHA, not just `status == "built"`, or you
   will conclude a change did not deploy when it simply has not built yet.
8. **The Pyodide runtime still comes from a CDN** (`cdn.jsdelivr.net`).
   `panel convert` hardcodes it. The demo is not usable fully offline.

Added 2026-09-02 by the corestone build, each one paid for:

9. **Ship artesian itself.** `-p` is repeatable, and section 1's command omits
   it. An app that imports `artesian.live` -- which is what section 2 tells you
   to write -- dies in the browser with `ModuleNotFoundError: No module named
   'artesian'`. Build with `-p <your model> -p /path/to/artesian`. (GRLP does
   not hit this: it predates the helpers and hand-rolls its own toggle.)
10. **The browser's SciPy is newer than yours, and the difference is silent.**
    Pyodide 0.29 ships SciPy 1.14 and NumPy 2.2. `bicgstab`'s relative
    tolerance was `tol` until SciPy 1.12, `rtol` after, and `tol` was REMOVED
    in 1.14 -- so `tol=` works on a workstation with 1.11 and raises
    `TypeError` in the browser. Resolve such names from the signature rather
    than hardcoding them, and *check any keyword argument to SciPy* against
    the version Pyodide ships:
    `curl -s https://cdn.jsdelivr.net/pyodide/v<VER>/full/pyodide-lock.json`.
11. **A worker's exceptions never reach the page console.** In the default
    `pyodide-worker` mode the model runs in a Web Worker, so the traceback
    above appeared NOWHERE: no error in devtools, no error in Selenium's
    browser log. The demo rendered its first frame perfectly and then refused
    to advance, silently. Two things localise this quickly -- run the same app
    under `panel serve` (Python in-process, exceptions visible), and build
    artesian's `examples/hillslope.py` as a known-good control.
12. **Your browser caches the wheel by filename.** A rebuilt wheel keeps its
    name and version, so a reloaded page happily reinstalls the OLD one. I
    spent a round concluding a fix had not worked when it had. Test with a
    fresh profile and caching off, and hard-reload after any rebuild.
13. **Verify by pixels, not by the DOM.** Panel renders into shadow DOM and
    bokeh draws into a canvas, so `document.body.innerText` finds nothing and
    an ordinary Selenium selector finds no buttons. Walk `shadowRoot`
    recursively to click, and decide whether the model is actually running by
    screenshotting before and after and counting changed pixels.

---

## 7. What a reader downloads

For GRLP, measured: ~62 MB on first load -- 11.6 MB Pyodide core, ~16 MB of
numpy/scipy, and 35 MB of self-hosted panel and bokeh wheels. It takes 10-30 s
the first time and is cached afterwards. Say so on the page; the GRLP one does.

A 30-student class is roughly 2 GB against a 100 GB/month soft bandwidth limit.
Comfortable.

`panel` alone is 28.9 MB and is the single largest item. It is already shared, so
there is nothing further to win there.

---

## 8. Checklist

- [ ] `artesian check` on every import your app makes
- [ ] App written, runs natively, `sim` not `state`
- [ ] `responsive(fig)`; `DESIGN_WIDTH` / `SLIDER_WIDTH` declared
- [ ] Source committed and pushed; note the commit
- [ ] Built into `exercises/apps/` (shared wheels; check `grlp_panel.*` survived)
- [ ] Provenance row added to `exercises/apps/README.md`
- [ ] Page created from the GRLP one, front matter kept, script copied
- [ ] Single-dollar inline maths
- [ ] Pushed; Pages build watched **by SHA**; page loaded and the model driven

---

## 9. Known wart

`DESIGN_WIDTH` currently lives in two places: the app declares it, and the page
script repeats it. Nothing enforces that they agree. With a second exercise
arriving, the right fix is for artesian to emit the design width into the page so
the app is the single source of truth. Worth doing before there are three.

## People and pointers

- artesian: <https://github.com/MNiMORPH/artesian> (README covers the build tool,
  the helpers, and the aspect-ratio trap; `docs/validation.md` records how it was
  checked against the hand-rolled original)
- Worked example: `~/models/GRLP/interactive_demo/grlp_panel.py`
- Small example: `examples/hillslope.py` in the artesian repository
- Live result: <https://geomorphonline.github.io/exercises/gravel-river-long-profile/>
