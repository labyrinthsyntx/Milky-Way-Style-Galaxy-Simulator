# GALAXY\_SIM — Milky‑Way‑Style Galaxy Simulator

**Version:** rev‑B (2025‑06‑17)
**Author:** Derek Martinez

---

## 1 . Overview

`GALAXY_SIM.m` is a single‑file MATLAB demo that renders and animates a Milky‑Way‑like spiral galaxy in real time.  It balances visual realism with reasonable performance on a desktop workstation and includes multiple physics modes, interactive controls, and optional video output.

Key features

* **Visual fidelity**: black‑body star colours, magnitude‑scaled marker sizes, thin‑disk dust layer.
* **Disk physics**: analytic rotation curve (pseudo‑isothermal halo), vertical “breathing,” optional rotating central bar.
* **Interactivity**: keyboard camera controls, pause/speed toggle, snapshot hot‑key, live MP4 recording.
* **Physics modes**: `"analytic"` (default, fast), `"gravity"` (direct N‑body, slow), `"barneshut"` (placeholder), `"satellite"` (placeholder for tidal‑interaction demo).

---

## 2 . Requirements

* MATLAB R2019b or newer (tested on R2024a).
* No toolboxes are strictly required, but **Parallel Computing Toolbox** and a CUDA‑capable GPU are recommended for future Barnes–Hut/GPU extensions.

---

## 3 . Quick start

1. Place `GALAXY_SIM.m` in your MATLAB path.
2. `>> run GALAXY_SIM` (or open and press **Run**).
   The galaxy appears face‑on with two spiral arms.
3. Use the keyboard shortcuts below to explore.

### Default hot‑keys

| Key       | Action                                 |
| --------- | -------------------------------------- |
| **← / →** | Rotate camera azimuth                  |
| **↑ / ↓** | Zoom camera in / out                   |
| **Space** | Pause / resume integration             |
| **+ / –** | Speed up / slow down simulation        |
| **S**     | Save a PNG snapshot to current folder  |
| **V**     | Toggle MP4 recording (`GalaxySim.mp4`) |

> **Tip:** Edit parameters under `%% 0. USER SETTINGS` to change star count, arm pitch, dust, etc.  A restart is required after editing.

---

## 4 . Customisation

| Parameter     | Purpose (typical values)                                |
| ------------- | ------------------------------------------------------- |
| `Nstars`      | Total star particles (2 k – 50 k; higher = slower)      |
| `Narms`       | Spiral arm count (2 ≈ Milky Way, 3–4 for grand‑design)  |
| `pitch`       | Arm tightness, rad / kpc (0.2–0.4)                      |
| `orbitModel`  | `"analytic"`, `"gravity"`, `"barneshut"`, `"satellite"` |
| `enableBar`   | `true`/`false` central bar toggle                       |
| `simSpeed`    | Multiplier for time step (0.1 slow‑mo, 1 real‑time)     |
| `recordVideo` | Start recording at launch (`true`)                      |

Performance tips

* Keep `Nstars ≤ 20 000` for smooth 60 fps in analytic mode.
* Use `simSpeed < 0.5` if motion feels too fast.
* `orbitModel = "gravity"` is educational but scales *O(N²)* — keep `Nstars ≤ 1 000`.

---

## 5 . Extending the simulator

Planned but not yet implemented:

1. **Barnes–Hut force solver** (`orbitModel="barneshut"`) for \~50 k stars.
2. **GPU acceleration** using `gpuArray` for pair‑wise gravity.
3. **Satellite galaxy fly‑through** (`orbitModel="satellite"`) to showcase tidal arms.
4. **GUI wrapper** in App Designer for live parameter tuning.

Feel free to fork and contribute!

---

## 6 . License

MIT License — see `LICENSE.txt` (to be added).
