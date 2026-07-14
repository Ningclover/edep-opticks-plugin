# Update: Simphony GPU photon hits now land in TG4Event.PhotonDetectors

*Xuyang Ning, 2026-07-14.  Follow-up to Brett's 7/8 email ("plugin writes
Simphony photons directly to the TFile") and Clark's reply (HitSurface SD →
PhotonDetectors, different SD names per source).  Implemented, validated, and
pushed — details below.*

## TL;DR

The plugin no longer needs `EDepSim::RootPersistencyManager` for its photon
output.  GPU hits are now converted to `EDepSim::HitSurface` objects and
inserted as a normal Geant4 hit collection, so the edep-sim
`PersistencyManager` — **any** concrete backend, including the Phlex thin
one — automatically stores them as

```
TG4Event.PhotonDetectors["SimphonyPhotonDetector"]     (GPU / Simphony photons)
TG4Event.PhotonDetectors["<gdml-sd-name>"]             (CPU / Geant4 photons)
```

Different SD names per source, exactly as proposed in the email thread.  No
TFile, no downcast, no Phlex-specific code.

Commits:
- **edep-simphony-plugin** `482905d` on `master`
  (github.com/Ningclover/edep-simphony-plugin)
- **edep-sim** `716a1ee` on `integrate_eicoptics`
  (github.com/Ningclover/edep-sim) — two small patches, PR-able to upstream

## The problem (recap)

The plugin used to persist ALL its output by downcasting the persistency
manager to `EDepSim::RootPersistencyManager`, taking its `TFile`, and writing
private TTrees next to `EDepSimEvents`.  Under Phlex the concrete manager is a
thin TG4Event-only class → the downcast returns null → the plugin bailed out
and every GPU photon was silently lost.  Even in CLI mode, GPU hits lived in a
side tree, never inside `TG4Event`.

## How the fix works

`EDepSim::PersistencyManager::SummarizePhotonDetectors` already walks **every**
hit collection registered in the event and accepts any whose hits are
`EDepSim::HitSurface`, keyed by SD name.  That is the standard doorway into
`PhotonDetectors` — so the fix makes the GPU look like one more sensitive
detector:

1. **A pseudo-SD** (`SimphonyPhotonSD`, name `SimphonyPhotonDetector`,
   collection `SimphonyHits`) is registered with `G4SDManager` at begin-of-run.
   It is never attached to a volume (GPU photons don't step); it exists so the
   hit-collection table has a slot.
2. **At end of event**, after `G4CXOpticks::simulate()`, the plugin converts
   each SEvt hit (an `sphoton`) into a `HitSurface` and inserts the collection
   into `G4HCofThisEvent`.  The persistency manager picks it up when it stores
   the event (user actions run before `Store`, so ordering just works).
3. Field mapping (semantics match the CPU-tracked photons):
   - position/time: verbatim (Opticks works in G4 units, mm/ns);
   - energy: `h*c / wavelength` with Geant4's own constants, so
     `TG4PhotonHit::GetWavelength()` round-trips exactly;
   - `PrimaryId`: the charged parent G4 track, recovered from the genstep
     provenance map (photon index → genstep offset range → recorded TrackID);
     the persistency manager then remaps it onto the stored trajectory, same
     as for CPU photons;
   - creator process: decoded from the photon's Opticks flag history, and
     WLS-aware — a re-emitted photon (`BULK_REEMIT`) reports `OpWLS` (34),
     matching Geant4's convention that the arriving photon is a new track
     created by the WLS process; else Cerenkov (21) / Scintillation (22);
   - `Start`: the true creation point when the record buffer is available
     (debug runs); otherwise zero — the case `TG4PhotonHit` already documents
     for offloaded photon tracking.

## What edep-sim needed (2 small patches, on `integrate_eicoptics`)

1. **`HitSurface` value constructor** (`EDepSimHitSurface.{hh,cc}`): the class
   has private fields and its only filling constructor takes a `const
   G4Step*`; GPU photons have no step (and no live track, no creator-process
   instance).  A 9-line constructor taking the seven values directly.  No
   behavior change for existing code; the persistent `TG4PhotonHit` / ROOT
   I/O layer is untouched.
2. **Null-guards in the summarizers** (`EDepSimPersistencyManager.cc`,
   `SummarizePhotonDetectors` + `SummarizeSegmentDetectors`):
   `HCofEvent->GetHC(id)` returns null for a registered-but-unfilled slot,
   and the old code dereferenced it unconditionally.  Volume-attached SDs
   never hit this (their `Initialize()` always inserts a collection), but a
   detached SD does — we found it as a real segfault on an event with zero
   GPU hits.  The plugin also now inserts an *empty* collection on every
   event, so both sides are safe independently.

These two are the PR-able unit for upstream (they sit on top of the earlier
`268ddde` "external actions before the no-hits early return" fix, which the
plugin also needs).

## Complete edep-sim patch set (branch `integrate_eicoptics`, based on Clark's `master`)

Everything the plugin needs from edep-sim is contained in **two commits, six
file changes**, on `integrate_eicoptics` at
github.com/Ningclover/edep-sim (HEAD = `716a1ee`).

Base: the branch is on the **ClarkMcGrew/edep-sim `master` lineage** (NOT the
DUNE/edep-sim fork, which has diverged separately). **Update (2026-07-14
evening):** Clark's current master (`4548701`, v4.3.0, incl. his
trajectory-rule rework) has been **merged into the branch** (merge commit
`4b2a99b`, conflict-free — all three of our patches survived intact).
edep-sim and the plugin were rebuilt against v4.3.0 and the full validation
suite re-passed with identical deterministic results (10-event: 5518/2219,
PhotonDetectors ≡ legacy 1:1; dualtraj: 575 CPU + 408 GPU). The branch is
now directly PR-able against Clark's master with no rebase needed.

Checked against Clark's `4548701` (2026-07-14): **neither patch is upstream
yet** — `HitSurface` still has only the default and `G4Step` constructors,
and both summarizers still dereference `GetHC()` unguarded. His PR #95
(`bea3aae`, 7/7) is the *complement*, not a duplicate: it adds public read
accessors (`GetEventSummary()`, `GetPhotonDetectors()`, …) for a derived thin
persistency manager to hand `TG4Event` content out. So: his PR = the read
side, our two patches = the write side (getting offloaded hits *in*).

### Commit `268ddde` (2026-04-28) — "Fix optical photon integration for GPU plugin support"

| File | Change | Why necessary |
|---|---|---|
| `src/EDepSimUserEventAction.cc` | call the **external user actions before** the `if (!HCofEvent) return` early exit | the plugin's `EndOfEventAction` is where the GPU launch happens and (now) where the PhotonDetectors collection is inserted — it must run on **every** event, including events with no ionisation hits. Without this, such events lose their GPU photons and (post-migration) leave the hit-collection slot unfilled |
| `src/EDepSimSurfaceSD.cc` | null-check the stacking action before calling `SetKillOpticalPhotons` | `SurfaceSD` assumed edep-sim's own stacking action is always installed; in the plugin's setup that assumption can fail → crash guard |
| `src/EDepSimPersistencyManager.cc` | use `CreateAttValues()` instead of the deprecated `GetAttValues()` | fixes the trajectory-error diagnostic path that plugin-era running exercised |

### Commit `716a1ee` (2026-07-14, latest) — "Support externally filled photon hits (GPU-offloaded transport)"

| File | Change | Why necessary |
|---|---|---|
| `src/EDepSimHitSurface.hh/.cc` | new **value constructor** `HitSurface(primaryId, energyDeposit, position, start, pdg, creatorType, creatorSubtype)` | `HitSurface` is the **only** hit class `SummarizePhotonDetectors` accepts into `PhotonDetectors`, its fields are private with no setters, and its only filling constructor needs a `const G4Step*` — which offloaded (GPU) photons don't have. Pure initializer list; zero impact on existing callers and on the persistent `TG4PhotonHit`/ROOT I/O layer |
| `src/EDepSimPersistencyManager.cc` | null-guards `if (!g4Hits \|\| g4Hits->GetSize()<1)` in `SummarizePhotonDetectors` **and** `SummarizeSegmentDetectors` | the `G4HCtable` lists every collection that *could* exist, but `GetHC(id)` returns **null** for a slot not filled this event. Volume-attached SDs always fill theirs in `Initialize()`; a detached SD (like the plugin's pseudo-SD) doesn't get `Initialize()` called — observed as a real segfault on an event with zero GPU hits. Makes a missing collection equivalent to an empty one |

Everything else in our edep-sim working tree (DokeBirks debug prints,
`.gitignore`, tutorial/test files) is debug or local housekeeping —
deliberately **not** committed, so the branch is exactly the upstream-PR
unit.

## Debug output (full trajectories) — kept, made portable

We checked whether the full photon paths could move into `TG4Event` too: they
can't, structurally.  edep-sim's trajectory store collapses any photon with
<1 eV SD deposit to 2 points *before* any accuracy/rule option applies,
excludes optical steps from point selection, and caps the densifier — and GPU
photons are never G4 tracks at all.  (Measured on a 1 MeV e⁻: internal store
median 2 points/photon vs the plugin's full-path record median 8 / mean 118 /
max 50k.)

So the full-trajectory debug trees (`CPU/GPUPhotonTracks`, `CPU/GPUPhotonSteps`)
stay plugin-owned; what changed is where they go:

- ROOT persistency manager present (edep-sim CLI): same file as always,
  next to `EDepSimEvents` — zero change for existing workflows;
- otherwise (Phlex): a plugin-owned sidecar file,
  `EDEP_SIMPHONY_DEBUG_FILE` (default `simphony_debug.root`).

The legacy flat `GPUPhotonHits` tree is still written by default
(`EDEP_SIMPHONY_LEGACY_HITTREE=0` disables it) until our analysis scripts
finish migrating to `PhotonDetectors`.

## Validation (all pass, deterministic seeds)

| Check | Result |
|---|---|
| 10-event benchmark, GPU mode | physics identical to pre-change baseline (CPU 5518 / GPU 2219 hits); `PhotonDetectors` ≡ legacy tree **1:1** (count, position, time, wavelength) |
| DUAL run (1 MeV e⁻) | ONE `TG4Event` carries both populations: `PhotonDetector` (CPU, 575) + `SimphonyPhotonDetector` (GPU, 408); creator-process codes agree (34 = OpWLS for the TPB-shifted detected photons on both sides) |
| Pure-photon regression (1000 × 128 nm) | bit-for-bit with the June result: CPU 39/1000, GPU 36/1000 detected |
| Full-trajectory debug | `compare_traj_methods.py` unchanged |
| No-ROOT-manager emulation (`EDEP_SIMPHONY_FORCE_NO_ROOTPM=1`) | `PhotonDetectors` still filled in the event output; all five plugin trees land in `simphony_debug.root`; no crash |

Reading the result needs nothing new:

```python
f = ROOT.TFile.Open("output.root"); t = f.Get("EDepSimEvents")
t.GetEntry(0)
for sd, hits in t.Event.PhotonDetectors:
    print(sd, len(hits))     # 'SimphonyPhotonDetector' 236   (+ CPU SD in DUAL)
```

## What this means for Phlex

- Nothing to implement on the Phlex side for the hits: if the thin persistency
  manager reuses `UpdateSummaries` (or calls `SummarizePhotonDetectors`), the
  GPU photons appear in `TG4Event.PhotonDetectors["SimphonyPhotonDetector"]`
  automatically.
- Build order: edep-sim needs the two-patch branch (`716a1ee`), then rebuild
  the plugin (`482905d`).
- The SD name string is a single constant (`SimphonyPhotonSD::kSDName`) — easy
  to change if a different convention is preferred.
- Full docs: `doc/photondetectors-conduit.md` in the plugin repo (field
  mapping, the empty-collection invariant, every env switch).

## Remaining (phase 2, our side)

- Migrate the analysis scripts off the legacy `GPUPhotonHits` tree
  (analyze.py, compare_photon.py, extract_boundaries.py; the viz extractor
  turned out not to use it at all), then flip the legacy default off.
- Open the upstream edep-sim PR from `integrate_eicoptics`.
