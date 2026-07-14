# edep-Simphony plugin — edep-sim + Simphony GPU Optical Integration

## Project Summary

The **edep-Simphony plugin** (`libedep-simphony-plugin.so`) integrates **edep-sim** (Geant4 CPU simulation) with **Simphony** (GPU optical photon transport via NVIDIA OptiX) into a single pipeline. The result is a simulation where:

- edep-sim handles charged-particle physics (ionisation, EM showers) on CPU
- When Cerenkov/Scintillation photons would be generated, they are collected as "gensteps" and handed to Simphony
- Simphony ray-traces all photons in the event on GPU using OptiX
- The GPU photon hits are written **into the standard `TG4Event` object** under
  `PhotonDetectors["SimphonyPhotonDetector"]`, on equal footing with CPU-tracked
  photons — so any consumer of `TG4Event` (edep-sim CLI ROOT file, Phlex, …)
  sees them without knowing about the plugin. See
  [`doc/photondetectors-conduit.md`](doc/photondetectors-conduit.md) for the
  mechanism.

> **Naming**: Simphony was formerly called *eic-opticks* (upstream is
> [BNLNPPS/simphony](https://github.com/BNLNPPS/simphony)). The local checkout
> and install live under `<prefix>/simphony/`.

---

## Architecture

```
edep-sim process (one event)
  │
  │  ┌─────────────────────────────── CPU (Geant4) ───────────────────────────────┐
  ├─▶│ Geant4 tracks the charged particle step by step                            │
  │  │   Cerenkov/Scintillation fires at each step                                │
  │  │   → instrumented process records a genstep in the SEvt buffer (GPU-side)   │
  │  │   → optical-photon secondaries are KILLED on CPU   (GPU-only, default)     │
  │  │       …or, in DUAL mode, ALSO tracked on CPU by a stock G4Scintillation    │
  │  │         → CPU hits reach the SurfaceSD; per-photon fate + path recorded    │
  │  │ SimphonyStepAction records genstepIdx → G4 TrackID provenance map          │
  │  │ (input-photon mode: primary optical photons are injected into SEvt as      │
  │  │  Opticks input photons instead of relying on gensteps)                     │
  │  └────────────────────────────────────────────────────────────────────────────┘
  │            │  EndOfEventAction calls G4CXOpticks::simulate()  [blocking]
  │            ▼
  │  ┌─────────────────────────── GPU (Simphony / OptiX) ──────────────────────────┐
  │  │ OptiX ray-traces every photon in the event                                 │
  │  │   bounce / scatter / absorb / WLS-reemit until killed or detected          │
  │  │   photon hits a surface with EFFICIENCY>0 → SURFACE_DETECT flag            │
  │  │   max bounces per photon = EDEP_SIMPHONY_MAXBOUNCE                          │
  │  └────────────────────────────────────────────────────────────────────────────┘
  │            │  CPU resumes after the GPU finishes
  │            ▼
  │  ┌─────────────────────────────── CPU (readout) ──────────────────────────────┐
  ├─▶│ reads hits from SEvt (recovers TrackId via genstep provenance map)         │
  │  │   → EDepSim::HitSurface collection "SimphonyPhotonDetector/SimphonyHits"   │
  │  │     picked up by the edep-sim PersistencyManager (ANY backend) and stored  │
  │  │     as TG4Event.PhotonDetectors["SimphonyPhotonDetector"]                  │
  │  │   → also the legacy flat GPUPhotonHits tree   (EDEP_SIMPHONY_LEGACY_HITTREE)│
  │  │ DebugHeavy/input-photon: walk ALL photons + their record buffer            │
  │  │   → GPUPhotonTracks (per-photon fate) + GPUPhotonSteps (full path)         │
  │  └────────────────────────────────────────────────────────────────────────────┘
  │
  └─ edep-sim ROOT file (EDepSimEvents + up to 5 plugin trees)
       EDepSimEvents    ← TG4Event: ionisation, trajectories, primaries, AND
       │                   PhotonDetectors["SimphonyPhotonDetector"] (GPU hits)
       │                   + PhotonDetectors[<gdml-sd-name>] (CPU hits, DUAL)
       GPUPhotonHits    ← GPU detected hits, legacy flat tree   (default on)
       GPUPhotonTracks  ← every GPU photon + final fate             (DebugHeavy)
       GPUPhotonSteps   ← every GPU bounce point (full path)        (DebugHeavy)
       CPUPhotonTracks  ← every CPU photon + final fate             (DUAL)
       CPUPhotonSteps   ← every CPU step point (full path)          (DUAL)
```

> **No ROOT persistency manager?** (e.g. edep-sim embedded in Phlex): the GPU
> hits still reach `TG4Event.PhotonDetectors` — that path needs no ROOT file at
> all. The plugin TTrees then go to a plugin-owned fallback file
> (`EDEP_SIMPHONY_DEBUG_FILE`, default `simphony_debug.root`).

**CPU and GPU are strictly sequential**: the CPU blocks during `simulate()`, then resumes to collect hits. There is no overlap between events. The GPU box runs entirely inside the single `G4CXOpticks::simulate()` call.

### Run modes at a glance

The default pipeline (above) is **GPU-only**: CPU photons are killed, only the GPU
transports light, and only detected hits (`GPUPhotonHits`) are written. Several
runtime modes — all toggled by `EDEP_SIMPHONY_*` environment variables — layer
extra behaviour on top:

| Mode | Env | What it adds |
|---|---|---|
| **GPU-only** (default) | *(none)* | Kill CPU photons; GPU genstep transport → `GPUPhotonHits` |
| **DUAL** | `EDEP_SIMPHONY_DUAL=1` | Also keep CPU optical tracking on, so a stock `G4Scintillation` tracks photons on CPU **alongside** the GPU. CPU and GPU light land in the **same** ROOT file for a per-event comparison |
| **All-photon / trajectory** | `EDEP_SIMPHONY_DEBUGHEAVY=1` | Put Opticks in `DebugHeavy` event mode → save **every** GPU photon (detected or not) + its full bounce-by-bounce path to `GPUPhotonTracks` / `GPUPhotonSteps` |
| **Input-photon** | `EDEP_SIMPHONY_INPUT_PHOTONS=1` | Capture the event's **primary optical photons** and inject them into Opticks as input photons, so the GPU transports the *identical* photons the CPU tracks (clean transport comparison). Implies `DebugHeavy` |

See **GPU Controls** below for the full variable list, and the
[`tests_benchmark` driver](../tests_benchmark/README.md) for ready-made
`run.sh` modes (`both`/`cpu`/`gpu`/`dual`/`photon`/`photon1`/`dualtraj`) that
set these for you.

---

## Key Components

### 1. Simphony (rebuilt from source)
- Source: `<prefix>/simphony/`
- Install: `<prefix>/simphony/install/` (headers under `include/simphony/`)
- Built with **OptiX 8.1.0** headers (ABI 93), compatible with driver 555.42.06
  - OptiX 9.0.0 (ABI 105) is NOT compatible with this driver
  - OptiX 8.1.0 headers live at: `<prefix>/optix810-sdk/`
- PTX kernel: `simphony/install/lib/CSGOptiX7.ptx`

### 2. edep-sim (three source changes)
- Source: `<prefix>/edep-sim/`
- Install: `<prefix>/edep-sim/install/`
- `src/EDepSimUserEventAction.cc` — the external action loop was moved before
  the `if (!HCofEvent) return` early exit, so the GPU plugin runs even on
  events with no ionisation hits.
- `src/EDepSimHitSurface.{hh,cc}` — added a **value constructor**
  (`HitSurface(primaryId, energyDeposit, position, start, pdg, creatorType,
  creatorSubtype)`). The class has private fields and its only filling
  constructor takes a `const G4Step*`; GPU photons have no step, so this is
  what lets the plugin build the hits that reach
  `TG4Event.PhotonDetectors`.
- `src/EDepSimPersistencyManager.cc` — null-guards in
  `SummarizePhotonDetectors` / `SummarizeSegmentDetectors`
  (`if (!g4Hits || …)`). Fixes a segfault when a registered hit-collection
  slot is not filled in a given event (possible for SDs that are not attached
  to a volume, like the plugin's pseudo-SD).

### 3. Plugin library
- Source: this repository
- Build output: `build/libedep-simphony-plugin.so`

| File | Role |
|---|---|
| `src/SimphonyRunAction.cc` | Reads the `EDEP_SIMPHONY_*` knobs and configures `SEventConfig` (DebugHeavy, MaxSlot/Photon/Record/Bounce) **before** `SEvt` is created; sets `KillOpticalPhotons` (off in DUAL); registers `LArTPCSensorIdentifier` **and the `SimphonyPhotonSD` pseudo-SD**; calls `SEvt::CreateOrReuse()` then `G4CXOpticks::SetGeometry()`; acquires the TTree file via `AcquireOutputFile()` (RootPersistencyManager's file when available, else a plugin-owned `EDEP_SIMPHONY_DEBUG_FILE`) and creates the plugin TTrees + their static branch buffers |
| `src/SimphonyEventAction.cc` | Calls `simulate()`; converts detected hits to `EDepSim::HitSurface` and inserts the `SimphonyPhotonDetector/SimphonyHits` collection into the event (`FillPhotonDetectorHits`, **every** event — empty when no hits) so the persistency manager stores them in `TG4Event.PhotonDetectors`; also fills the legacy `GPUPhotonHits` tree; in DebugHeavy/input-photon mode walks **all** GPU photons + their record buffer → `GPUPhotonTracks`/`GPUPhotonSteps`; recovers TrackId via genstep provenance; in input-photon mode captures primary optical photons and injects them as Opticks input photons |
| `src/SimphonyPhotonSD.hh` | Pseudo sensitive detector (never attached to a volume, `ProcessHits` unused): registers the `SimphonyPhotonDetector/SimphonyHits` hit-collection slot that `FillPhotonDetectorHits` fills and `SummarizePhotonDetectors` reads. The SD **name** is the provenance key that separates GPU hits from CPU hits inside `PhotonDetectors` |
| `src/SimphonyStepAction.cc` | Records genstep index → G4 TrackID provenance map; logs DokeBirks dE/dx + photon-yield sampling; fills `CPUPhotonSteps` (every step of every CPU optical photon) |
| `src/SimphonyCpuPhotonTracker.cc/.hh` | External `G4UserTrackingAction` (DUAL/CPU modes): records the **fate** of every CPU-tracked optical photon (detected, absorbed, WLS, escaped, …) → `CPUPhotonTracks`. Reads `G4OpBoundaryProcess::GetStatus()` only at a real geometry boundary to avoid stale `Detection` leaking onto bulk steps |
| `src/SimphonyPhysicsSwap.cc` | Replaces G4Cerenkov/G4Scintillation with instrumented versions; in DUAL mode also installs a stock `G4Scintillation` (+ DokeBirks `AddSaturation`, stacking on) so the CPU produces real trackable photons |
| `src/LArTPCSensorIdentifier.h` | Custom sensor identifier: finds volumes with a G4 SensitiveDetector (works for any geometry, not just PMT-named volumes) |
| `src/plugin_entry.cc` | `extern "C"` factory entry points (RunAction, EventAction, StepAction, PhysicsConstructor, **TrackAction**) |
| `macro/run_3gev_electron.mac` | Run macro (particle gun, geometry update, plugin load) |
| `macro/simphony_plugin.mac` | Loads plugin actions via `/edep/actions/load*` |
| `setup_env.example.sh` | Template for `setup_env.sh` — copy and edit the paths for your machine |
| `setup_env.sh` (gitignored) | User-local copy of the template; sets all required environment variables |

> **Geometry**: the plugin is geometry-agnostic — point `-g` at any GDML whose
> LAr has proper optical material properties (RINDEX etc.) and whose photon
> sensors carry a G4 SensitiveDetector. `LArTPCSensorIdentifier` finds the
> sensor volumes automatically (no PMT naming required).

---

## Environment Setup

### Dependencies

The build finds all of its dependencies through the standard CMake
`find_package` mechanism and `CMAKE_PREFIX_PATH` shell environment variable.  In
cases where dependencies can not be found in this way, `cmake` can be given
command line options to help in their location.

Direct dependencies, each located via its own exported CMake config:

| Package | Provides | Config directory (`<prefix>` = its install prefix) |
|---|---|---|
| `ROOT`     | Core, Tree, RIO                | `<prefix>/share/root/cmake/` |
| `Geant4`   | Geant4 libraries               | `<prefix>/lib/cmake/Geant4/` |
| `EDepSim`  | `EDepSim::edepsim`, `edepsim_io` | `<prefix>/lib/cmake/EDepSim/` |
| `glm`      | `glm::glm`                     | `<prefix>/share/glm/` |
| `simphony` | `simphony::{G4CX,U4,CSGOptiX,QUDARap,CSG,SysRap}` | `<prefix>/lib/cmake/simphony/` |

The `simphony::*` imported targets carry their own include directories and pull
in **plog, OptiX and the CUDA runtime** transitively, so the plugin never names
those paths itself.  Simphony 0.8.0 is expected include commit 042a282 and then
"gml" becomes a transitive dependency.

The **C++ standard** is inherited from ROOT's build unless you explicitly pass
`-DCMAKE_CXX_STANDARD=NN` to force a value.

### Building the plugin

First activate the environment that provides the dependencies, e.g. the Spack
view used to build them:

```bash
source <prefix>/.envrc         
## or: 
# cd <prefix> && direnv allow
```

**Minimal build config** — every dependency is found automatically because the
activated environment put its prefixes on `CMAKE_PREFIX_PATH`, and the C++
standard is taken from ROOT:

```bash
cmake -S . -B build
cmake --build build -j
```

**Maximal build config** — nothing is on `CMAKE_PREFIX_PATH`, so each dependency
is pointed at explicitly and the C++ standard is pinned. Each `*_DIR` is the
directory holding that package's `*Config.cmake` (see the table above):

```bash
cmake -S . -B build \
  -DCMAKE_CXX_STANDARD=23 \
  -DROOT_DIR=/opt/root/share/root/cmake \
  -DGeant4_DIR=/opt/geant4/lib/cmake/Geant4 \
  -DEDepSim_DIR=/opt/edep-sim/lib/cmake/EDepSim \
  -Dglm_DIR=/opt/glm/share/glm \
  -Dsimphony_DIR=/opt/simphony/lib/cmake/simphony \
  -DCUDAToolkit_ROOT=/usr/local/cuda
cmake --build build -j
```

> Simphony's transitive dependencies (CUDA, OptiX, plog, …) are normally resolved
> from the same prefixes. If one lives somewhere unusual, hint it the same way —
> e.g. `-DCUDAToolkit_ROOT=…` — or just prepend its prefix to
> `-DCMAKE_PREFIX_PATH="/a;/b;/c"`.

### Runtime environment

Running the plugin needs three *distinct kinds* of shell variable, kept separate
below.

**1. Placeholders** — local shell variables that are **not** part of any
interface. Nothing reads them by these names; they exist only to build the real
variables in groups 2 and 3. Edit them to match your machine.

```bash
PREFIX=/path/to/deps                 # spack view (or dependency install) prefix
SIMPHONY_INST=${PREFIX}              # simphony prefix (holds lib/CSGOptiX7.ptx)
PLUGIN_BUILD=/path/to/this/repo/build
```

**2. Standard PATH-like variables** — the loader search path understood by the
OS. Order matters: the correctly-versioned libraries must come first.

```bash
export LD_LIBRARY_PATH=${PLUGIN_BUILD}:${PREFIX}/lib:${PREFIX}/lib64:${LD_LIBRARY_PATH}
# add the CUDA runtime libdir here too if it is not already under ${PREFIX}
```

> **Warning**: do NOT add any directory containing old, conflicting copies of
> Geant4, edep-sim or Simphony to `LD_LIBRARY_PATH` — they will shadow the
> correct libraries and cause hard-to-diagnose symbol errors.

**3. Plugin run-time options** — the actual runtime interface, read by the plugin
and Simphony to select behavior:

```bash
# Shared library edep-sim loads, plus the physics-constructor hook.
export PLUGIN_LIB=${PLUGIN_BUILD}/libedep-simphony-plugin.so
export EXTRAPHYSICS="EXTERN:${PLUGIN_LIB}:CreatePhysicsConstructor"

# Opticks / Simphony integration
export OPTICKS_INTEGRATION_MODE=1                        # 1 = GPU-only
export CSGOptiX__ptxpath=${SIMPHONY_INST}/lib/CSGOptiX7.ptx
export OPTICKS_MAX_SLOT=M1                               # GPU photon-slot cap (M1 = 1e6)
export OPTICKS_OUT_FOLD=/path/to/scratch/opticks_output  # Simphony .npy output folder

# Plugin knob (see "GPU Controls" below for the full EDEP_SIMPHONY_* set)
export EDEPSIM_DOKEBIRKS_VISE=1                          # drift-field-aware photon yield
```

`macro/simphony_plugin.mac` refers to `$(PLUGIN_LIB)`, so that variable must be
exported before the run. See **GPU Controls** for every `EDEP_SIMPHONY_*` knob.

### Run

```bash
edep-sim -p QGSP_BERT \
         -g <path-to-geometry>.gdml \
         -o output.root \
         -e 1000 \
         macro/run_3gev_electron.mac
```

The particle gun is configured in `run_3gev_electron.mac`. Edit `/gps/energy` and
`/gps/position` there to change the beam.

> **Note**: despite its name, `run_3gev_electron.mac` currently fires a **1 MeV**
> electron (`/gps/energy 1 MeV`, `/gps/position 0 0 -400 mm`). Change `/gps/energy`
> if you want a different beam energy.

---

## GPU Controls (`EDEP_SIMPHONY_*` environment variables)

These are read by the plugin at runtime (in `SimphonyRunAction::BeginOfRunAction`,
before `SEvt` is created). All are optional — unset means the default GPU-only,
hit-only path.

| Variable | Meaning | Default |
|---|---|---|
| `EDEP_SIMPHONY_DUAL` | `1` = also track optical photons on **CPU** (stock `G4Scintillation`) next to the GPU genstep path, so both write to the same ROOT file. Keeps `KillOpticalPhotons` **off** | unset (GPU-only, CPU photons killed) |
| `EDEP_SIMPHONY_DEBUGHEAVY` | `1` = Opticks `DebugHeavy` event mode: gather photon + record + hit, so **all** GPU photons + full per-photon trajectory are saved (`GPUPhotonTracks`/`GPUPhotonSteps`). Memory-hungry — few-photon runs only | unset (hit-only `Minimal` mode) |
| `EDEP_SIMPHONY_INPUT_PHOTONS` | `1` = inject the event's **primary optical photons** into Opticks as input photons (instead of relying on scintillation gensteps). Implies DebugHeavy | unset |
| `EDEP_SIMPHONY_MAXBOUNCE` | Max GPU bounces per photon before it is killed (the transport bounce cap). Applies in **all** modes. 128 nm photons Rayleigh-scatter heavily in LAr and need a high budget; too low → photons hit the cap before being absorbed/detected | Opticks default (~31) |
| `EDEP_SIMPHONY_MAXRECORD` | Length of the saved per-photon trajectory buffer. Clamps `MaxBounce` to `MaxRecord-1` unless `MAXBOUNCE` is set after | auto (DebugHeavy) |
| `EDEP_SIMPHONY_MAXSLOT` | GPU photon-slot / photon budget for DebugHeavy (record alloc = MaxSlot × MaxRecord). Capped small for debug runs to avoid CUDA OOM | `200000` (DebugHeavy) |
| `EDEP_SIMPHONY_CERENKOV` | `0` = disable Cerenkov (scintillation-only comparison) | — |
| `EDEP_SIMPHONY_SCINT` | Scint process: `thin` (`SimphonyScintProcess`, default) or `fork` (`Local_DsG4Scintillation`) | `thin` |
| `EDEP_SIMPHONY_LEGACY_HITTREE` | `0` = do **not** write the legacy flat `GPUPhotonHits` tree. The same hits are always available in `TG4Event.PhotonDetectors["SimphonyPhotonDetector"]`; the tree is kept while old analysis scripts migrate | `1` (tree written) |
| `EDEP_SIMPHONY_DEBUG_FILE` | Path for the plugin TTrees when the edep-sim ROOT persistency manager is **not** available (e.g. under Phlex). Ignored in normal CLI runs, where the trees go into the edep-sim output file | `simphony_debug.root` |
| `EDEP_SIMPHONY_FORCE_NO_ROOTPM` | `1` = pretend the ROOT persistency manager is absent (test hook for the fallback path above) | unset |

> **Record-length hard cap**: the per-photon GPU *flag-history* record is also
> bounded by Opticks `sseq::SLOTS` (= `16 × NSEQ`, baked into
> `simphony/sysrap/sseq.h`; currently **112** with `NSEQ = 7`). To record paths
> longer than that you must raise **both** `EDEP_SIMPHONY_MAXBOUNCE` (env) **and**
> `NSEQ` (recompile simphony). Note `MAXBOUNCE` is the *transport* kill — a photon
> may be allowed to keep bouncing well past the record length (only the first
> `SLOTS` points are stored). The CPU side has **no** bounce cap — it propagates
> until physically absorbed/detected.

### Loading the CPU tracking action (DUAL / CPU fate)

To fill `CPUPhotonTracks` (per-photon CPU fate), the tracking-action factory must
be loaded in the macro alongside the other plugin actions:

```
/edep/actions/loadUserTrackAction $(PLUGIN_LIB) CreateUserTrackAction ""
```

`CPUPhotonSteps` (full CPU path) is filled by the step action and additionally
needs edep-sim's own trajectory saving turned on in the macro **before**
`/edep/update`:

```
/edep/db/set/savePhotonTraj  true
/edep/db/set/saveAllPrimTraj true
```

---

## Reading the Output ROOT File

### GPU hits inside TG4Event (the primary record)

The detected GPU hits are stored **inside each event** as `TG4PhotonHit`
objects, keyed by sensitive-detector name — GPU photons under
`SimphonyPhotonDetector`, CPU-tracked photons (DUAL mode) under the
geometry's own SD name. This is the persistency-agnostic record: it exists in
any context that stores `TG4Event`, with or without a ROOT file.

```python
import ROOT
f = ROOT.TFile.Open("output.root")
t = f.Get("EDepSimEvents")
for i in range(t.GetEntries()):
    t.GetEntry(i)
    for sd, hits in t.Event.PhotonDetectors:   # sd = SD name string
        for h in hits:
            h.GetStop()        # detection position+time (TLorentzVector; mm, ns)
            h.GetStart()       # creation point (zero unless DebugHeavy run)
            h.GetWavelength()  # detected (post-WLS) wavelength
            h.GetPrimaryId()   # stored-trajectory id of the charged parent
            h.GetProcess()     # creator subtype: 21 Cerenkov / 22 Scint / 34 WLS
```

`tests_benchmark/check_photondet_migration.py` verifies these hits are 1:1
identical with the legacy `GPUPhotonHits` tree.

### Plugin TTrees

Depending on the run mode, the file additionally holds up to five plugin
TTrees alongside edep-sim's own `EDepSimEvents` (in a fallback
`simphony_debug.root` when there is no ROOT persistency manager):

| Tree | Written when | One row per | Contents |
|---|---|---|---|
| `GPUPhotonHits` | `EDEP_SIMPHONY_LEGACY_HITTREE=1` (default) | detected GPU photon | **legacy** flat copy of the PhotonDetectors hits, kept while analysis scripts migrate |
| `GPUPhotonTracks` | DebugHeavy / input-photon | every GPU photon | per-photon final **fate** (detected or not + why) |
| `GPUPhotonSteps` | DebugHeavy / input-photon | every GPU step point | full bounce-by-bounce GPU path |
| `CPUPhotonTracks` | DUAL + `CreateUserTrackAction` | every CPU optical photon | per-photon CPU **fate** |
| `CPUPhotonSteps` | DUAL + photon-traj saving | every CPU step point | full CPU path |

> Full per-photon **trajectories** live only in the `*PhotonSteps` trees —
> `TG4Event` cannot hold them (`TG4PhotonHit` is start/stop only, and edep-sim's
> internal trajectory store sparsifies photon paths to ~2 points by design).

```python
import uproot
import numpy as np

f = uproot.open("output.root")

# GPU optical photon hits (legacy flat tree; same content as PhotonDetectors)
hits = f["GPUPhotonHits"]
print(hits.keys())
# ['EventId', 'TrackId', 'Process', 'Wavelength', 'HitPos', 'StartPos']

df = hits.arrays(library="pd")
print(df.head())

# All-photon fate + trajectory (DebugHeavy / input-photon modes)
tracks = f["GPUPhotonTracks"].arrays(library="pd")  # one row per photon + fate
steps  = f["GPUPhotonSteps"].arrays(library="pd")   # full GPU bounce path

# CPU side (DUAL mode)
cpu_tracks = f["CPUPhotonTracks"].arrays(library="pd")
cpu_steps  = f["CPUPhotonSteps"].arrays(library="pd")

# CPU ionisation + trajectory info
events = f["EDepSimEvents"]
```

With ROOT/PyROOT:

```python
import ROOT
f = ROOT.TFile("output.root")
t = f.Get("GPUPhotonHits")
for entry in t:
    print(entry.EventId, entry.TrackId, entry.Wavelength,
          entry.HitPos.X(), entry.HitPos.Y(), entry.HitPos.Z())
```

**`GPUPhotonHits` branches:**

| Branch | Type | Description |
|---|---|---|
| `EventId` | int | Geant4 event number |
| `TrackId` | int | G4 TrackID of the charged parent (1 = primary) |
| `Process` | int | 0 = Cerenkov, 2 = Scintillation, -1 = unknown |
| `Wavelength` | float | Photon wavelength in nm |
| `HitPos` | TLorentzVector | Hit position (x,y,z in mm; t in ns) |
| `StartPos` | TLorentzVector | Genstep origin (reserved; currently zero-filled) |

**`GPUPhotonTracks` branches** (one row per GPU photon, detected or not):

| Branch | Type | Description |
|---|---|---|
| `EventId` | int | Geant4 event number |
| `PhotonId` | int | Photon index within the event |
| `TrackId` | int | GPU photon global index (one per photon; persists across WLS re-emission) |
| `ParentId` | int | Charged G4 parent track of the genstep, or **-1** for input/primary photons |
| `Process` | int | 0 = Cerenkov, 2 = Scintillation, 3 = TORCH (primary), -1 |
| `Wavelength` | float | nm |
| `Flag` | unsigned | Terminating Opticks flag bit |
| `FlagMask` | unsigned | OR of all history flags |
| `Detected` | int | 1 if SURFACE_DETECT / EFFICIENCY_COLLECT |
| `Fate` | char[8] | Short reason: `SD`, `AB`, `SA`, `MI`, … |
| `NStep` | int | Number of recorded trajectory points |
| `StartPos` / `EndPos` | TLorentzVector | First / last record point (mm, ns) |

**`GPUPhotonSteps` branches** (one row per GPU step point):

| Branch | Type | Description |
|---|---|---|
| `EventId`, `PhotonId` | int | Identify the photon |
| `Step` | int | Step index along the photon |
| `Flag` | unsigned | Opticks flag at this step |
| `FlagAbbr` | char[8] | Abbreviated flag |
| `Pos` | TLorentzVector | Step position (mm, ns) |

**`CPUPhotonTracks` branches** (one row per CPU-tracked optical photon):

| Branch | Type | Description |
|---|---|---|
| `EventId`, `TrackId`, `ParentId` | int | G4 identifiers |
| `Wavelength` | float | nm |
| `Energy` | float | eV |
| `Detected` | int | 1 if boundary status was `Detection` |
| `Fate` | char[16] | `SD`, `SA`, `AB`, `WLS`, `MI`, `Reflect`, `Rayleigh`, or the process name |
| `EndVolume` | char[64] | Volume where the photon ended (`OutOfWorld` if it escaped) |
| `NStep` | int | Number of geometry steps |
| `StartPos` / `EndPos` | TLorentzVector | Creation / termination point (mm, ns) |

**`CPUPhotonSteps` branches** (one row per CPU optical-photon step):

| Branch | Type | Description |
|---|---|---|
| `EventId`, `TrackId` | int | G4 identifiers |
| `Step` | int | G4 step number along the photon |
| `Process` | char[24] | Process that ended the step |
| `Volume` | char[48] | Post-step volume |
| `Pos` | TLorentzVector | Post-step position (mm, ns) |

> **Fate codes** (shared CPU/GPU vocabulary): `SD` detected on a sensor ·
> `SA` absorbed at a surface · `AB` bulk absorption in LAr · `WLS` absorbed by a
> wavelength shifter · `MI` left the world (escaped) · `Reflect`/`Rayleigh` ended
> on a reflection/scatter.

---

## Relevant Source

| Item | Path |
|---|---|
| Plugin source | this repository |
| Plugin build | `build/` |
| PhotonDetectors conduit doc | [`doc/photondetectors-conduit.md`](doc/photondetectors-conduit.md) — SEvt→HitSurface field mapping, required edep-sim patches, rationale |
| Simphony repo | [github.com/Ningclover/simphony](https://github.com/Ningclover/simphony) — source `<prefix>/simphony/`, install `<prefix>/simphony/install/` |
| edep-sim repo | [github.com/ClarkMcGrew/edep-sim](https://github.com/ClarkMcGrew/edep-sim) — source `<prefix>/edep-sim/`, install `<prefix>/edep-sim/install/` |
| gegede (geometry tool) | [github.com/brettviren/gegede](https://github.com/brettviren/gegede) |
| Benchmark CPU-vs-GPU driver | [github.com/Ningclover/tests_benchmark](https://github.com/Ningclover/tests_benchmark) — example geometry, `.mac` files, and visual tools that drive this plugin |
