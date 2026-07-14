# GPU photon hits → `TG4Event.PhotonDetectors`: the HitSurface conduit

*2026-07-11.  Covers `SimphonyEventAction::FillPhotonDetectorHits()`,
`SimphonyPhotonSD`, and the two small edep-sim patches this feature requires.*

## 1. Why this exists

Historically the plugin persisted GPU hits by downcasting the active
persistency manager to `EDepSim::RootPersistencyManager`, grabbing its
`TFile`, and writing a private `GPUPhotonHits` TTree next to `EDepSimEvents`.
That breaks in any context where the concrete persistency manager is not the
ROOT one (e.g. Phlex embeds edep-sim with a thin TG4Event-only manager: the
downcast returns null and all GPU output used to be silently dropped). It also
left GPU photons *outside* `TG4Event`, on unequal footing with CPU-tracked
photons.

The fix routes GPU hits through edep-sim's standard mechanism instead, so
they end up **inside the event object**:

```
TG4Event.PhotonDetectors["SimphonyPhotonDetector"]   (GPU / Simphony photons)
TG4Event.PhotonDetectors["<gdml-sd-name>"]           (CPU / Geant4 photons)
```

Distinct sensitive-detector names keep the two populations distinguishable
when both run in one job (DUAL mode).

## 2. The conduit: how hits get into TG4Event at all

`EDepSim::PersistencyManager::UpdateSummaries` fills `TG4Event` from exactly
four sources; there is no generic "append to event" API. Photon hits have one
doorway:

- `SummarizePhotonDetectors` (EDepSimPersistencyManager.cc) walks **every**
  hit collection registered in the event's `G4HCofThisEvent`.
- For each collection it `dynamic_cast`s **hit 0** to `EDepSim::HitSurface`.
  Match → each hit is copied into a `TG4PhotonHit` under
  `PhotonDetectors[SDname]`. No match → the collection is skipped.

So the hit **class is the routing key**: to reach `PhotonDetectors` you must
produce `EDepSim::HitSurface` objects in a registered hit collection —
nothing else is accepted. This runs in the *base* persistency manager, so any
backend that reuses `UpdateSummaries` (CLI ROOT manager, Phlex thin manager)
gets the GPU photons with no code change.

Ordering works without any hook: the plugin inserts the collection during its
`EndOfEventAction`, and `PersistencyManager::Store` runs afterwards inside
`G4RunManager::AnalyzeEvent`.

### The pseudo-SD (`SimphonyPhotonSD.hh`)

A hit collection needs a slot in the global `G4HCtable`, and slots come from
registered sensitive detectors. `SimphonyPhotonSD` exists solely for that:

- SD name `SimphonyPhotonDetector`, collection name `SimphonyHits`;
- registered with `G4SDManager` in `SimphonyRunAction::BeginOfRunAction`;
- **never attached to a logical volume** — GPU photons do not step, so
  `ProcessHits()` is never called and just returns false.

Being detached has one important consequence (§5).

## 3. `FillPhotonDetectorHits()` field by field

Input: after `gx->simulate()`, the SEvt hit array holds one `sphoton` per
detected photon (terminating flag `SURFACE_DETECT`) — a 4×4-float snapshot
*at the moment of detection*: position (mm), time (ns), wavelength (nm),
photon slot `index`, terminating `flag`, and `flagmask` = OR of every flag in
the photon's history. No G4 track or step exists for these photons.

Output: one `EDepSim::HitSurface` per hit, built with the value constructor
(§6), semantically equivalent to what `HitSurface(const G4Step*)` produces
for a CPU photon.

| HitSurface field | GPU source | CPU equivalent / convention |
|---|---|---|
| `fPosition` | `(p.pos, p.time)` verbatim — Opticks works in G4 units (mm, ns) | post-step point at detection |
| `fEnergyDeposit` | `h_Planck*c_light / (p.wavelength*nm)` | `step->GetTotalEnergyDeposit()` = full photon energy at a photocathode |
| `fPrimaryId` | genstep walk, see §3.1 | `track->GetParentID()` = charged emitter track |
| `fStart` | record buffer row `index`, entry 0 (DebugHeavy runs); else `(0,0,0,0)` | track vertex; TG4PhotonHit documents the zero case as expected for offloaded tracking |
| `fPDGEncoding` | `G4OpticalPhoton::OpticalPhotonDefinition()->GetPDGEncoding()` | same (−22 in this Geant4) |
| `fCreatorType/Subtype` | decoded from `flagmask`, see §3.2 | `track->GetCreatorProcess()` type/subtype |

Notes on individual fields:

- **Energy**: the sphoton stores wavelength, not energy, so we invert with
  Geant4's own `h_Planck`, `c_light`, `nm` constants. `TG4PhotonHit`
  re-derives wavelength as `1.23984e-9 / EnergyDeposit` (hc in MeV·mm) —
  same constants both directions, so the round-trip is exact (validated 1:1
  against the legacy tree). This is the *detected* wavelength: a TPB-shifted
  photon reports ~450 nm, exactly as the CPU sensor would see it.
- **Start**: DebugHeavy runs gather the `record` array shaped
  `(n_photon, max_record, 4, 4)` — a full sphoton snapshot per bounce. Row
  `p.get_index()`, entry 0 is the creation point; `flag()!=0` distinguishes a
  written slot from zero padding. Production (hit-only) runs don't gather it,
  and Start legitimately stays zero.

### 3.1 PrimaryId: the two-hop genstep walk

The parent linkage is recorded at *emission* time, before photons ever reach
the GPU. The instrumented scintillation process emits one **genstep** per
charged step ("G4 track T emits N photons here"), and the GPU materializes
photons contiguously per genstep: genstep 0 → indices `[0,N₀)`, genstep 1 →
`[N₀,N₀+N₁)`, …

`RecoverTrackId(pIdx)` inverts this:

1. **photon → genstep**: walk `sev->gs` (per-genstep summaries with
   `offset`/`photons`) to find the genstep whose range contains `pIdx`;
2. **genstep → G4 track**: look up `fGenstepTrackIds`, a provenance map the
   plugin's stepping action filled the moment the genstep was collected
   (while the emitting track was alive); fallback: the track id the
   scintillation process stamped into the genstep `quad6` payload itself.

Input-photon mode (pure-photon tests) has no charged parent → `-1`.

Downstream, `SummarizePhotonHits` passes this raw track id through the
persistency manager's `fTrackIdMap`, mapping it to the nearest *stored
trajectory* — so the final `TG4PhotonHit.PrimaryId` points into
`TG4Event.Trajectories` identically for CPU and GPU photons. The plugin only
supplies a valid raw G4 track id; edep-sim does the remap.

### 3.2 Creator process: matching Geant4's WLS convention

CPU hits report the creator process of the *arriving* photon. Geant4 makes a
NEW track at WLS re-emission, so a TPB-shifted photon arrives with creator
`OpWLS`, not `Scintillation`. Opticks instead keeps one continuous photon,
but its `flagmask` remembers re-emission. The decode therefore checks
re-emission first:

```
flagmask & BULK_REEMIT    → (fOptical=3,         fOpWLS=34)
flagmask & CERENKOV       → (fElectromagnetic=2, fCerenkov=21)
flagmask & SCINTILLATION  → (fElectromagnetic=2, fScintillation=22)
otherwise (input/torch)   → (fUserDefined=6,     0)
```

Verified in a DUAL run: both CPU (575 hits) and GPU (408 hits) detected
populations report `TG4PhotonHit.Process == 34`.

### 3.3 Insertion mechanics

```cpp
G4Event* ev = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
G4HCofThisEvent* hce = ev->GetHCofThisEvent();      // created if absent
hce->AddHitsCollection(hcid, coll);
```

- `hcid` is resolved once per run via
  `GetCollectionID("SimphonyPhotonDetector/SimphonyHits")` — the fixed slot
  index that `SummarizePhotonDetectors` later reads with `GetHC(HCId)`.
- `EndOfEventAction` gets a `const G4Event*`; inserting a collection mutates
  the event, hence `GetNonconstCurrentEvent()`. Legal: the event is still
  live; the persistency manager reads it only afterwards.
- Ownership of the collection and hits transfers to the event.

## 4. The empty-collection invariant (crash avoidance)

`FillPhotonDetectorHits` is called on **every** path out of
`EndOfEventAction` — including the "no GPU work" early returns — and inserts
the collection **even when it is empty** (the fill loop just doesn't run).

Why: Geant4's volume-attached SDs re-create and insert their collection every
event from their `Initialize(HCE)` callback, so edep-sim could historically
assume "registered slot ⇒ non-null collection". Our pseudo-SD is detached, so
`Initialize()` is never called and nothing inserts a collection
automatically. An event without the insert leaves `GetHC(hcid)` returning
**null**, and the summarizer's `g4Hits->GetSize()` dereferences it →
segfault (observed on the first `photon1` validation run, exit 139, on an
event with zero GPU hits). An empty collection instead yields
`GetSize()==0`, which the summarizer skips gracefully.

## 5. Required edep-sim updates and why

Both are in edep-sim proper (rebuild with `make -j4 install`), are minimal,
and are safe for all existing users — upstream-PR material.

### 5.1 `HitSurface` value constructor (`src/EDepSimHitSurface.hh/.cc`)

```cpp
HitSurface(int primaryId, double energyDeposit,
           const G4LorentzVector& position, const G4LorentzVector& start,
           int pdgEncoding, int creatorType, int creatorSubtype);
```

`HitSurface` had **private fields, getters only, and two constructors**:
default and `HitSurface(const G4Step*)`. The step constructor reaches
*through* the step into a live `G4Track`, its particle definition, and its
creator-process pointer — objects that do not exist for a GPU photon
(transport happened outside Geant4; by end-of-event the per-track G4Step has
been overwritten anyway, and the honest creator-process instance may not even
be constructed in GPU-only mode). Fabricating a synthetic
`G4DynamicParticle`+`G4Track`+`G4Step` per hit just so the constructor can
copy seven values back out is possible but fragile and was rejected.

The value constructor is a pure initializer list: no logic, no behaviour
change for existing callers, no change to the persistent `TG4PhotonHit` or
ROOT I/O. It matches intent already documented in `TG4PhotonHit` ("PrimaryId
/ Start … may not be available if photon tracking is offloaded"). Without it
the plugin cannot construct the one object type the persistency manager
accepts for `PhotonDetectors` (subclassing is impossible — the fields are
private, and the only friend is `EDepSim::PersistencyManager`).

### 5.2 Null-guards in the summarizers (`src/EDepSimPersistencyManager.cc`)

```cpp
G4VHitsCollection* g4Hits = HCofEvent->GetHC(HCId);
if (!g4Hits || g4Hits->GetSize()<1) continue;    // was: g4Hits->GetSize()<1
```

Applied in **both** `SummarizePhotonDetectors` and
`SummarizeSegmentDetectors` (same pattern). The `G4HCtable` lists every
collection that *could* exist; `GetHC()` returns what was *actually inserted
this event*, which can be null for any registered-but-unfilled slot. The old
code dereferenced unconditionally — a latent crash for any externally filled
collection, first triggered by our detached pseudo-SD (§4). The guard makes a
missing collection equivalent to an empty one.

§4 and §5.2 are belt-and-suspenders for the same failure: either alone stops
the crash; together, neither side depends on the other's discipline.

## 6. Runtime switches

| Env var | Default | Effect |
|---|---|---|
| `EDEP_SIMPHONY_LEGACY_HITTREE` | `1` | Also fill the flat `GPUPhotonHits` TTree (kept while analysis scripts migrate; set `0` to disable). |
| `EDEP_SIMPHONY_DEBUG_FILE` | `simphony_debug.root` | Where the plugin TTrees go when the ROOT persistency manager is not available (tier-2 fallback in `SimphonyRunAction::AcquireOutputFile()`). |
| `EDEP_SIMPHONY_FORCE_NO_ROOTPM` | unset | `1` forces the tier-2 fallback even in the CLI — test hook that emulates the Phlex situation. |

PhotonDetectors filling itself has **no switch**: it is always on and needs
no ROOT file at all.

## 7. Reading the result

```python
import ROOT
f = ROOT.TFile.Open("output/gpu.root")
t = f.Get("EDepSimEvents")
t.GetEntry(0); ev = t.Event
for sd, hits in ev.PhotonDetectors:          # sd is the SD name string
    for h in hits:
        h.GetStop()        # detection position+time (TLorentzVector, mm/ns)
        h.GetStart()       # creation point (zero unless DebugHeavy)
        h.GetWavelength()  # nm-equivalent (derived from EnergyDeposit)
        h.GetPrimaryId()   # stored-trajectory id of the charged parent
        h.GetProcess()     # creator subtype: 21 Cerenkov / 22 Scint / 34 WLS
```

`tests_benchmark/check_photondet_migration.py` validates
`PhotonDetectors["SimphonyPhotonDetector"]` ≡ legacy `GPUPhotonHits` 1:1
(count, position, time, wavelength) per event.
