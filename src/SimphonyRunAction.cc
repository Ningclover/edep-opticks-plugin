#include "SimphonyRunAction.hh"
#include "SimphonyEventAction.hh"

#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4TransportationManager.hh>
#include <G4VPersistencyManager.hh>
#include <G4UserStackingAction.hh>
#include <G4ParticleGun.hh>
#include <G4VPrimaryGenerator.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4SystemOfUnits.hh>
#include "EDepSimUserStackingAction.hh"

#include "G4CXOpticks.hh"
#include "SEvt.hh"
#include "SEventConfig.hh"
#include "OpticksPhoton.hh"
#include "LArTPCSensorIdentifier.h"
#include "SimphonyPhotonSD.hh"
#include <G4SDManager.hh>
// OPTICKS_INTEGRATION_MODE=1 must be set in the environment for GPU-only mode.
// G4CXOpticks::SetGeometry() calls SEvt::CreateOrReuse() internally.

// edep-sim ROOT persistency
#include "EDepSimRootPersistencyManager.hh"

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <iostream>
#include <cstdlib>
#include <string>

// Per-event branch buffers (filled by SimphonyEventAction, stored here as statics
// so TTree branches keep valid pointers across events).
static int              gEventId   = -1;
static int              gTrackId   = -1;
static int              gProcess   = -1;  // 0=Cerenkov, 2=Scintillation
static float            gWavelength = 0;  // nm
static TLorentzVector   gHitPos;          // mm, ns
static TLorentzVector   gStartPos;        // mm, ns (genstep step position)

// ── ALL-photon tree buffers (GPUPhotonTracks): one row per photon, detected
//    or not, with the final fate (terminating flag) so we know WHY it was not
//    detected.
static int            gtEventId   = -1;
static int            gtPhotonId  = -1;  // photon index within the event
static int            gtTrackId   = -1;  // GPU photon track id (== PhotonId; one
                                         // per photon, persists through WLS)
static int            gtParentId  = -1;  // charged parent G4 track that made the
                                         // genstep, or -1 for input photons
static int            gtProcess   = -1;  // 0=Cerenkov, 2=Scintillation
static float          gtWavelength = 0;  // nm
static unsigned       gtFlag      = 0;   // final/terminating flag bit
static unsigned       gtFlagMask  = 0;   // OR of all history flags
static int            gtDetected  = 0;   // 1 if SURFACE_DETECT/EFFICIENCY_COLLECT
static char           gtFate[8]   = {0}; // short reason: "SD","AB","SA","MI",...
static int            gtNStep     = 0;   // number of recorded trajectory points
static TLorentzVector gtStartPos;        // mm, ns (first record point)
static TLorentzVector gtEndPos;          // mm, ns (last record point / hit)

// ── Per-step trajectory tree buffers (GPUPhotonSteps): one row per step point
//    of every photon, so the full bounce-by-bounce path is saved.
static int            gsEventId  = -1;
static int            gsPhotonId = -1;
static int            gsStep     = -1;   // step index along this photon
static unsigned       gsFlag     = 0;    // flag at this step point
static char           gsFlagAbbr[8] = {0};
static TLorentzVector gsPos;             // mm, ns

// ── CPU photon fate tree buffers (CPUPhotonTracks): one row per CPU-tracked
//    optical photon, detected or not, with the reason it ended.
static int            gcEventId  = -1;
static int            gcTrackId  = -1;
static int            gcParentId = -1;
static float          gcWavelength = 0;  // nm
static float          gcEnergy   = 0;    // eV
static int            gcDetected = 0;    // 1 if Detection
static char           gcFate[16] = {0};  // "SD","AB","SA","WLS","MI","Reflect"
static char           gcEndVolume[64] = {0};
static int            gcNStep    = 0;     // # geometry steps the photon took
static TLorentzVector gcStartPos;        // mm, ns
static TLorentzVector gcEndPos;          // mm, ns

// ── CPU per-step trajectory tree buffers (CPUPhotonSteps): one row per step
//    point of every CPU-tracked optical photon (full bounce-by-bounce path,
//    which edep-sim's Trajectories branch collapses to start+end).
static int            gpsEventId  = -1;
static int            gpsTrackId  = -1;
static int            gpsStep     = -1;   // G4 step number along this photon
static char           gpsProc[24] = {0}; // process that ended the step
static char           gpsVolume[48] = {0};
static TLorentzVector gpsPos;            // post-step position (mm, ns)

// Static debug accumulators (declared in SimphonyRunAction.hh).
double   SimphonyRunAction::sTotalEdepSimPhotons        = 0.0;
int64_t  SimphonyRunAction::sTotalSimphonyGenstepPhotons = 0;
int64_t  SimphonyRunAction::sStepCountInLAr             = 0;
int64_t  SimphonyRunAction::sNewGenstepsThisRun         = 0;
SimphonyRunAction* SimphonyRunAction::sInstance         = nullptr;

SimphonyRunAction::SimphonyRunAction(const char* /*option*/)
{
    sInstance = this;
}

void SimphonyRunAction::CreateCPUTrackTree(TFile* f)
{
    f->cd();
    fCPUTrackTree = new TTree("CPUPhotonTracks",
                              "CPU (edep-sim) optical photons (all) + final fate");
    fCPUTrackTree->Branch("EventId",    &gcEventId);
    fCPUTrackTree->Branch("TrackId",    &gcTrackId);
    fCPUTrackTree->Branch("ParentId",   &gcParentId);
    fCPUTrackTree->Branch("Wavelength", &gcWavelength);  // nm
    fCPUTrackTree->Branch("Energy",     &gcEnergy);       // eV
    fCPUTrackTree->Branch("Detected",   &gcDetected);     // 1=detected
    fCPUTrackTree->Branch("Fate",       gcFate, "Fate/C");
    fCPUTrackTree->Branch("EndVolume",  gcEndVolume, "EndVolume/C");
    fCPUTrackTree->Branch("NStep",      &gcNStep);
    fCPUTrackTree->Branch("StartPos",   "TLorentzVector", &gcStartPos);
    fCPUTrackTree->Branch("EndPos",     "TLorentzVector", &gcEndPos);

    // Full CPU per-photon trajectory (every step), mirroring GPUPhotonSteps.
    fCPUStepTree = new TTree("CPUPhotonSteps",
                             "CPU optical photon trajectory points (all steps)");
    fCPUStepTree->Branch("EventId",  &gpsEventId);
    fCPUStepTree->Branch("TrackId",  &gpsTrackId);
    fCPUStepTree->Branch("Step",     &gpsStep);
    fCPUStepTree->Branch("Process",  gpsProc, "Process/C");
    fCPUStepTree->Branch("Volume",   gpsVolume, "Volume/C");
    fCPUStepTree->Branch("Pos",      "TLorentzVector", &gpsPos);
}

void SimphonyRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    // Reset run-level debug accumulators
    sTotalEdepSimPhotons        = 0.0;
    sTotalSimphonyGenstepPhotons = 0;
    sStepCountInLAr             = 0;
    sNewGenstepsThisRun         = 0;

    // ── 0. Report the primary generator particle/energy (request item 1) ──
    auto* pgAction = G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
    if (pgAction) {
        // Try to introspect a G4ParticleGun if that's what's in use.
        // GPS (/gps/...) doesn't expose a public "current particle" without
        // running an event, so we fall back to a generic note.
        G4VPrimaryGenerator* gen = nullptr;
        // No portable accessor across edep-sim / Geant4 versions; print what we can
        std::cout << "[SimphonyPlugin][DBG] BeginOfRun: primary-generator action installed; "
                  << "actual particle/energy printed on first event.\n";
        (void)gen;
    } else {
        std::cout << "[SimphonyPlugin][DBG] BeginOfRun: no primary-generator action found\n";
    }

    // ── 1. Translate Geant4 geometry → CSGFoundry and init OptiX ──────────
    G4Navigator* nav = G4TransportationManager::GetTransportationManager()
                           ->GetNavigatorForTracking();
    const G4VPhysicalVolume* world = nav->GetWorldVolume();
    if (!world) {
        std::cerr << "[SimphonyPlugin] ERROR: world volume not available in BeginOfRunAction\n";
        return;
    }

    // Must call SEvt::CreateOrReuse() BEFORE SetGeometry so that scontext
    // initializes CUDA and sets SEventConfig::_DeviceName (HasDevice()=true).
    // Without this, SetGeometry skips CSGOptiX::Create and SEvt stays null.
    // This mirrors what G4CXOpticks::SetGeometry_JUNO does.
    // OPTICKS_INTEGRATION_MODE=1 (GPU-only) creates the EGPU SEvt instance.
    //
    // To save ALL photons + the full per-photon step record (trajectory), put
    // Opticks in "DebugHeavy" event mode BEFORE SEvt is instantiated (which runs
    // SEventConfig::Initialize, settable only once). DebugHeavy makes the gather
    // mask include photon + record + hit and auto-sets MaxRecord to the record
    // limit. The default "Minimal" mode gathers only hits, so getPhoton() and
    // the record buffer would be empty. Only do this when the all-photon /
    // input-photon debug features are wanted (DUAL or INPUT_PHOTONS), to leave
    // the fast hit-only production path untouched.
    {
        // DebugHeavy (all photons + full trajectory) is only enabled for the
        // input-photon test, which uses few photons. The scintillation DUAL
        // mode makes ~millions of photons/event, where per-photon trajectory
        // recording is impractical (OOM) — it keeps the hit-only path.
        const char* inp_c0  = std::getenv("EDEP_SIMPHONY_INPUT_PHOTONS");
        // EDEP_SIMPHONY_DEBUGHEAVY=1 forces the all-photon + full-trajectory GPU
        // recording on for the scintillation DUAL path too (NOT just input
        // photons). Only safe for FEW-photon runs (e.g. a single low-energy
        // primary) — a multi-event high-stats DUAL run would OOM. Size the
        // device buffers with EDEP_SIMPHONY_MAXSLOT accordingly.
        const char* dbgh_c  = std::getenv("EDEP_SIMPHONY_DEBUGHEAVY");
        auto on = [](const char* c){ return c && (std::string(c)=="1" ||
                          std::string(c)=="true" || std::string(c)=="on"); };
        if (on(inp_c0) || on(dbgh_c)) {
            SEventConfig::SetDebugHeavy();   // gather photon,record,hit + MaxRecord
            // DebugHeavy allocates max_slot * max_record (and *max_prd/*max_aux)
            // device buffers. With the default VRAM-sized max_slot (~millions)
            // that is tens of GB -> CUDA OOM. These debug runs use very few
            // photons, so cap the slot/photon budget. Override with
            // EDEP_SIMPHONY_MAXSLOT if a run needs more.
            const char* ms_c = std::getenv("EDEP_SIMPHONY_MAXSLOT");
            int maxSlot = ms_c ? std::atoi(ms_c) : 200000;  // plenty for the
            if (maxSlot < 1000) maxSlot = 1000;             // pure-photon test
            SEventConfig::SetMaxSlot(maxSlot);
            SEventConfig::SetMaxPhoton(maxSlot);
            const char* mr_c = std::getenv("EDEP_SIMPHONY_MAXRECORD");
            if (mr_c) {
                int mr = std::atoi(mr_c);
                if (mr >= 2) {
                    SEventConfig::SetMaxRecord(mr);
                    if (SEventConfig::MaxBounce() > mr - 1)
                        SEventConfig::SetMaxBounce(mr - 1);
                }
            }
            // Decouple the TRANSPORT bounce limit from the RECORD length. The
            // per-photon trajectory buffer is hard-capped at sseq::SLOTS (32),
            // but the photon can be allowed to keep propagating well past 32
            // bounces (only the first ~32 points are stored). Set
            // EDEP_SIMPHONY_MAXBOUNCE to raise the bounce kill independently of
            // MaxRecord, e.g. =1000 to effectively remove the cap so we can see
            // whether the 31-bounce guillotine (not physics) is what kills
            // photons mid-flight. Applied AFTER the MaxRecord clamp above so it
            // wins.
            const char* mb_c = std::getenv("EDEP_SIMPHONY_MAXBOUNCE");
            if (mb_c) {
                int mb = std::atoi(mb_c);
                if (mb >= 1) SEventConfig::SetMaxBounce(mb);
            }
            std::cout << "[SimphonyPlugin] OPTICKS EventMode=DebugHeavy"
                      << " MaxRecord=" << SEventConfig::MaxRecord()
                      << " MaxBounce=" << SEventConfig::MaxBounce()
                      << " (all photons + full trajectory recording ON)\n";
        }

        // EDEP_SIMPHONY_MAXBOUNCE applies in ALL modes (not just the
        // input-photon debug path above): the bounce kill is a transport
        // parameter independent of whether per-photon trajectories are
        // recorded. In the hit-only DUAL scintillation path it controls how
        // many bounces each GPU photon may take before being killed, e.g.
        // =1000 to remove the default 31-bounce guillotine.
        const char* mb_all = std::getenv("EDEP_SIMPHONY_MAXBOUNCE");
        if (mb_all) {
            int mb = std::atoi(mb_all);
            if (mb >= 1) {
                SEventConfig::SetMaxBounce(mb);
                std::cout << "[SimphonyPlugin] MaxBounce set to "
                          << SEventConfig::MaxBounce()
                          << " (EDEP_SIMPHONY_MAXBOUNCE)\n";
            }
        }
    }
    SEvt::CreateOrReuse();

    // The lighttrap.gdml SurfaceDetector registration calls
    // SetKillOpticalPhotons(false) so Geant4 would track photons on CPU.
    // Override that: we kill photons on CPU and let the GPU handle transport.
    // EDEP_SIMPHONY_DUAL=1 keeps CPU optical-photon tracking ON so stock
    // G4Scintillation's photons reach the SurfaceSD (→ PhotonDetectors),
    // alongside the GPU genstep path. Default (unset/0): kill on CPU.
    const char* dual_c = std::getenv("EDEP_SIMPHONY_DUAL");
    const std::string dual = dual_c ? std::string(dual_c) : std::string("0");
    const bool useDual = (dual == "1" || dual == "true" || dual == "on");

    auto* sa = dynamic_cast<EDepSim::UserStackingAction*>(
                   const_cast<G4UserStackingAction*>(
                       G4RunManager::GetRunManager()->GetUserStackingAction()));
    if (sa) {
        sa->SetKillOpticalPhotons(!useDual);
        std::cout << "[SimphonyPlugin] KillOpticalPhotons=" << (!useDual)
                  << (useDual ? " (DUAL: CPU also tracks photons)"
                              : " (GPU handles transport)") << "\n";
    }

    // Register custom sensor identifier before geometry translation so that
    // eic-opticks recognises our SiPM volumes (which have SDs but aren't
    // named "PMT" as the default identifier requires).
    G4CXOpticks::SetSensorIdentifier(new LArTPCSensorIdentifier());

    G4CXOpticks::SetGeometry(world);
    std::cout << "[SimphonyPlugin] G4CXOpticks geometry set, OptiX ready\n";

    // ── 1b. Register the pseudo-SD for the GPU photon hits ────────────────
    // Gives the G4 hit-collection table a "SimphonyPhotonDetector/SimphonyHits"
    // slot; SimphonyEventAction fills it at end of event and the edep-sim
    // persistency manager folds it into TG4Event.PhotonDetectors (works with
    // ANY persistency backend, no RootPersistencyManager needed).
    if (!G4SDManager::GetSDMpointer()->FindSensitiveDetector(
            SimphonyPhotonSD::kSDName, /*warning=*/false)) {
        G4SDManager::GetSDMpointer()->AddNewDetector(new SimphonyPhotonSD());
        std::cout << "[SimphonyPlugin] Registered pseudo-SD '"
                  << SimphonyPhotonSD::kSDName
                  << "' (GPU hits -> TG4Event.PhotonDetectors)\n";
    }

    // ── 2. Create the plugin TTrees (debug + legacy hits) ─────────────────
    TFile* f = AcquireOutputFile();
    if (!f) {
        std::cout << "[SimphonyPlugin] No ROOT file available; plugin TTrees "
                     "disabled (GPU hits still go to TG4Event.PhotonDetectors)\n";
        return;
    }
    f->cd();

    // Legacy flat GPU hit tree. The same hits now also land in
    // TG4Event.PhotonDetectors["SimphonyPhotonDetector"]; keep the tree while
    // analysis scripts migrate. Disable with EDEP_SIMPHONY_LEGACY_HITTREE=0.
    const char* legacy_c = std::getenv("EDEP_SIMPHONY_LEGACY_HITTREE");
    const bool legacyTree = !legacy_c || !(std::string(legacy_c) == "0"
                             || std::string(legacy_c) == "false"
                             || std::string(legacy_c) == "off");
    if (legacyTree) {
        fGPUTree = new TTree("GPUPhotonHits",
                             "GPU optical photon hits (eic-opticks)");
        fGPUTree->Branch("EventId",    &gEventId);
        fGPUTree->Branch("TrackId",    &gTrackId);    // G4 TrackID of charged parent
        fGPUTree->Branch("Process",    &gProcess);    // 0=Cerenkov, 2=Scintillation
        fGPUTree->Branch("Wavelength", &gWavelength); // nm
        fGPUTree->Branch("HitPos",    "TLorentzVector", &gHitPos);    // mm, ns
        fGPUTree->Branch("StartPos",  "TLorentzVector", &gStartPos);  // mm, ns
        std::cout << "[SimphonyPlugin] GPUPhotonHits TTree created in ROOT file\n";
    }

    // ── ALL photons (detected or not) + fate ──────────────────────────────
    fGPUTrackTree = new TTree("GPUPhotonTracks",
                              "GPU optical photons (all) + final fate");
    fGPUTrackTree->Branch("EventId",    &gtEventId);
    fGPUTrackTree->Branch("PhotonId",   &gtPhotonId);
    fGPUTrackTree->Branch("TrackId",    &gtTrackId);     // one per photon (no WLS split)
    fGPUTrackTree->Branch("ParentId",   &gtParentId);    // charged parent, -1 if input
    fGPUTrackTree->Branch("Process",    &gtProcess);     // 0=Cerenkov 2=Scint
    fGPUTrackTree->Branch("Wavelength", &gtWavelength);  // nm
    fGPUTrackTree->Branch("Flag",       &gtFlag);        // terminating flag bit
    fGPUTrackTree->Branch("FlagMask",   &gtFlagMask);    // OR of history flags
    fGPUTrackTree->Branch("Detected",   &gtDetected);    // 1=detected
    fGPUTrackTree->Branch("Fate",       gtFate, "Fate/C"); // "SD","AB","SA","MI"
    fGPUTrackTree->Branch("NStep",      &gtNStep);       // # trajectory points
    fGPUTrackTree->Branch("StartPos",   "TLorentzVector", &gtStartPos);
    fGPUTrackTree->Branch("EndPos",     "TLorentzVector", &gtEndPos);

    // ── full per-photon trajectory points ─────────────────────────────────
    fGPUStepTree = new TTree("GPUPhotonSteps",
                             "GPU optical photon trajectory points (all bounces)");
    fGPUStepTree->Branch("EventId",  &gsEventId);
    fGPUStepTree->Branch("PhotonId", &gsPhotonId);
    fGPUStepTree->Branch("Step",     &gsStep);
    fGPUStepTree->Branch("Flag",     &gsFlag);
    fGPUStepTree->Branch("FlagAbbr", gsFlagAbbr, "FlagAbbr/C");
    fGPUStepTree->Branch("Pos",      "TLorentzVector", &gsPos);

    std::cout << "[SimphonyPlugin] GPUPhotonTracks + GPUPhotonSteps TTrees created\n";

    // ── CPU photon fate tree (filled by the external tracking action) ──────
    // The full CPU trajectory path is saved by edep-sim itself into the
    // Trajectories branch when /edep/db/set/savePhotonTrajectories true is set;
    // here we add the per-photon FATE (why it was / wasn't detected).
    CreateCPUTrackTree(f);
    std::cout << "[SimphonyPlugin] CPUPhotonTracks TTree created\n";
}

TFile* SimphonyRunAction::AcquireOutputFile()
{
    // Tier 1: the edep-sim CLI ROOT persistency manager's own file, so the
    // plugin trees sit next to EDepSimEvents (historical behaviour). The
    // downcast fails for other persistency managers (e.g. Phlex's thin
    // TG4Event-only one) — that is expected, not an error.
    // EDEP_SIMPHONY_FORCE_NO_ROOTPM=1 skips this tier (for testing tier 2).
    const char* force_c = std::getenv("EDEP_SIMPHONY_FORCE_NO_ROOTPM");
    const bool forceNoPM = force_c && (std::string(force_c) == "1"
                             || std::string(force_c) == "true"
                             || std::string(force_c) == "on");
    if (!forceNoPM) {
        auto* pm = dynamic_cast<EDepSim::RootPersistencyManager*>(
                       G4VPersistencyManager::GetPersistencyManager());
        if (pm && pm->IsOpen()) return pm->GetTFile();
        std::cout << "[SimphonyPlugin] edep-sim RootPersistencyManager "
                  << (pm ? "has no open file" : "not in use")
                  << "; falling back to a plugin-owned debug file\n";
    }

    // Tier 2: a plugin-owned file. Written + closed in EndOfRunAction.
    const char* path_c = std::getenv("EDEP_SIMPHONY_DEBUG_FILE");
    const std::string path = path_c ? path_c : "simphony_debug.root";
    fOwnedFile = TFile::Open(path.c_str(), "RECREATE");
    if (!fOwnedFile || fOwnedFile->IsZombie()) {
        std::cerr << "[SimphonyPlugin] WARNING: cannot open '" << path
                  << "' for the plugin TTrees\n";
        delete fOwnedFile;
        fOwnedFile = nullptr;
        return nullptr;
    }
    std::cout << "[SimphonyPlugin] plugin TTrees -> " << path << "\n";
    return fOwnedFile;
}

void SimphonyRunAction::EndOfRunAction(const G4Run* /*run*/)
{
    // ── Run-level debug summary (request items 4 & 5) ─────────────────────
    std::cout << "\n[SimphonyPlugin][DBG] ===== Run summary =====\n"
              << "[SimphonyPlugin][DBG]   Charged steps in LAr (sampled population)  : "
              << sStepCountInLAr << "\n"
              << "[SimphonyPlugin][DBG]   New gensteps emitted to eic-opticks        : "
              << sNewGenstepsThisRun << "\n"
              << "[SimphonyPlugin][DBG] (4) Total photons generated by edep-sim      : "
              << sTotalEdepSimPhotons
              << "  [DokeBirks visE / 19.5 eV summed over charged steps]\n"
              << "[SimphonyPlugin][DBG] (5) Total photons accepted by eic-opticks    : "
              << sTotalSimphonyGenstepPhotons
              << "  [sum of sgs.photons across all collected gensteps]\n"
              << "[SimphonyPlugin][DBG] =========================\n";

    // Plugin-owned fallback file (tier 2): the edep-sim persistency manager
    // knows nothing about it, so write + close it here. (In tier 1 the
    // RootPersistencyManager writes the file on /edep/db/close.)
    if (fOwnedFile) {
        fOwnedFile->Write();
        fOwnedFile->Close();
        std::cout << "[SimphonyPlugin] plugin debug file written: "
                  << fOwnedFile->GetName() << "\n";
        delete fOwnedFile;
        fOwnedFile = nullptr;
        fGPUTree = fGPUTrackTree = fGPUStepTree = nullptr;
        fCPUTrackTree = fCPUStepTree = nullptr;
    }

    G4CXOpticks::Finalize();
    std::cout << "[SimphonyPlugin] G4CXOpticks finalized\n";
}

// ── Accessors used by SimphonyEventAction ────────────────────────────────────
int&            GetGEventId()    { return gEventId;   }
int&            GetGTrackId()    { return gTrackId;   }
int&            GetGProcess()    { return gProcess;   }
float&          GetGWavelength() { return gWavelength; }
TLorentzVector& GetGHitPos()     { return gHitPos;    }
TLorentzVector& GetGStartPos()   { return gStartPos;  }

// GPUPhotonTracks (all photons + fate)
int&            GetGtEventId()    { return gtEventId;   }
int&            GetGtPhotonId()   { return gtPhotonId;  }
int&            GetGtTrackId()    { return gtTrackId;   }
int&            GetGtParentId()   { return gtParentId;  }
int&            GetGtProcess()    { return gtProcess;   }
float&          GetGtWavelength() { return gtWavelength; }
unsigned&       GetGtFlag()       { return gtFlag;      }
unsigned&       GetGtFlagMask()   { return gtFlagMask;  }
int&            GetGtDetected()   { return gtDetected;  }
char*           GetGtFate()       { return gtFate;      }
int&            GetGtNStep()      { return gtNStep;     }
TLorentzVector& GetGtStartPos()   { return gtStartPos;  }
TLorentzVector& GetGtEndPos()     { return gtEndPos;    }

// GPUPhotonSteps (trajectory points)
int&            GetGsEventId()  { return gsEventId;  }
int&            GetGsPhotonId() { return gsPhotonId; }
int&            GetGsStep()     { return gsStep;     }
unsigned&       GetGsFlag()     { return gsFlag;     }
char*           GetGsFlagAbbr() { return gsFlagAbbr; }
TLorentzVector& GetGsPos()      { return gsPos;      }

// CPUPhotonTracks (CPU photon fate)
int&            GetGcEventId()    { return gcEventId;   }
int&            GetGcTrackId()    { return gcTrackId;   }
int&            GetGcParentId()   { return gcParentId;  }
float&          GetGcWavelength() { return gcWavelength; }
float&          GetGcEnergy()     { return gcEnergy;    }
int&            GetGcDetected()   { return gcDetected;  }
char*           GetGcFate()       { return gcFate;      }
char*           GetGcEndVolume()  { return gcEndVolume; }
int&            GetGcNStep()      { return gcNStep;     }
TLorentzVector& GetGcStartPos()   { return gcStartPos;  }
TLorentzVector& GetGcEndPos()     { return gcEndPos;    }

// CPUPhotonSteps (CPU trajectory points)
int&            GetGpsEventId()  { return gpsEventId; }
int&            GetGpsTrackId()  { return gpsTrackId; }
int&            GetGpsStep()     { return gpsStep;    }
char*           GetGpsProc()     { return gpsProc;    }
char*           GetGpsVolume()   { return gpsVolume;  }
TLorentzVector& GetGpsPos()      { return gpsPos;     }
