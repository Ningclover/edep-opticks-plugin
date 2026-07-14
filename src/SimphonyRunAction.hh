#pragma once
#include <G4UserRunAction.hh>
#include <cstdint>

class TTree;
class TBranch;
class TFile;
struct GPUHitRecord;

/// Initialises G4CXOpticks at begin-of-run, opens a parallel GPU hit TTree
/// in the same ROOT file as edep-sim, and tears everything down at end-of-run.
class SimphonyRunAction : public G4UserRunAction
{
public:
    explicit SimphonyRunAction(const char* option);
    ~SimphonyRunAction() override = default;

    void BeginOfRunAction(const G4Run* run) override;
    void EndOfRunAction(const G4Run* run) override;

    /// Accessed by SimphonyEventAction to fill the trees each event.
    TTree* GetGPUTree()      const { return fGPUTree; }
    TTree* GetGPUTrackTree() const { return fGPUTrackTree; }
    TTree* GetGPUStepTree()  const { return fGPUStepTree; }

    /// CPU photon fate tree (filled by SimphonyCpuPhotonTracker). Created here
    /// so it lands in the same ROOT file.
    void   CreateCPUTrackTree(TFile* f);
    TTree* GetCPUTrackTree()  const { return fCPUTrackTree; }
    TTree* GetCPUStepTree()   const { return fCPUStepTree; }
    /// Singleton lookup so the dlopen'd tracking-action factory can find the
    /// active run action (and thus the CPU tree) without extra plumbing.
    static SimphonyRunAction* Instance() { return sInstance; }

    // ── Debug accumulators (process-wide, reset at BeginOfRun) ───────────
    // (4) total photon count from edep-sim's DokeBirks visE / 19.5 eV
    // (5) total photons accepted by eic-opticks (sum of sgs.photons)
    static double   sTotalEdepSimPhotons;     // photon-equivalents from DokeBirks
    static int64_t  sTotalSimphonyGenstepPhotons;
    static int64_t  sStepCountInLAr;          // for sampling stats
    static int64_t  sNewGenstepsThisRun;

private:
    /// Get the ROOT file the plugin trees should live in.  Two-tier fallback:
    /// (1) the edep-sim RootPersistencyManager's file when that concrete
    ///     manager is in use (edep-sim CLI) — trees land next to EDepSimEvents;
    /// (2) otherwise (e.g. Phlex, where the persistency manager is a thin
    ///     TG4Event-only class) a plugin-owned file, path from
    ///     EDEP_SIMPHONY_DEBUG_FILE (default ./simphony_debug.root), closed in
    ///     EndOfRunAction.
    /// Returns nullptr if neither is available.
    TFile* AcquireOutputFile();

    TTree* fGPUTree      = nullptr;
    TTree* fGPUTrackTree = nullptr;  // all photons + fate
    TTree* fGPUStepTree  = nullptr;  // full trajectory points
    TTree* fCPUTrackTree = nullptr;  // CPU photon fate
    TTree* fCPUStepTree  = nullptr;  // CPU full trajectory points

    TFile* fOwnedFile    = nullptr;  // non-null only on fallback tier (2)

    static SimphonyRunAction* sInstance;
};
