#ifndef SimphonyPhotonSD_hh
#define SimphonyPhotonSD_hh

#include <G4VSensitiveDetector.hh>

/// A pseudo sensitive detector for the GPU (Simphony) optical photon hits.
///
/// This SD is never attached to a logical volume and its ProcessHits() is
/// never called: GPU photons are transported outside Geant4, so there is no
/// G4Step to process.  It exists so that the G4SDManager hit-collection table
/// has a registered slot named "SimphonyPhotonDetector/SimphonyHits".  The
/// plugin fills a HitSurfaceCollection under that slot at EndOfEventAction
/// (see SimphonyEventAction::FillPhotonDetectorHits), and the edep-sim
/// PersistencyManager::SummarizePhotonDetectors then copies the hits into
/// TG4Event.PhotonDetectors["SimphonyPhotonDetector"] — for ANY persistency
/// backend (edep-sim CLI ROOT file, Phlex, ...).
///
/// The SD name is deliberately different from the GDML sensor SD name so
/// that CPU-tracked (Geant4) photons and GPU (Simphony) photons stay
/// distinguishable inside the same TG4Event.
class SimphonyPhotonSD : public G4VSensitiveDetector {
public:
    static constexpr const char* kSDName = "SimphonyPhotonDetector";
    static constexpr const char* kHCName = "SimphonyHits";

    SimphonyPhotonSD() : G4VSensitiveDetector(kSDName) {
        collectionName.insert(kHCName);
    }

    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override { return false; }
};

#endif
