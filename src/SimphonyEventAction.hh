#pragma once
#include <G4UserEventAction.hh>
#include <vector>
#include <cstdint>

class SimphonyRunAction;

/// Runs the GPU optical photon simulation at end of each event and fills
/// the GPUPhotonHits TTree with the resulting hits + provenance.
class SimphonyEventAction : public G4UserEventAction
{
public:
    explicit SimphonyEventAction(SimphonyRunAction* runAction);
    ~SimphonyEventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    /// Called by SimphonyStepAction to record genstep→TrackID provenance.
    void RecordGenstep(int64_t genstepIdx, int trackId);

    /// Register the step action so BeginOfEventAction can reset its baseline.
    void SetStepAction(class SimphonyStepAction* sa) { fStepAction = sa; }

private:
    SimphonyRunAction*   fRunAction;
    class SimphonyStepAction* fStepAction = nullptr;

    std::vector<std::pair<int64_t, int>> fGenstepTrackIds;
    int64_t fGenstepCountAtEventStart = 0;

    int RecoverTrackId(int photonIndex) const;

    /// Save ALL photons (detected or not) + final fate to GPUPhotonTracks and
    /// the full per-photon trajectory points to GPUPhotonSteps.
    void SaveAllPhotonsAndTrajectories(const G4Event* event, int64_t nphoton);

    /// Convert the GPU hits to EDepSim::HitSurface objects and insert them as
    /// the "SimphonyPhotonDetector/SimphonyHits" collection of this event, so
    /// the edep-sim persistency manager (ANY backend) copies them into
    /// TG4Event.PhotonDetectors["SimphonyPhotonDetector"].
    void FillPhotonDetectorHits(const G4Event* event, int64_t nhit);

    /// Pure-photon mode (EDEP_SIMPHONY_INPUT_PHOTONS=1): capture the event's
    /// primary optical photons and inject them into Opticks as input photons.
    void InjectPrimaryPhotons(const G4Event* event);
};
