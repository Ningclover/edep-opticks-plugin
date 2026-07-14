#include "SimphonyEventAction.hh"
#include "SimphonyRunAction.hh"
#include "SimphonyStepAction.hh"

#include <G4Event.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>

#include "G4CXOpticks.hh"
#include "SEvt.hh"
#include "SEventConfig.hh"
#include "OpticksPhoton.hh"
#include "sphoton.h"
#include "sevent.h"
#include "sgs.h"
#include "squad.h"
#include "NP.hh"
#include "NPFold.h"
#include "SComp.h"
#include "SEvent.hh"
#include "OpticksGenstep.h"

#include <G4PhysicalConstants.hh>
#include <cstring>
#include <cstdlib>
#include <vector>

// GPU hits -> TG4Event.PhotonDetectors conduit
#include "SimphonyPhotonSD.hh"
#include "EDepSimHitSurface.hh"
#include <G4EventManager.hh>
#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4LorentzVector.hh>
#include <G4OpticalPhoton.hh>
#include <G4ProcessType.hh>
#include <G4EmProcessSubType.hh>
#include <G4OpProcessSubType.hh>

// EDEP_SIMPHONY_INPUT_PHOTONS=1 -> capture the event's PRIMARY optical photons
// and inject them into Opticks as input photons (instead of relying on
// scintillation gensteps). Used by the pure-128nm-photon comparison so the GPU
// transports the SAME photons the CPU tracks. Cached at first use.
static bool InputPhotonMode()
{
    static int cached = -1;
    if (cached < 0) {
        const char* c = std::getenv("EDEP_SIMPHONY_INPUT_PHOTONS");
        cached = (c && (std::string(c) == "1" || std::string(c) == "true"
                        || std::string(c) == "on")) ? 1 : 0;
    }
    return cached == 1;
}

#include <TTree.h>
#include <TLorentzVector.h>

#include <iostream>

// Branch buffer accessors defined in SimphonyRunAction.cc
int&            GetGEventId();
int&            GetGTrackId();
int&            GetGProcess();
float&          GetGWavelength();
TLorentzVector& GetGHitPos();
TLorentzVector& GetGStartPos();

int&            GetGtEventId();
int&            GetGtPhotonId();
int&            GetGtTrackId();
int&            GetGtParentId();
int&            GetGtProcess();
float&          GetGtWavelength();
unsigned&       GetGtFlag();
unsigned&       GetGtFlagMask();
int&            GetGtDetected();
char*           GetGtFate();
int&            GetGtNStep();
TLorentzVector& GetGtStartPos();
TLorentzVector& GetGtEndPos();

int&            GetGsEventId();
int&            GetGsPhotonId();
int&            GetGsStep();
unsigned&       GetGsFlag();
char*           GetGsFlagAbbr();
TLorentzVector& GetGsPos();

// Geant4 mm → ROOT mm (same), Geant4 ns → ROOT ns (same)
static constexpr double kMM = 1.0;
static constexpr double kNS = 1.0;
// Photon energy → wavelength: E [MeV] → wavelength [nm]
// hc = 1.23984e-3 MeV·nm
static constexpr double kHC_MeV_nm = 1.23984e-3;

// OpticksGenstep gentype codes
static constexpr int kGenCerenkov      = 1;  // OpticksGenstep_G4Cerenkov_1042
static constexpr int kGenScintillation = 2;  // OpticksGenstep_DsG4Scintillation_r4695

SimphonyEventAction::SimphonyEventAction(SimphonyRunAction* runAction)
    : fRunAction(runAction)
{}

void SimphonyEventAction::BeginOfEventAction(const G4Event* event)
{
    fGenstepTrackIds.clear();

    SEvt* sev = SEvt::Get_EGPU();
    fGenstepCountAtEventStart = sev ? sev->getNumGenstepCollected() : 0;

    // Sync the step action's baseline to match
    if (fStepAction) fStepAction->ResetGenstepBaseline(fGenstepCountAtEventStart);

    // ── (1) Report primary particle for this event ───────────────────────
    int eid = event->GetEventID();
    int nVtx = event->GetNumberOfPrimaryVertex();
    for (int iv = 0; iv < nVtx; ++iv) {
        G4PrimaryVertex* v = event->GetPrimaryVertex(iv);
        if (!v) continue;
        for (G4PrimaryParticle* p = v->GetPrimary(); p != nullptr; p = p->GetNext()) {
            const G4ParticleDefinition* pdef = p->GetParticleDefinition();
            G4String pname = pdef ? pdef->GetParticleName() : G4String("(unknown)");
            double  ke_MeV = p->GetKineticEnergy() / MeV;
            std::cout << "[SimphonyPlugin][DBG] Event " << eid
                      << " primary: " << pname
                      << "  KE=" << ke_MeV << " MeV"
                      << "  dir=(" << p->GetMomentumDirection().x() << ","
                                   << p->GetMomentumDirection().y() << ","
                                   << p->GetMomentumDirection().z() << ")\n";
        }
    }

    // Pure-photon mode: capture primaries and hand them to Opticks as input
    // photons, then run the SEvt input-photon lifecycle.
    if (InputPhotonMode()) InjectPrimaryPhotons(event);
}

// Build an sphoton input-photon array from the event's primary optical photons
// and register it with the GPU SEvt, then start the input-photon lifecycle.
void SimphonyEventAction::InjectPrimaryPhotons(const G4Event* event)
{
    SEvt* sev = SEvt::Get_EGPU();
    if (!sev) return;

    std::vector<sphoton> ph;
    int nVtx = event->GetNumberOfPrimaryVertex();
    for (int iv = 0; iv < nVtx; ++iv) {
        G4PrimaryVertex* v = event->GetPrimaryVertex(iv);
        if (!v) continue;
        const G4ThreeVector pos = v->GetPosition();   // mm
        const double t = v->GetT0();                  // ns
        for (G4PrimaryParticle* p = v->GetPrimary(); p; p = p->GetNext()) {
            const G4ParticleDefinition* pdef = p->GetParticleDefinition();
            if (!pdef || pdef->GetPDGEncoding() != -22) continue;  // optical only

            G4ThreeVector dir = p->GetMomentumDirection();
            G4ThreeVector pol = p->GetPolarization();
            double E = p->GetKineticEnergy();          // MeV
            double wl_nm = (E > 0) ? (h_Planck * c_light / E) / nm : 0.0;

            sphoton s;
            s.zero();
            s.pos = make_float3(pos.x(), pos.y(), pos.z());
            s.time = static_cast<float>(t);
            s.mom = make_float3(dir.x(), dir.y(), dir.z());
            s.pol = make_float3(pol.x(), pol.y(), pol.z());
            s.wavelength = static_cast<float>(wl_nm);
            ph.push_back(s);
        }
    }

    if (ph.empty()) return;

    // NP array shaped (N,4,4) as required by SEvt::setInputPhoton.
    NP* arr = NP::Make<float>(static_cast<int>(ph.size()), 4, 4);
    std::memcpy(arr->bytes(), ph.data(), ph.size() * sizeof(sphoton));

    // setInputPhoton + a placeholder (identity) frame (the frame set also runs
    // transformInputPhoton so input_photon_transformed is ready).
    sev->setInputPhoton(arr);
    sev->setFramePlaceholder();   // identity model->world transform

    // Explicitly build the INPUT_PHOTON genstep and add it so the GPU sim
    // (QSim) has a genstep to consume — the embedded GPU path does not run the
    // standalone input-photon lifecycle that would otherwise create it.
    NP* igs = SEvent::MakeInputPhotonGenstep(sev->getInputPhoton(),
                                             OpticksGenstep_INPUT_PHOTON,
                                             &sev->fr);
    SEvt::AddGenstep(igs);

    std::cout << "[SimphonyPlugin] InputPhotonMode: injected " << ph.size()
              << " primary optical photons (wl~" << ph[0].wavelength
              << " nm), genstep added\n";
}

void SimphonyEventAction::RecordGenstep(int64_t genstepIdx, int trackId)
{
    fGenstepTrackIds.emplace_back(genstepIdx, trackId);
}

int SimphonyEventAction::RecoverTrackId(int photonIndex) const
{
    // Walk the genstep summary vector to find which genstep produced photon[photonIndex]
    SEvt* sev = SEvt::Get_EGPU();
    if (!sev) return -1;

    const std::vector<sgs>& gs_vec = sev->gs;
    for (const sgs& gs : gs_vec) {
        if (photonIndex >= gs.offset && photonIndex < gs.offset + gs.photons) {
            // Prefer the provenance map recorded by SimphonyStepAction
            for (const auto& [idx, tid] : fGenstepTrackIds) {
                if (idx == gs.index) return tid;
            }
            // Fallback: read trackid directly from the quad6 genstep data
            const quad6* gsdata = sev->getGenstepVecData();
            if (gsdata) return static_cast<int>(gsdata[gs.index].trackid());
            return -1;
        }
    }
    return -1;
}

void SimphonyEventAction::EndOfEventAction(const G4Event* event)
{
    G4CXOpticks* gx = G4CXOpticks::Get();
    if (!gx) {
        std::cout << "[SimphonyPlugin] EndOfEvent: G4CXOpticks::Get() = null\n";
        FillPhotonDetectorHits(event, 0);
        return;
    }

    SEvt* sev = SEvt::Get_EGPU();
    if (!sev) {
        std::cout << "[SimphonyPlugin] EndOfEvent: SEvt::Get_EGPU() = null\n";
        FillPhotonDetectorHits(event, 0);
        return;
    }

    // In input-photon mode the genstep is created inside simulate()'s
    // beginOfEvent(), so getNumGenstepCollected() is 0 here — gate on the
    // input photon instead.
    int64_t ngenstep = sev->getNumGenstepCollected();
    bool haveWork = InputPhotonMode() ? sev->hasInputPhoton() : (ngenstep > 0);
    std::cout << "[SimphonyPlugin] EndOfEvent: ngenstep=" << ngenstep
              << (InputPhotonMode() ? " (input-photon mode)" : "") << "\n";
    if (!haveWork) {
        FillPhotonDetectorHits(event, 0);   // empty collection, see below
        gx->reset(event->GetEventID());
        return;
    }

    // ── 1. Run GPU optical photon simulation ──────────────────────────────
    // simulate() internally runs beginOfEvent (builds the input-photon
    // genstep), transports, gathers each slice and concats the gathered
    // photon/record/hit arrays up to topfold — so getPhoton()/getNumHit()/
    // topfold->get(record) are readable afterwards (no extra gather() needed).
    int64_t nphoton = sev->getNumPhotonCollected();
    std::cout << "[SimphonyPlugin] EndOfEvent: nphoton=" << nphoton << "\n";
    gx->simulate(event->GetEventID(), /*reset=*/false);

    // ── 2. Collect hits ───────────────────────────────────────────────────
    int64_t nhit = SEvt::GetNumHit_EGPU();
    std::cout << "[SimphonyPlugin] EndOfEvent: nhit=" << nhit << "\n";
    TTree* tree = fRunAction->GetGPUTree();

    if (nhit > 0 && tree) {
        int& evId    = GetGEventId();
        int& trackId = GetGTrackId();
        int& proc    = GetGProcess();
        float& wl    = GetGWavelength();
        TLorentzVector& hitPos   = GetGHitPos();
        TLorentzVector& startPos = GetGStartPos();

        evId = event->GetEventID();

        for (int64_t i = 0; i < nhit; ++i) {
            sphoton p;
            sev->getHit(p, static_cast<unsigned>(i));

            // Photon index within the event (used to find its genstep)
            int pIdx = static_cast<int>(p.get_index());

            trackId = RecoverTrackId(pIdx);

            // Process type from flagmask
            unsigned fm = p.flagmask;
            if      (fm & 0x0001) proc = 0;  // Cerenkov
            else if (fm & 0x0002) proc = 2;  // Scintillation
            else                  proc = -1;

            // Wavelength [nm] (stored directly in sphoton)
            wl = p.wavelength; // already in nm

            // Hit position (Geant4 mm, ns)
            hitPos.SetXYZT(p.pos.x * kMM, p.pos.y * kMM, p.pos.z * kMM,
                           p.time * kNS);

            // Start position: not stored per-hit; zero-fill
            startPos.SetXYZT(0, 0, 0, 0);

            tree->Fill();
        }
    }

    // ── 2a. GPU hits → TG4Event.PhotonDetectors["SimphonyPhotonDetector"] ──
    // Persistency-agnostic path (works in edep-sim CLI AND embedded contexts
    // like Phlex): the edep-sim PersistencyManager picks this collection up in
    // SummarizePhotonDetectors when it stores the event (which runs after this
    // user action).
    FillPhotonDetectorHits(event, nhit);

    // ── 2b. Save ALL photons (detected or not) + final fate, and the full
    //        per-photon trajectory (step record buffer). ────────────────────
    SaveAllPhotonsAndTrajectories(event, nphoton);

    // ── 3. Reset SEvt for next event ─────────────────────────────────────
    gx->reset(event->GetEventID());
    fGenstepTrackIds.clear();
}

void SimphonyEventAction::FillPhotonDetectorHits(const G4Event* event,
                                                 int64_t nhit)
{
    static int hcid = -2;   // resolve once; -1 = collection not registered
    if (hcid == -2) {
        hcid = G4SDManager::GetSDMpointer()->GetCollectionID(
                   G4String(SimphonyPhotonSD::kSDName) + "/"
                   + SimphonyPhotonSD::kHCName);
    }
    if (hcid < 0) return;

    // Always insert the collection, even when empty: the pseudo-SD is not
    // attached to a volume, so nothing else creates it, and consumers that
    // walk the hit-collection table expect a registered slot to be non-null.
    SEvt* sev = SEvt::Get_EGPU();
    if (!sev) nhit = 0;

    // The record buffer (present in DebugHeavy runs) provides each photon's
    // creation point for TG4PhotonHit.Start; without it Start stays (0,0,0,0),
    // which TG4PhotonHit documents as expected for offloaded photon tracking.
    const NP* recNP = sev->topfold ? sev->topfold->get(SComp::RECORD_) : nullptr;
    const sphoton* rec = recNP ? reinterpret_cast<const sphoton*>(recNP->bytes())
                               : nullptr;
    const int     maxRec     = recNP ? static_cast<int>(recNP->shape[1]) : 0;
    const int64_t nRecPhoton = recNP ? recNP->shape[0] : 0;

    const int pdgOptical
        = G4OpticalPhoton::OpticalPhotonDefinition()->GetPDGEncoding();

    auto* coll = new EDepSim::HitSurface::HitSurfaceCollection(
                     SimphonyPhotonSD::kSDName, SimphonyPhotonSD::kHCName);

    for (int64_t i = 0; i < nhit; ++i) {
        sphoton p;
        sev->getHit(p, static_cast<unsigned>(i));
        const int pIdx = static_cast<int>(p.get_index());

        // Creator process, mirroring what a CPU-tracked photon's HitSurface
        // carries (process type + subtype of the process that created the
        // arriving photon). Geant4 makes a NEW track at WLS re-emission whose
        // creator is OpWLS, so a GPU photon that was re-emitted in flight
        // (BULK_REEMIT in its history) reports OpWLS too, not its original
        // generation process.
        int creatorType = fUserDefined, creatorSubtype = 0;   // input/torch
        if (p.flagmask & BULK_REEMIT) {
            creatorType = fOptical;          creatorSubtype = fOpWLS;
        } else if (p.flagmask & CERENKOV) {
            creatorType = fElectromagnetic;  creatorSubtype = fCerenkov;
        } else if (p.flagmask & SCINTILLATION) {
            creatorType = fElectromagnetic;  creatorSubtype = fScintillation;
        }

        // Photon energy in G4 units from the wavelength (nm).
        const double e = (p.wavelength > 0)
                       ? h_Planck * c_light / (p.wavelength * nm) : 0.0;

        G4LorentzVector position(p.pos.x, p.pos.y, p.pos.z, p.time);  // mm, ns
        G4LorentzVector start(0, 0, 0, 0);
        if (rec && maxRec > 0 && pIdx >= 0 && pIdx < nRecPhoton) {
            const sphoton& r0 = rec[static_cast<int64_t>(maxRec) * pIdx];
            if (r0.flag() != 0)
                start.set(r0.pos.x, r0.pos.y, r0.pos.z, r0.time);
        }

        // Charged parent track that made the genstep; input/primary photons
        // have no charged parent (same convention as GPUPhotonTracks).
        const int primaryId = InputPhotonMode() ? -1 : RecoverTrackId(pIdx);
        coll->insert(new EDepSim::HitSurface(
            primaryId, e, position, start, pdgOptical,
            creatorType, creatorSubtype));
    }

    G4Event* ev = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
    G4HCofThisEvent* hce = ev ? ev->GetHCofThisEvent() : nullptr;
    if (!hce) {
        if (!ev) { delete coll; return; }
        hce = new G4HCofThisEvent(
                  G4SDManager::GetSDMpointer()->GetCollectionCapacity());
        ev->SetHCofThisEvent(hce);
    }
    hce->AddHitsCollection(hcid, coll);

    if (nhit > 0)
        std::cout << "[SimphonyPlugin] EndOfEvent: " << nhit
                  << " GPU hits -> PhotonDetectors['"
                  << SimphonyPhotonSD::kSDName << "']\n";
}

// Decode the terminating-flag of a GPU photon into a short fate string and a
// "detected?" boolean, mirroring the CPU-side fate labels.
static int decode_proc(unsigned flagmask)
{
    if      (flagmask & CERENKOV)      return 0;
    else if (flagmask & SCINTILLATION) return 2;
    else if (flagmask & TORCH)         return 3;   // primary "torch" photons
    return -1;
}

void SimphonyEventAction::SaveAllPhotonsAndTrajectories(const G4Event* event,
                                                        int64_t nphoton)
{
    SEvt* sev = SEvt::Get_EGPU();
    if (!sev || nphoton <= 0) return;

    TTree* trackTree = fRunAction->GetGPUTrackTree();
    TTree* stepTree  = fRunAction->GetGPUStepTree();
    if (!trackTree) return;

    // Pull buffer refs once.
    int& tEv  = GetGtEventId();   int& tPid = GetGtPhotonId();
    int& tTid = GetGtTrackId();   int& tPar = GetGtParentId();
    int& tPr  = GetGtProcess();   float& tWl = GetGtWavelength();
    unsigned& tFlag = GetGtFlag(); unsigned& tMask = GetGtFlagMask();
    int& tDet = GetGtDetected();  char* tFate = GetGtFate();
    int& tNs  = GetGtNStep();
    TLorentzVector& tStart = GetGtStartPos();
    TLorentzVector& tEnd   = GetGtEndPos();

    int& sEv  = GetGsEventId();   int& sPid = GetGsPhotonId();
    int& sStp = GetGsStep();      unsigned& sFlag = GetGsFlag();
    char* sAbbr = GetGsFlagAbbr();
    TLorentzVector& sPos = GetGsPos();

    tEv = sEv = event->GetEventID();

    // After simulate(), the gathered HOST arrays have been concat'd up to the
    // top fold. getPhoton() reads topfold; the record buffer is at
    // topfold->get("record"), shaped (num_photon, max_record, 4, 4) of float.
    const NP* phoNP = sev->getPhoton();
    const NP* recNP = sev->topfold ? sev->topfold->get(SComp::RECORD_) : nullptr;
    if (!phoNP || phoNP->shape[0] <= 0) return;

    const sphoton* pho = reinterpret_cast<const sphoton*>(phoNP->bytes());
    const int64_t  nP  = phoNP->shape[0];

    const sphoton* rec = recNP ? reinterpret_cast<const sphoton*>(recNP->bytes())
                               : nullptr;
    const int maxRec   = recNP ? static_cast<int>(recNP->shape[1]) : 0;
    const int64_t nRecPhoton = recNP ? recNP->shape[0] : 0;

    for (int64_t i = 0; i < nP; ++i) {
        const sphoton& p = pho[i];

        tPid = static_cast<int>(i);
        // GPU TrackId = the photon's own global index (one per photon; it
        // persists across WLS re-emission — Opticks keeps the same photon).
        // ParentId = the charged G4 track that produced the genstep, or -1 for
        // input/primary photons (which have no charged parent).
        tTid = static_cast<int>(p.get_index());
        // Input/primary photons have no charged parent track -> -1. For
        // genstep (scintillation) photons, recover the charged parent track.
        tPar = InputPhotonMode() ? -1
                                 : RecoverTrackId(static_cast<int>(p.get_index()));
        tPr  = decode_proc(p.flagmask);
        tWl  = p.wavelength;                 // nm
        tFlag = p.flag();                    // terminating flag
        tMask = p.flagmask;
        tDet = OpticksPhoton::IsDetectFlag(tFlag) ? 1 : 0;
        const char* ab = OpticksPhoton::Abbrev(tFlag);
        std::strncpy(tFate, ab ? ab : "??", 7);
        tFate[7] = '\0';

        // Final position of the photon (last record point if available,
        // else the photon's own end position).
        tEnd.SetXYZT(p.pos.x, p.pos.y, p.pos.z, p.time);

        // Walk this photon's record points (the full bounce-by-bounce path).
        int nstep = 0;
        if (rec && maxRec > 0 && i < nRecPhoton) {
            for (int s = 0; s < maxRec; ++s) {
                const sphoton& rp = rec[maxRec * i + s];
                // An unused record slot has flag()==0 (no step yet).
                if (s > 0 && rp.flag() == 0 &&
                    rp.pos.x == 0.f && rp.pos.y == 0.f && rp.pos.z == 0.f)
                    break;
                if (nstep == 0)
                    tStart.SetXYZT(rp.pos.x, rp.pos.y, rp.pos.z, rp.time);
                if (stepTree) {
                    sPid = static_cast<int>(i);
                    sStp = s;
                    sFlag = rp.flag();
                    const char* sab = OpticksPhoton::Abbrev(rp.flag());
                    std::strncpy(sAbbr, sab ? sab : "??", 7);
                    sPos.SetXYZT(rp.pos.x, rp.pos.y, rp.pos.z, rp.time);
                    stepTree->Fill();
                }
                nstep++;
            }
        }
        if (nstep == 0) {  // no record buffer -> start == end
            tStart = tEnd;
        }
        tNs = nstep;
        trackTree->Fill();
    }
}
