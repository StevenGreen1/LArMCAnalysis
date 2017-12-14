/**
 *  @file   MCAnalyser/MCAnalysis_module.cc
 *
 *  @brief  Analysis module for MC particles
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"

#include "TTree.h"                                       
#include "TFile.h"                                                                                                             
#include "art/Framework/Services/Optional/TFileService.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief MCAnalysis class
 */

class MCAnalysis : public art::EDAnalyzer 
{
public:
    /**
     *  @brief Constructor
     *
     *  @param pset
     */
    MCAnalysis(fhicl::ParameterSet const & p);

    /**
     *  @brief Destructor
     */
    virtual ~MCAnalysis();

    /**
     *  @brief Read information from fhicl
     */
    void reconfigure(fhicl::ParameterSet const & parameterSet);

    /**
     *  @brief Setup root tree for filling
     */
    void beginJob();

    /**
     *  @brief Save root tree 
     */
    void endJob();

    /**
     *  @brief Reset member variables for next MCParticle
     */
    void reset();

    /**
     *  @brief Extract event information 
     */
    void analyze(art::Event const & event);

private:
    /**
     *  @brief Is the point (x,y,z) in the LArTPC
     *
     *  @param x coordinate
     *  @param y coordinate
     *  @param z coordinate
     */ 
    bool IsInTPC(double x, double y, double z);

//    TFile           *m_pTFile;          ///< TFile for saving event information 
    TTree           *m_pTTree;          ///< TTree for saving event information

    int              m_eventNumber;     ///< Event number

    bool             m_isPrimary;       ///< Does the MCParticle have any parents
    bool             m_isBeam;          ///< Is the MCParticle from the beam, or is it a cosmic ray
    bool             m_entersTPC;       ///< Does the MCParticle enter the TPC
    int              m_pdg;             ///< PDG code for the MCParticle
    int              m_trackID;         ///< Geant4 track ID for the MCParticle
    int              m_parent;          ///< TrackID of parent MCParticle
    std::vector<int> m_daughters;       ///< TrackIDs of daughter MCParticles

    double           m_startEnergy;     ///< MCParticle energy at start
    double           m_startPX;         ///< MCParticle momentum along x direction at start
    double           m_startPY;         ///< MCParticle momentum along y direction at start
    double           m_startPZ;         ///< MCParticle momentum along z direction at start
    double           m_startT;          ///< MCParticle time at start
    double           m_startX;          ///< MCParticle x position at start 
    double           m_startY;          ///< MCParticle y position at start
    double           m_startZ;          ///< MCParticle z position at start

    double           m_endEnergy;       ///< MCParticle energy at end
    double           m_endPX;           ///< MCParticle momentum along x direction at end
    double           m_endPY;           ///< MCParticle momentum along y direction at end
    double           m_endPZ;           ///< MCParticle momentum along z direction at end
    double           m_endT;            ///< MCParticle time at end
    double           m_endX;            ///< MCParticle x position at end
    double           m_endY;            ///< MCParticle y position at end
    double           m_endZ;            ///< MCParticle z position at end

    double           m_centreX;         ///< ProtoDUNE x coordinate centre
    double           m_centreY;         ///< ProtoDUNE y coordinate centre
    double           m_centreZ;         ///< ProtoDUNE z coordinate centre
    double           m_widthX;          ///< ProtoDUNE x coordinate width
    double           m_widthY;          ///< ProtoDUNE y coordinate width
    double           m_widthZ;          ///< ProtoDUNE z coordinate width
    double           m_lowXTPC;         ///< Lowest x coordinate extent of ProtoDUNE TPC 
    double           m_highXTPC;        ///< Highest x coordinate extent of ProtoDUNE TPC 
    double           m_lowYTPC;         ///< Lowest y coordinate extent of ProtoDUNE TPC 
    double           m_highYTPC;        ///< Highest y coordinate extent of ProtoDUNE TPC 
    double           m_lowZTPC;         ///< Lowest z coordinate extent of ProtoDUNE TPC 
    double           m_highZTPC;        ///< Highest z coordinate extent of ProtoDUNE TPC 
};

//------------------------------------------------------------------------------------------------------------------------------------------

MCAnalysis::MCAnalysis(fhicl::ParameterSet const & p) : art::EDAnalyzer(p)
{
    this->reset();
    this->reconfigure(p);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCAnalysis::~MCAnalysis()
{                                                                                                                                 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCAnalysis::reconfigure(fhicl::ParameterSet const & p)
{
    m_centreX = p.get<double>("ProtoDUNE_CentreX",0.0);
    m_centreY = p.get<double>("ProtoDUNE_CentreY",303.749359131);
    m_centreZ = p.get<double>("ProtoDUNE_CentreZ",347.396240234);
    m_widthX = p.get<double>("ProtoDUNE_WidthX",759.439697266);
    m_widthY = p.get<double>("ProtoDUNE_WidthY",599.88873291);
    m_widthZ = p.get<double>("ProtoDUNE_WidthZ",695.780029297);

    m_lowXTPC = m_centreX - 0.5 * m_widthX;
    m_highXTPC = m_centreX + 0.5 * m_widthX;
    m_lowYTPC = m_centreY - 0.5 * m_widthY;
    m_highYTPC = m_centreY + 0.5 * m_widthY;
    m_lowZTPC = m_centreZ - 0.5 * m_widthZ;
    m_highZTPC = m_centreZ + 0.5 * m_widthZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCAnalysis::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
//    m_pTFile = new TFile("MCAnalysis.root", "RECREATE");
    m_pTTree = tfs->make<TTree>("MCAnalysisTree","MCAnalysisTree");

    m_pTTree->Branch("eventNumber", &m_eventNumber, "eventNumber/I");
    m_pTTree->Branch("isPrimary", &m_isPrimary, "isPrimary/O");
    m_pTTree->Branch("isBeam", &m_isBeam, "isBeam/O");
    m_pTTree->Branch("entersTPC", &m_entersTPC, "entersTPC/O");
    m_pTTree->Branch("pdg", &m_pdg, "pdg/I");
    m_pTTree->Branch("trackID", &m_trackID, "trackID/I");
    m_pTTree->Branch("parent", &m_parent, "parent/I");
    m_pTTree->Branch("daughters", &m_daughters);

    m_pTTree->Branch("startEnergy", &m_startEnergy, "startEnergy/D");
    m_pTTree->Branch("startPX", &m_startPX, "startPX/D");
    m_pTTree->Branch("startPY", &m_startPY, "startPY/D");
    m_pTTree->Branch("startPZ", &m_startPZ, "startPZ/D");
    m_pTTree->Branch("startT", &m_startT, "startT/D");
    m_pTTree->Branch("startX", &m_startX, "startX/D");
    m_pTTree->Branch("startY", &m_startY, "startY/D");
    m_pTTree->Branch("startZ", &m_startZ, "startZ/D");

    m_pTTree->Branch("endEnergy", &m_endEnergy, "endEnergy/D");
    m_pTTree->Branch("endPX", &m_endPX, "endPX/D");
    m_pTTree->Branch("endPY", &m_endPY, "endPY/D");
    m_pTTree->Branch("endPZ", &m_endPZ, "endPZ/D");
    m_pTTree->Branch("endT", &m_endT, "endT/D");
    m_pTTree->Branch("endX", &m_endX, "endX/D");
    m_pTTree->Branch("endY", &m_endY, "endY/D");
    m_pTTree->Branch("endZ", &m_endZ, "endZ/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCAnalysis::endJob()
{
//    m_pTFile->cd();
//    m_pTTree->Write("MCAnalysisTree");
//    m_pTTree->Write();
//    m_pTFile->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCAnalysis::reset()
{
    m_eventNumber = -999;
    m_isPrimary = false;
    m_isBeam = false;
    m_entersTPC = false;
    m_pdg = 0;
    m_trackID = -999;
    m_parent = -999;
    m_daughters.clear();

    m_startEnergy = -1.0;
    m_startPX = -1.0;
    m_startPY = -1.0;
    m_startPZ = -1.0;
    m_startT = -1.0;
    m_startX = -1.0;
    m_startY = -1.0;
    m_startZ = -1.0;

    m_endEnergy = -1.0;
    m_endPX = -1.0;
    m_endPY = -1.0;
    m_endPZ = -1.0;
    m_endT = -1.0;
    m_endX = -1.0;
    m_endY = -1.0;
    m_endZ = -1.0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCAnalysis::analyze(art::Event const &event)
{
    art::ServiceHandle<cheat::BackTracker> backTracker;
    const std::string producerName("largeant");
    art::Handle< std::vector<simb::MCParticle> > mcParticles;
    event.getByLabel(producerName, mcParticles);

    for (unsigned int mcPart = 0; mcPart < mcParticles->size(); mcPart++)
    {
        const art::Ptr<simb::MCParticle> mcParticle(mcParticles, mcPart);
        this->reset();

        m_eventNumber = event.id().event();
        m_pdg = mcParticle->PdgCode();
        m_trackID = mcParticle->TrackId();
        m_parent = mcParticle->Mother();

        if ("primary" == backTracker->TrackIDToParticle(m_trackID)->Process())
        {
            m_isPrimary = true;
        }

        if (simb::kSingleParticle == backTracker->TrackIDToMCTruth(m_trackID)->Origin())
        {
            m_isBeam = true;
        }

        for (int daughter = 0; daughter < mcParticle->NumberDaughters(); daughter++)
        {
            m_daughters.push_back(mcParticle->Daughter(daughter));
        }

        m_startEnergy = mcParticle->E();
        m_startPX = mcParticle->Px();
        m_startPY = mcParticle->Py();
        m_startPZ = mcParticle->Pz();
        m_startT = mcParticle->T();
        m_startX = mcParticle->Vx();
        m_startY = mcParticle->Vy();
        m_startZ = mcParticle->Vz();

        m_endEnergy = mcParticle->EndE();
        m_endPX = mcParticle->EndPx();
        m_endPY = mcParticle->EndPy();
        m_endPZ = mcParticle->EndPz();
        m_endT = mcParticle->EndT();
        m_endX = mcParticle->EndX();
        m_endY = mcParticle->EndY();
        m_endZ = mcParticle->EndZ();

        for (unsigned int step = 0; step < mcParticle->NumberTrajectoryPoints(); step++)
        {
            double stepX(mcParticle->Vx(step));
            double stepY(mcParticle->Vy(step));
            double stepZ(mcParticle->Vz(step));

            bool inTPC(this->IsInTPC(stepX, stepY, stepZ));

            if (inTPC)
            {
                m_entersTPC = true;
                break;
            }
        }

        m_pTTree->Fill();
    }

/*
    // Get MC Truth
std::cout << "GENERATOR" << std::endl;
    const std::string producerNameGen("generator");
    art::Handle< std::vector<simb::MCTruth> > mcTruth;
    e.getByLabel(producerNameGen, mcTruth);
    std::cout << "mcTruth->size() = " << mcTruth->size() << std::endl;
    const art::Ptr<simb::MCTruth> pMCTruth(mcTruth, 0);
    const int nMCTruthParticles(pMCTruth->NParticles());
    simb::Origin_t nMCTruthOrigin(pMCTruth->Origin());
    std::cout << "MC Truth" << std::endl;
    std::cout << "nMCTruthParticles = " << nMCTruthParticles << std::endl;
    std::cout << "nMCTruthOrigin = " << nMCTruthOrigin << std::endl;

    std::cout << "TrackId | PDG | (E,  Px, Py, Pz) | Start(t,x,y,z) | End(t,x,y,z) | NumberTrajectoryPoints | nDaughters [Daughter TrackIDs] | Mother" << std::endl;

    std::cout.precision(3);

    for (int mcPart = 0; mcPart < nMCTruthParticles; mcPart++)
    {
        const simb::MCParticle mcParticle(pMCTruth->GetParticle(mcPart));
        std::cout << mcParticle.TrackId() << " | " << mcParticle.PdgCode();
        std::cout << std::scientific;
        std::cout << " | (" << mcParticle.E() << ", " << mcParticle.Px() << ", " << mcParticle.Py() << ", " << mcParticle.Pz() << " )";
        std::cout << " | ( " << mcParticle.T() << ", " << mcParticle.Vx() << ", " << mcParticle.Vy() << ", " << mcParticle.Vz() << " )";
        std::cout << " | ( " << mcParticle.EndT() << ", " << mcParticle.EndX() << ", " << mcParticle.EndY() << ", " << mcParticle.EndZ() << " )";
        std::cout << " | " << mcParticle.NumberTrajectoryPoints() << " | " << mcParticle.NumberDaughters() << " [ ";

        for (int daughter = 0; daughter < mcParticle.NumberDaughters(); daughter++)
        {
            std::cout << mcParticle.Daughter(daughter);
            if (daughter != mcParticle.NumberDaughters() - 1)
                std::cout << ", ";
        }
        std::cout << " ] | " << mcParticle.Mother();
        std::cout << std::endl;
    }

std::cout << "COSMICGENERATOR" << std::endl;

    const std::string producerNameCRGen("cosmicgenerator");
    art::Handle< std::vector<simb::MCTruth> > mcCRTruth;
    e.getByLabel(producerNameCRGen, mcCRTruth);
    std::cout << "mcCRTruth->size() = " << mcCRTruth->size() << std::endl;
    const art::Ptr<simb::MCTruth> pMCCRTruth(mcCRTruth, 0);
    const int nMCCRTruthParticles(pMCCRTruth->NParticles());
    simb::Origin_t nMCCRTruthOrigin(pMCCRTruth->Origin());
    std::cout << "MC CR Truth" << std::endl;
    std::cout << "nMCCRTruthParticles = " << nMCCRTruthParticles << std::endl;
    std::cout << "nMCCRTruthOrigin = " << nMCCRTruthOrigin << std::endl;

    std::cout << "TrackId | PDG | (E,  Px, Py, Pz) | Start(t,x,y,z) | End(t,x,y,z) | NumberTrajectoryPoints | nDaughters [Daughter TrackIDs] | Mother" << std::endl;

    std::cout.precision(3);

    for (int mcPart = 0; mcPart < nMCCRTruthParticles; mcPart++)
    {
        const simb::MCParticle mcParticle(pMCCRTruth->GetParticle(mcPart));
        std::cout << mcParticle.TrackId() << " | " << mcParticle.PdgCode();
        std::cout << std::scientific;
        std::cout << " | (" << mcParticle.E() << ", " << mcParticle.Px() << ", " << mcParticle.Py() << ", " << mcParticle.Pz() << " )";
        std::cout << " | ( " << mcParticle.T() << ", " << mcParticle.Vx() << ", " << mcParticle.Vy() << ", " << mcParticle.Vz() << " )";
        std::cout << " | ( " << mcParticle.EndT() << ", " << mcParticle.EndX() << ", " << mcParticle.EndY() << ", " << mcParticle.EndZ() << " )";
        std::cout << " | " << mcParticle.NumberTrajectoryPoints() << " | " << mcParticle.NumberDaughters() << " [ ";

        for (int daughter = 0; daughter < mcParticle.NumberDaughters(); daughter++)
        {
            std::cout << mcParticle.Daughter(daughter);
            if (daughter != mcParticle.NumberDaughters() - 1)
                std::cout << ", ";
        }
        std::cout << " ] | " << mcParticle.Mother();
        std::cout << std::endl;
    }

std::cout << "GEANT4" << std::endl;

    // Get MC Particles
    const std::string producerName("largeant");
    art::Handle< std::vector<simb::MCParticle> > mcParticles;
    e.getByLabel(producerName, mcParticles);
    const int nMCParticles(mcParticles->size());

    std::cout << "TrackId | PDG | (E,  Px, Py, Pz) | Start(t,x,y,z) | End(t,x,y,z) | NumberTrajectoryPoints | nDaughters [Daughter TrackIDs] | Mother" << std::endl;

    std::cout.precision(3);

    for (unsigned int mcPart = 0; mcPart < mcParticles->size(); mcPart++)
    {
        const art::Ptr<simb::MCParticle> mcParticle(mcParticles, mcPart);
        std::cout << mcParticle->TrackId() << " | " << mcParticle->PdgCode();
        std::cout << std::scientific;
        std::cout << " | (" << mcParticle->E() << ", " << mcParticle->Px() << ", " << mcParticle->Py() << ", " << mcParticle->Pz() << " )";
        std::cout << " | ( " << mcParticle->T() << ", " << mcParticle->Vx() << ", " << mcParticle->Vy() << ", " << mcParticle->Vz() << " )";
        std::cout << " | ( " << mcParticle->EndT() << ", " << mcParticle->EndX() << ", " << mcParticle->EndY() << ", " << mcParticle->EndZ() << " )";
        std::cout << " | " << mcParticle->NumberTrajectoryPoints() << " | " << mcParticle->NumberDaughters() << " [ ";

        for (int daughter = 0; daughter < mcParticle->NumberDaughters(); daughter++)
        {
            std::cout << mcParticle->Daughter(daughter);
            if (daughter != mcParticle->NumberDaughters() - 1)
                std::cout << ", ";
        }
        std::cout << " ] | " << mcParticle->Mother();
        std::cout << std::endl;
    }

    std::cout << "Found a total of " << nMCParticles << " MC Particles from the " << producerName << std::endl;
*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MCAnalysis::IsInTPC(double x, double y, double z)
{
    if (m_lowXTPC < x && x < m_highXTPC && m_lowYTPC < y && y < m_highYTPC && m_lowZTPC < z && z < m_highZTPC)
        return true;
    else
        return false;
}

DEFINE_ART_MODULE(MCAnalysis)
