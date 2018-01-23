/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  Rafał Lalik <rafal.lalik@@ph.tum.de>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef PP35_LAMBDA_H
#define PP35_LAMBDA_H

#include "henergylosscorrpar.h"

#include "KAbstractAnalysis.h"
#include "KTifiniAnalysis.h"
#include "KBeamCalibration.h"

// constants
static const Float_t NULL_DATA = -10000.;
static const Float_t MASS_OK = 10000.;
static const Float_t ERR_NOT_IN_VERTEX = -10001.;
static const Float_t ERR_SAME_TRACKS = -10002.;
static const Float_t ERR_NO_LAMBDA = -10003.;
static const Float_t ERR_VERTEX_Z_MISSMATCH = -10004.;
static const Float_t ERR_MASS_OUT_OF_RANGE = -10005.;
static const Float_t ERR_SIM_VERTEX_MISSMATCH = -10006.;

//#define SHOWREJECTED

struct AnaDataSet
{
    // Kinematic properties
    Float_t fE;
    Float_t fM;
    Float_t fM_miss;
    Float_t fMt;
    Float_t fTheta;
    Float_t fCosTheta;
    Float_t fCosTheta_cms;
    Float_t fPhi;
    Float_t fP;
    Float_t fP_cms;
    Float_t fPt;
    Float_t fY;
    Float_t fY_cms;

    // track-reco properties
    Float_t fMTD;
    Float_t fVertDistX;
    Float_t fPVA;
    Float_t fPVA_miss;
    Float_t fPolX, fPolY, fPolZ;
    Float_t fZetaX, fZetaY, fZetaZ;
    Float_t fXf;
    Float_t fPolXi, fPolYi, fPolZi;
    Float_t fZetaXi, fZetaYi, fZetaZi;
    Float_t fXfi;

    // vertex properties
    Float_t fPrimVertexX;
    Float_t fPrimVertexY;
    Float_t fPrimVertexZ;
    Float_t fPrimVertexR;

    Float_t fFitVertexX;
    Float_t fFitVertexY;
    Float_t fFitVertexZ;

    Float_t fDecayVertexX;
    Float_t fDecayVertexY;
    Float_t fDecayVertexZ;
    Float_t fDecayVertexR;

    Float_t fEventVertexX;
    Float_t fEventVertexY;
    Float_t fEventVertexZ;

    Float_t fDVres;

    // daughter properties
    Float_t fEA, fEB;
    Float_t fMassA, fMassB;
    Float_t fMomA, fMomB;
    Float_t fMomAx, fMomAy, fMomAz, fMomBx, fMomBy, fMomBz;
    Float_t fMomA_cms, fMomB_cms;
    Float_t fMDCdEdxA, fMDCdEdxB;
    Float_t fAngleAB;
    Float_t fRelAngleA, fRelAngleB;
    Float_t fChiA, fChiB;
    Float_t fMetaMatchQA, fMetaMatchQB;
    Float_t fVertDistA, fVertDistB;
    Int_t fMultA, fMultB;
    Int_t fAccA, fAccB;

    // event stats
    Int_t fHadesTracksNum;
    Int_t fFwDetTracksNum;
    Int_t fIsFwDetData;
    Int_t fEventLCounter;
    Int_t fBestTry;
    Int_t fRealLambda;
    Int_t fPrimLambda;

    // other
    Int_t fGeantInfoNum;
    Float_t fGeantWeight;

    Float_t ret;
    Int_t fSortOrderMtd;
#ifdef SHOWREJECTED
    Int_t fRejected;
#endif

    // Wal data
    Float_t fWallT;
    Float_t fWallX;
    Float_t fWallY;
    Float_t fWallR;
    Float_t fWallCharge;
    Float_t fWallP;
    Float_t fWallPx;
    Float_t fWallPy;
    Float_t fWallPz;
    Float_t fWallBeta;
    Float_t fWallGamma;
    Int_t fWallClusterSize;
    Int_t fWallClustersNum;

    // Geant
    Float_t fGeaP, fGeaPx, fGeaPy, fGeaPz;
    Float_t fGeaTheta, fGeaPhi;
    Float_t fGeaPolX, fGeaPolY, fGeaPolZ;
    Float_t fGeaZetaX, fGeaZetaY, fGeaZetaZ;
    Float_t fGeaXf;
    Float_t fGeaAngleAB;

    Int_t fPVtype;

    void clear()
    {
        fHadesTracksNum = 0;
        fFwDetTracksNum = 0;
        fIsFwDetData = 0;

        fGeantWeight = 1.0;
        fEventVertexX = 0.0;
        fEventVertexY = 0.0;
        fEventVertexZ = 0.0;

        fEventLCounter = 0;

        fAccA = -1;
        fAccB = -1;

        fMultA = 0;
        fMultB = 0;
    }

    void init()
    {
        fE = NULL_DATA;
        fM = NULL_DATA;
        fM_miss = NULL_DATA;
        fMt = NULL_DATA;
        fTheta = NULL_DATA;
        fCosTheta = NULL_DATA;
        fCosTheta_cms = NULL_DATA;
        fPhi = NULL_DATA;
        fP = NULL_DATA;
        fP_cms = NULL_DATA;
        fPt = NULL_DATA;
        fY = NULL_DATA;
        fY_cms = NULL_DATA;
        fPolX = NULL_DATA;
        fPolY = NULL_DATA;
        fPolZ = NULL_DATA;
        fZetaX = NULL_DATA;
        fZetaY = NULL_DATA;
        fZetaZ = NULL_DATA;
        fXf = NULL_DATA;
        fPolXi = NULL_DATA;
        fPolYi = NULL_DATA;
        fPolZi = NULL_DATA;
        fZetaXi = NULL_DATA;
        fZetaYi = NULL_DATA;
        fZetaZi = NULL_DATA;
        fXfi = NULL_DATA;

        fMTD = NULL_DATA;
        fVertDistX = NULL_DATA;
        fPVA = NULL_DATA;
        fPVA_miss = NULL_DATA;

        fPrimVertexX = NULL_DATA;
        fPrimVertexY = NULL_DATA;
        fPrimVertexZ = NULL_DATA;
        fPrimVertexR = NULL_DATA;

        fFitVertexX = NULL_DATA;
        fFitVertexY = NULL_DATA;
        fFitVertexZ = NULL_DATA;

        fDecayVertexX = NULL_DATA;
        fDecayVertexY = NULL_DATA;
        fDecayVertexZ = NULL_DATA;
        fDecayVertexR = NULL_DATA;

        fDVres = NULL_DATA;

        fEA = fEB = NULL_DATA;
        fMassA = fMassB = NULL_DATA;
        fMomA = fMomB = NULL_DATA;
        fMomAx = fMomAy = fMomAz = fMomBx = fMomBy = fMomBz = NULL_DATA;
        fMomA_cms = fMomB_cms = NULL_DATA;
        fMDCdEdxA = fMDCdEdxB = NULL_DATA;
        fAngleAB = NULL_DATA;
        fRelAngleA = fRelAngleB = NULL_DATA;
        fChiA = fChiB = NULL_DATA;
        fMetaMatchQA = fMetaMatchQB = NULL_DATA;
        fVertDistA = fVertDistB = NULL_DATA;

        fGeantInfoNum = 0;
        fRealLambda = NULL_DATA;
        fPrimLambda = NULL_DATA;

        fWallT = 0.0;
        fWallX = 0.0;
        fWallY = 0.0;
        fWallR = 0.0;
        fWallCharge = 0.0;
        fWallP = 0.0;
        fWallPx = 0.0;
        fWallPy = 0.0;
        fWallPz = 0.0;
        fWallBeta = 0.0;
        fWallGamma = 0.0;
        fWallClusterSize = 0;
        fWallClustersNum = 0;

        fGeaP = fGeaPx = fGeaPy = fGeaPz = NULL_DATA;
        fGeaTheta = fGeaPhi = NULL_DATA;
        fGeaPolX = fGeaPolY = fGeaPolZ = NULL_DATA;
        fGeaZetaX = fGeaZetaY = fGeaZetaZ = NULL_DATA;
        fGeaXf = NULL_DATA;
        fGeaAngleAB = NULL_DATA;

        fPVtype = NULL_DATA;

#ifdef SHOWREJECTED
        fRejected = 0;
#endif /*SHOWREJECTED*/

        ret = 0.;
        fSortOrderMtd = 0;
    }
};

class pp45_Lambda : public KAbstractAnalysis
{
public:
    pp45_Lambda(const TString & analysisName, const TString & treeName);
    virtual bool analysis(Int_t event_num, HCategory * pcand, Int_t cand_size, HCategory * vcand, Int_t vect_size);

    void initAnalysis(KT::Experiment exp, KT::AnalysisType analysisType);
    void finalizeAnalysis();

    int Configure(int argc, char ** argv);
    void Usage() const;

protected:
    void configureTree(TTree * tree);
    void configureGraphicalCuts(KTrackInspector & cuts);

    AnaDataSet singlePairAnalysis(Int_t event_num, HCategory * pcand, int trackA_num, int trackB_num, bool quick_run = false);
    AnaDataSet singleFwDetPairAnalysis(Int_t event_num, HCategory * pcand, HCategory * vcand, int trackA_num, int trackB_num, bool quick_run = false);

    HEnergyLossCorrPar * eLossCorr;
    KBeamCalibration * beamCal;
    HGeomVector refBeamVector;

    // opts
    int flag_nosecvertcuts;
    int flag_elosscorr;
    int flag_nosigmas;
    int flag_useeventvertex;
    int flag_usewall;

    int par_Mtd;			// CLI
    int par_VertDistX;		// CLI
    int par_VertDistA;		// CLI
    int par_VertDistB;		// CLI
};

#endif // PP35_LAMBDA_H
