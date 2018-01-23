/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2011  Rafał Lalik <rafal.lalik@ph.tum.de>

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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <getopt.h>

#include "hphysicsconstants.h"
#include "hparticlecandsim.h"

#include "KTools.h"
#include "KTrackInspector.h"
#include "KBeamCalibration.h"

#include "pp45_Lambda.h"

#ifdef SHOWREJECTED
//int fRejected;
#endif /*SHOWREJECTED*/

// double TARGET_X_MEAN = 2.0;
// double TARGET_Y_MEAN = 0.0;
// double TARGET_Z_MIN = -70.0;
// double TARGET_Z_MAX = 80.0;

static const float E_kin_beam       = 4500.0;
static const float E_kin_target     = 0.0;
static const float E_total_beam     = E_kin_beam + HPhysicsConstants::mass(14);
static const float E_total_target   = E_kin_target + HPhysicsConstants::mass(14);
static const float pz_beam          = sqrt(E_total_beam*E_total_beam-HPhysicsConstants::mass(14)*HPhysicsConstants::mass(14));

static const float LAMBDA_MASS = 1115.9;
static const float XI_MASS = 1115.9;

// static const float D2R = TMath::DegToRad();
static const float R2D = TMath::RadToDeg();

static AnaDataSet g_ads;

float calcAngleVar(HGeomVector & v1, HGeomVector & v2)
{
    TVector3 _v1; _v1.SetXYZ(v1.X(), v1.Y(), v1.Z());
    TVector3 _v2; _v2.SetXYZ(v2.X(), v2.Y(), v2.Z());

    return _v1.Angle(_v2);
}

float calcAngleVar(HGeomVector & v1, TLorentzVector & v2)
{
    TVector3 _v1; _v1.SetXYZ(v1.X(), v1.Y(), v1.Z());
    TVector3 _v2; _v2.SetXYZ(v2.Px(), v2.Py(), v2.Pz());

    return _v1.Angle(_v2);
}

pp45_Lambda::pp45_Lambda(const TString& analysisName, const TString& treeName) : KAbstractAnalysis(analysisName, treeName),
    flag_nosecvertcuts(0), flag_elosscorr(0), flag_nosigmas(0), flag_useeventvertex(0), par_Mtd(10), par_VertDistX(45), par_VertDistA(5), par_VertDistB(15)
{
//     setGoodEventSelector(Particle::kGoodTRIGGER |
//         Particle::kGoodVertexClust |
// //         Particle::kGoodVertexCand |
//         Particle::kGoodSTART |
//         Particle::kNoPileUpSTART);

    setGoodEventSelector(0);

    eLossCorr = new HEnergyLossCorrPar("eLossCorr", "eLossCorr", "eLossCorr");
    eLossCorr->setDefaultPar("jan04");
}

bool pp45_Lambda::analysis(Int_t event_num, HCategory * pcand, Int_t cand_size, HCategory * vcand, Int_t vect_size)
{
    printf("%d  %d\n", cand_size, vect_size);
    if ( (cand_size + vect_size) < 2 )
        return false;

    int A_PID = KT::p;    // proton
    int B_PID = KT::pim;    // pi-

    g_ads.fHadesTracksNum = cand_size;
    g_ads.fFwDetTracksNum = vect_size;
    g_ads.fIsFwDetData = 0;

    g_ads.fGeantWeight = 0;

//     g_ads.fEventVertexX = event->getVertexX();
//     g_ads.fEventVertexY = event->getVertexY();
//     g_ads.fEventVertexZ = event->getVertexZ();

    g_ads.fEventLCounter = 0;

    float mass_diff = MASS_OK;
    g_ads.fAccA = -1;
    g_ads.fAccB = -1;

    g_ads.fMultA = 0;
    g_ads.fMultB = 0;

    size_t combo_cnt = 0;

    for(int i = 0; i < cand_size; ++i)
    {
        printf("PID_a = %d\n", hades_tracks[i].pid[A_PID][KT::Charge]);
        if (!hades_tracks[i].pid[A_PID][KT::Charge])
            continue;

        // count A particle
        ++g_ads.fMultA;

        for(int j = 0; j < cand_size; ++j)
        {
            printf("PID_b = %d\n", hades_tracks[i].pid[B_PID][KT::Charge]);
            if (!hades_tracks[j].pid[B_PID][KT::Charge])
                continue;
printf("pair accepted!\n");
            // count B particle, to avoid double counting, only when A is 1 so will not count B if no A available
            if (g_ads.fMultA == 1)
                ++g_ads.fMultB;

            AnaDataSet ads_ret = singlePairAnalysis(event_num, pcand, i, j, true);

            if (ads_ret.ret < NULL_DATA)
                continue;

            ++combo_cnt;

            if (fabs(ads_ret.ret) < mass_diff)
            {
                mass_diff = fabs(ads_ret.ret);
                g_ads.fAccA = i;
                g_ads.fAccB = j;
            }
        }
    }

    AnaDataSet * ads_arr = new AnaDataSet[combo_cnt];
    size_t combo_pos = 0;
    std::vector<float> mtd_list;

    for(int i = 0; i < /*fMultA*/g_ads.fHadesTracksNum; ++i)
    {
        if (!hades_tracks[i].pid[A_PID])
            continue;

        for(int j = 0; j < /*fMultB*/g_ads.fHadesTracksNum; ++j)
        {
            if (!hades_tracks[j].pid[B_PID])
                continue;

            AnaDataSet ads_ret = singlePairAnalysis(event_num, pcand, i, j);
PR(ads_ret.ret);

            if (ads_ret.ret == ERR_MASS_OUT_OF_RANGE)
                continue;

            if (ads_ret.ret == ERR_SAME_TRACKS)
                continue;

            if (ads_ret.ret == ERR_NO_LAMBDA)
                continue;

            if (ads_ret.ret <= NULL_DATA)
                continue;

            if (g_ads.fAccA == i and g_ads.fAccB == j)
                ads_ret.fBestTry = 1;
            else
                ads_ret.fBestTry = 0;
PR("pass");
            mtd_list.push_back(ads_ret.fMTD);

            ads_arr[combo_pos] = ads_ret;
            ++combo_pos;
        }
    }

    std::sort(mtd_list.begin(), mtd_list.end());

    for (size_t i = 0; i < combo_pos; ++i)
    {
        for (size_t j = 0; j < combo_pos; ++j)
        {
            if (ads_arr[i].fMTD == mtd_list[j])
                ads_arr[i].fSortOrderMtd = j;
        }
    }

    for (size_t i = 0; i < combo_pos; ++i)
    {//PR(222);
        g_ads = ads_arr[i];
        getTree()->Fill();
    }

    delete [] ads_arr;
    return true;
    return false;
}

AnaDataSet pp45_Lambda::singlePairAnalysis(int event_num, HCategory * pcand, int trackA_num, int trackB_num, bool quick_run)
{
    HGeomVector beamVector; // FIXME
//     if (analysisType == KT::Exp)
//     {
//         beamVector = beamCal->calculateBeamOffset(event->getRunId());
//     } else {
        beamVector = refBeamVector;
//     }

    TVector3 p_beam_vec(beamVector.X(), beamVector.Y(), 3000.0);
    p_beam_vec.SetMag(pz_beam);

    const TLorentzVector Vec_pp35_beam      = TLorentzVector(p_beam_vec.X(), p_beam_vec.Y(), p_beam_vec.Z(), E_total_beam);
    const TLorentzVector Vec_pp35_target    = TLorentzVector(0.0, 0.0, 0.0, E_total_target);
    const TLorentzVector Vec_pp35_sum       = Vec_pp35_beam + Vec_pp35_target;
    const float cmrap                       = Vec_pp35_sum.Rapidity();

    TLorentzVector vec_beam_cms             = Vec_pp35_beam;
    vec_beam_cms.Boost(-Vec_pp35_sum.BoostVector());


    AnaDataSet ads = g_ads;
    ads.init();
//     ads.fGeantWeight = pcand->getGeantGenweight(); FIXME

    HGeomVector baseA, dirA, baseB, dirB;
    HGeomVector dirMother, PrimVertexMother;

    int A_PID = KT::p;    // proton
    int B_PID = KT::pim;    // pi-

    // Get the tracks and calculate the direction and base vectors
    if (trackA_num == trackB_num)
    {
    #ifdef SHOWREJECTED
        ads.fRejected = ERR_SAME_TRACKS;
    #else
        ads.ret = ERR_SAME_TRACKS;
        return ads;
    #endif /*SHOWREJECTED*/
    }

    HParticleCandSim trackA = *(HParticleCandSim*)pcand->getObject(trackA_num);
    HParticleCandSim trackB = *(HParticleCandSim*)pcand->getObject(trackB_num);

    trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));
    trackB.calc4vectorProperties(HPhysicsConstants::mass(B_PID));

    KTifini::CalcSegVector(trackA.getZ(), trackA.getR(), trackA.getPhi(), trackA.getTheta(), baseA, dirA);
    KTifini::CalcSegVector(trackB.getZ(), trackB.getR(), trackB.getPhi(), trackB.getTheta(), baseB, dirB);

    float momentum_A_corr = 0;
    float momentum_B_corr = 0;

    ads.fMomA = trackA.getMomentum();
    ads.fMomB = trackB.getMomentum();

    if (flag_elosscorr and analysisType == KT::Exp)
    {
        // with corr
        momentum_A_corr = eLossCorr->getCorrMom(A_PID, ads.fMomA, trackA.getTheta());
        momentum_B_corr = eLossCorr->getCorrMom(B_PID, ads.fMomB, trackB.getTheta());
    }
    else
    {
        // no corr
        momentum_A_corr = ads.fMomA;
        momentum_B_corr = ads.fMomB;
        momentum_A_corr = eLossCorr->getCorrMom(A_PID, ads.fMomA, trackA.getTheta());
        momentum_B_corr = eLossCorr->getCorrMom(B_PID, ads.fMomB, trackB.getTheta());

    }

    trackA.setMomentum(momentum_A_corr);
    trackB.setMomentum(momentum_B_corr);

    trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));
    trackB.calc4vectorProperties(HPhysicsConstants::mass(B_PID));

    TLorentzVector trackAB = trackA + trackB;

    ads.fMomAx = trackA.Px();
    ads.fMomAy = trackA.Py();
    ads.fMomAz = trackA.Pz();

    ads.fMomBx = trackB.Px();
    ads.fMomBy = trackB.Py();
    ads.fMomBz = trackB.Pz();

    ads.fM = trackAB.M();
    ads.fE = trackAB.E();
    ads.fEA = trackA.E();
    ads.fEB = trackB.E();

//    printf("L: px=%f  py=%f  pz=%f  e=%f\n", lambda.Px(), lambda.Py(), lambda. Pz(), ads.fE);
//    printf("beta=%f p_p=%f p_l=%f\n", beta_vec.Mag(), trackA.P(), l_lambda.P());
//    PR(ads.fPol);

    // I do not need so many data!
    if (ads.fM > 1200)
    {
    #ifdef SHOWREJECTED
        ads.fRejected = ERR_MASS_OUT_OF_RANGE;
    #else
        ads.ret = ERR_MASS_OUT_OF_RANGE;
        return ads;
    #endif /*SHOWREJECTED*/
    }
//PR(ads.fM);
    if (quick_run)
    {
        ads.ret = (LAMBDA_MASS - ads.fM);
        return ads;// * fMinTrackDist;
    }

    ads.fMTD = KTifini::calculateMinimumDistance(baseA, dirA, baseB, dirB);       // minimum distance between the two tracks
//     if (quick_run) printf(" mtd: %f\n", fMinTrackDist);

    float GeantxVertexA    = 0;
    float GeantyVertexA    = 0;
    float GeantzVertexA    = 0;
    float GeantxVertexB    = 0;
    float GeantyVertexB    = 0;
    float GeantzVertexB    = 0;

    // extra checks for the simulation analysis
    if (analysisType == KT::Sim)
    {
        TLorentzVector geaA; geaA.SetXYZM(trackA.getGeantxMom(), trackA.getGeantyMom(), trackA.getGeantzMom(), HPhysicsConstants::mass(A_PID));
        TLorentzVector geaB; geaB.SetXYZM(trackB.getGeantxMom(), trackB.getGeantyMom(), trackB.getGeantzMom(), HPhysicsConstants::mass(B_PID));
        TLorentzVector geaAB = geaA + geaB;
        ads.fGeaP = geaAB.P();
        ads.fGeaPx = geaAB.Px();
        ads.fGeaPy = geaAB.Py();
        ads.fGeaPz = geaAB.Pz();
        ads.fGeaTheta = geaAB.Theta() * R2D;
        ads.fGeaPhi = geaAB.Phi() * R2D;
        ads.fGeaAngleAB = geaA.Angle(geaB.Vect());

        TLorentzVector geaAB_cms = geaAB;
        geaAB_cms.Boost(-Vec_pp35_sum.BoostVector());
        ads.fGeaXf = fabs(geaAB_cms.Pz()/vec_beam_cms.Pz());

        // find plane normal
        TVector3 beam = Vec_pp35_beam.Vect();
        TVector3 lambda = geaAB.Vect();

        TVector3 n_x = beam.Cross(lambda);
        n_x *= (1.0/n_x.Mag());

        TVector3 n_z = lambda;
        n_z *= (1.0/n_z.Mag());

        TVector3 n_y = n_z.Cross(n_x);

        // lambda boost
        TLorentzVector l_lambda = geaAB;

        // proton -boost
        TVector3 beta_vec = l_lambda.BoostVector();
        TLorentzVector l_proton = geaA;
        l_proton.Boost(-beta_vec);

        TVector3 p_proton = l_proton.Vect();
        ads.fGeaPolX = p_proton.Dot(n_x) / (p_proton.Mag() * n_x.Mag());
        ads.fGeaPolY = p_proton.Dot(n_y) / (p_proton.Mag() * n_y.Mag());
        ads.fGeaPolZ = p_proton.Dot(n_z) / (p_proton.Mag() * n_z.Mag());
        ads.fGeaZetaX = TMath::ACos(ads.fGeaPolX) * R2D;
        ads.fGeaZetaY = TMath::ACos(ads.fGeaPolY) * R2D;
        ads.fGeaZetaZ = TMath::ACos(ads.fGeaPolZ) * R2D;


        GeantxVertexA    = trackA.getGeantxVertex();
        GeantyVertexA    = trackA.getGeantyVertex();
        GeantzVertexA    = trackA.getGeantzVertex();
        GeantxVertexB    = trackB.getGeantxVertex();
        GeantyVertexB    = trackB.getGeantyVertex();
        GeantzVertexB    = trackB.getGeantzVertex();

        // the simulated vertex of particleA and particleB has to be the same
        ads.fRealLambda = (GeantxVertexA == GeantxVertexB && GeantyVertexA == GeantyVertexB && GeantzVertexA == GeantzVertexB);

        int GeantPIDA = trackA.getGeantPID();
        int GeantPIDB = trackB.getGeantPID();

        int GeantPIDAMother = trackA.getGeantParentPID();
        int GeantPIDBMother = trackB.getGeantParentPID();

        int GeantPIDAGparent = trackA.getGeantGrandParentPID();
        int GeantPIDBGparent = trackB.getGeantGrandParentPID();

//         ads.fGeantInfoNum = event->getGeantInfoNum(); FIXME

        // FIXME
//         if (!wasLambda or (flag_nosigmas and wasSigma))
//         {
// #ifdef SHOWREJECTED
//             ads.fRejected = ERR_NO_LAMBDA;
// #else
//             ads.ret = ERR_NO_LAMBDA;
// //            return ads;
// #endif /*SHOWREJECTED*/
//         }
    }

//     if (quick_run) printf(" : %f, : %f\n", trackA.P(), trackB.P());

//             float thetaA = trackA.getTheta();
//             float thetaB = trackB.getTheta();

//             if (fMomA > 0.0 and fMomB > 0.0)
//             {
//                 float momAscale = momentum_A_corr / fMomA;
//                 float momBscale = momentum_B_corr / fMomB;
//                 trackA.SetXYZM(trackA.Px()*momAscale, trackA.Py()*momAscale, trackA.Pz() * momAscale, trackA.M());
//                 trackB.SetXYZM(trackB.Px()*momBscale, trackB.Py()*momBscale, trackB.Pz() * momBscale, trackB.M());
//             }

    ads.fChiA = trackA.getChi2();
    ads.fChiB = trackB.getChi2();

    ads.fMetaMatchQA = trackA.getMetaMatchQuality();
    ads.fMetaMatchQB = trackB.getMetaMatchQuality();

    ads.fMassA = trackA.M();
    ads.fMassB = trackB.M();

    ads.fMDCdEdxA = trackA.getMdcdEdx();
    ads.fMDCdEdxB = trackB.getMdcdEdx();

    dirMother.setXYZ(trackAB.X(), trackAB.Y(), trackAB.Z());    // direction vector of the mother particle
    HGeomVector DecayVertex = KTifini::calcVertexAnalytical(baseA, dirA, baseB, dirB);       // vertex of the two tracks

    ads.fDecayVertexX = DecayVertex.getX();
    ads.fDecayVertexY = DecayVertex.getY();
    ads.fDecayVertexZ = DecayVertex.getZ();
    ads.fDecayVertexR = TMath::Sqrt(ads.fDecayVertexX*ads.fDecayVertexX + ads.fDecayVertexY*ads.fDecayVertexY);
    if (analysisType == KT::Sim)
    {
        ads.fDVres = TMath::Sqrt(
                TMath::Power(ads.fDecayVertexX - GeantxVertexA, 2) +
                TMath::Power(ads.fDecayVertexY - GeantyVertexA, 2) +
                TMath::Power(ads.fDecayVertexZ - GeantzVertexA, 2)
        );

    }


    hgvpair vvectors;
    if (flag_useeventvertex)
    {
        PrimVertexMother.setXYZ(ads.fEventVertexX, ads.fEventVertexY, ads.fEventVertexZ);
        //PrimVertexMother += beamVector;
    }
    else
    {
//         FIXME
//         PrimVertexMother = KTifini::calcPrimVertex_Track_Mother(event, beamVector, DecayVertex, dirMother, trackA_num, trackB_num, ads.fPVtype);
    }

    ads.fPrimVertexX = (float)PrimVertexMother.getX();
    ads.fPrimVertexY = (float)PrimVertexMother.getY();
    ads.fPrimVertexZ = (float)PrimVertexMother.getZ();
    ads.fPrimVertexR = TMath::Sqrt(ads.fPrimVertexX*ads.fPrimVertexX + ads.fPrimVertexY*ads.fPrimVertexY);

    if ((ads.fDecayVertexZ - ads.fPrimVertexZ) < 0)
    {
#ifdef SHOWREJECTED
        ads.fRejected = ERR_VERTEX_Z_MISSMATCH;
#else
        ads.ret = ERR_VERTEX_Z_MISSMATCH;
        return ads;
#endif /*SHOWREJECTED*/
    }

    ads.fFitVertexX = (float)vvectors.second.getX();
    ads.fFitVertexY = (float)vvectors.second.getY();
    ads.fFitVertexZ = (float)vvectors.second.getZ();

    HGeomVector v1 = DecayVertex - PrimVertexMother;
//    HGeomVector v2 = DecayVertex - trackAB;
    //vvectors.second;

    ads.fVertDistX = (v1).length();

    TVector3 _v1; _v1.SetXYZ(v1.X(), v1.Y(), v1.Z());
    TVector3 _v2; _v2.SetXYZ(trackAB.Px(), trackAB.Py(), trackAB.Pz());

    ads.fPVA = calcAngleVar(v1, trackAB);

    ads.fVertDistA = KTifini::calculateMinimumDistanceStraightToPoint(baseA,dirA,PrimVertexMother);
    ads.fVertDistB = KTifini::calculateMinimumDistanceStraightToPoint(baseB,dirB,PrimVertexMother);

    // HERE we assign PV-SV vactor for Lambda

//    TVector3 vLpsv = _v1;
//    vLpsv.SetMag(trackAB.P());

//    TParticleTrack Lpsv;
//    Lpsv.setMomentum(vLpsv.Mag());
//    Lpsv.setTheta(vLpsv.Theta() * R2D);
//    Lpsv.setPhi(vLpsv.Phi() * R2D);
//    Lpsv.calc4vectorProperties(HPhysicsConstants::mass(18));

    TLorentzVector lambdaAB = trackAB;
//    TLorentzVector lambdaAB = Lpsv;

    // find plane normal
    TVector3 beam = Vec_pp35_beam.Vect();
    TVector3 lambda = lambdaAB.Vect();
//    lambda = Lpsv.Vect();

    TVector3 n_x = beam.Cross(lambda);
    n_x *= (1.0/n_x.Mag());

    TVector3 n_z = lambda;
    n_z *= (1.0/n_z.Mag());

    TVector3 n_y = n_z.Cross(n_x);

    // lambda boost, m from IM
    {
        TLorentzVector l_lambda = lambdaAB;

        // proton -boost
        TVector3 beta_vec = l_lambda.BoostVector();
        TLorentzVector l_proton = trackA;
        l_proton.Boost(-beta_vec);

//        TLorentzVector try_la = l_lambda; try_la.Boost(-beta_vec); printf("P after boost = %f\n", try_la.P());

        TVector3 p_proton = l_proton.Vect();
        ads.fPolX = p_proton.Dot(n_x) / (p_proton.Mag() * n_x.Mag());
        ads.fPolY = p_proton.Dot(n_y) / (p_proton.Mag() * n_y.Mag());
        ads.fPolZ = p_proton.Dot(n_z) / (p_proton.Mag() * n_z.Mag());

        ads.fZetaX = TMath::ACos(ads.fPolX) * R2D;
        ads.fZetaY = TMath::ACos(ads.fPolY) * R2D;
        ads.fZetaZ = TMath::ACos(ads.fPolZ) * R2D;
    }
    // lambda boost, m from PDG
    {
        TLorentzVector lpc;
        lpc.SetRho(lambdaAB.P());
        lpc.SetTheta(lambdaAB.Theta());
        lpc.SetPhi(lambdaAB.Phi());
        lpc.SetE(sqrt(HPhysicsConstants::mass(18)*HPhysicsConstants::mass(18) + lambdaAB.P()*lambdaAB.P()));

        TLorentzVector l_lambda = lpc;

        // proton -boost
        TVector3 beta_vec = l_lambda.BoostVector();
        TLorentzVector l_proton = trackA;
        l_proton.Boost(-beta_vec);

        TVector3 p_proton = l_proton.Vect();
        ads.fPolXi = p_proton.Dot(n_x) / (p_proton.Mag() * n_x.Mag());
        ads.fPolYi = p_proton.Dot(n_y) / (p_proton.Mag() * n_y.Mag());
        ads.fPolZi = p_proton.Dot(n_z) / (p_proton.Mag() * n_z.Mag());

        ads.fZetaXi = TMath::ACos(ads.fPolXi) * R2D;
        ads.fZetaYi = TMath::ACos(ads.fPolYi) * R2D;
        ads.fZetaZi = TMath::ACos(ads.fPolZi) * R2D;
    }

//     double dist2 = pow(PrimVertexMother.getX() - beamVector.getX(), 2.0) +
//         pow(PrimVertexMother.getY() - beamVector.getY(), 2.0);

//     if ( !(dist2 < 100.0 and PrimVertexMother.getZ() < 0.0 and PrimVertexMother.getZ() > -90.0) )
//         return ERR_NOT_IN_VERTEX;
    if ( !(ads.fPrimVertexR < 10.0 and PrimVertexMother.getZ() < 0.0 and PrimVertexMother.getZ() > -90.0) )
    {
#ifdef SHOWREJECTED
        ads.fRejected = ERR_NOT_IN_VERTEX;
#else
        ads.ret = ERR_NOT_IN_VERTEX;
        return ads;
#endif /*SHOWREJECTED*/
    }

    ads.fP          = trackAB.P();              // Momentum
    ads.fPt         = trackAB.Pt();             // Transverse momentum
    ads.fY          = trackAB.Rapidity();       // Rapidity
    ads.fAngleAB    = trackA.Angle(trackB.Vect());
    ads.fRelAngleA  = trackAB.Angle(trackA.Vect());
    ads.fRelAngleB  = trackAB.Angle(trackB.Vect());

    //Boost in CMS: ***********************************************
    TLorentzVector Vec_beam_target = Vec_pp35_beam+Vec_pp35_target;
    double CMS_Beta     = Vec_beam_target.Beta();
    double CMS_Beta_x   = CMS_Beta*sin(Vec_beam_target.Theta())*cos(Vec_beam_target.Phi());     // x component of BetaVector
    double CMS_Beta_y   = CMS_Beta*sin(Vec_beam_target.Theta())*sin(Vec_beam_target.Phi());     // y component of BetaVector
    double CMS_Beta_z   = CMS_Beta*cos(Vec_beam_target.Theta());                  // z component of BetaVector

    TLorentzVector trackAB_CMS = trackAB;

    trackAB_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.fY_cms               = trackAB_CMS.Rapidity();
    ads.fP_cms               = trackAB_CMS.P();
    ads.fCosTheta_cms        = cos(trackAB_CMS.Theta());

    TLorentzVector trackA_CMS = trackA;
    trackA_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.fMomA_cms            = trackA_CMS.P();

    TLorentzVector trackB_CMS = trackB;
    trackB_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.fMomB_cms            = trackB_CMS.P();

    ads.fXf = fabs(trackAB_CMS.Pz()/vec_beam_cms.Pz());
    ads.fXfi = ads.fXf;

    //*************************************************************

    ads.fMt                = trackAB.Mt();               // Transverse mass
    ads.fTheta            = trackAB.Theta() * R2D;
    ads.fCosTheta        = cos(trackAB.Theta());
    ads.fPhi            = trackAB.Phi() * R2D;

    if (flag_nosecvertcuts == 0)
    {
    if (analysisType == KT::Sim)     //for simulated data
    {
//                 if( !(
//                      fMinTrackDist < par_Mtd/*kMTD*/
//                      && VerDistA > par_VertDistA/*kVDAB*/
//                      && VerDistB > par_VertDistB/*kVDAB*/
//                      && VerDistX > par_VertDistX/*kVDX*/
//                      &&
//                      GeantxVertexA == GeantxVertexB      // the simulated vertex of particleA and particleB has to be the same
//                      && GeantyVertexA == GeantyVertexB
//                      && GeantzVertexA == GeantzVertexB
//                     fRealLambda
//                           )
//                   )
//                 {
//                     continue;
//                 }
    }
    else       //for experimental data
    {
//                 if ( !(
//                     fMinTrackDist < par_Mtd/*kMTD*/
//                      && VerDistA > par_VertDistA/*kVDAB*/
//                      && VerDistB > par_VertDistB/*kVDAB*/
//                      && VerDistX > par_VertDistX/*kVDX*/
//                          )
//                       )
//                 {
//                     continue;
//                 }
    }
    }

//     if (analysisType == KT::Sim)     //for simulated data
//     {
//         if( !(GeantxVertexA == GeantxVertexB      // the simulated vertex of particleA and particleB has to be the same
//             && GeantyVertexA == GeantyVertexB
//             && GeantzVertexA == GeantzVertexB
//                 )
//             )
//         {
//             return ERR_SIM_VERTEX_MISSMATCH;
//         }
//     }

    ++ads.fEventLCounter;

    A_PID = KT::pip;
    trackA = *(HParticleCandSim*)pcand->getObject(trackA_num);
    trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));
    KTifini::CalcSegVector(trackA.getZ(), trackA.getR(), trackA.getPhi(), trackA.getTheta(), baseA, dirA);
    ads.fMomA = trackA.getMomentum();

    if (flag_elosscorr)
    {
        // with corr
        momentum_A_corr = eLossCorr->getCorrMom(A_PID, ads.fMomA, trackA.getTheta());
    }
    else
    {
        // no corr
        momentum_A_corr = ads.fMomA;
    }

    trackA.setMomentum(momentum_A_corr);

    trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));

    TLorentzVector trackAB_miss = trackA + trackB;

    ads.fM_miss = trackAB_miss.M();
    ads.fPVA_miss = calcAngleVar(v1, trackAB_miss);

    ads.ret = (LAMBDA_MASS - ads.fM);
    return ads;
}

AnaDataSet pp45_Lambda::singleFwDetPairAnalysis(Int_t event_num, HCategory * pcand, HCategory * vcand, int trackA_num, int trackB_num, bool quick_run)
{
    return AnaDataSet();
}

void pp45_Lambda::configureTree(TTree * tree)
{
    tree->Branch("fE",              &g_ads.fE,              "fE/F");
    tree->Branch("fM",              &g_ads.fM,              "fM/F");
    tree->Branch("fM_miss",         &g_ads.fM_miss,         "fM_miss/F");
    tree->Branch("fMt",             &g_ads.fMt,             "fMt/F");
    tree->Branch("fTheta",          &g_ads.fTheta,          "fTheta/F");
    tree->Branch("fCosTheta",       &g_ads.fCosTheta,       "fCosTheta/F");
    tree->Branch("fCosTheta_cms",   &g_ads.fCosTheta_cms,   "fCosTheta_cms/F");
    tree->Branch("fPhi",            &g_ads.fPhi,            "fPhi/F");
    tree->Branch("fP",              &g_ads.fP,              "fP/F");
    tree->Branch("fP_cms",          &g_ads.fP_cms,          "fP_cms/F");
    tree->Branch("fPt",             &g_ads.fPt,             "fPt/F");
    tree->Branch("fY",              &g_ads.fY,              "fY/F");
    tree->Branch("fY_cms",          &g_ads.fY_cms,          "fY_cms/F");
    tree->Branch("fPolX",           &g_ads.fPolX,           "fPolX/F");
    tree->Branch("fPolY",           &g_ads.fPolY,           "fPolY/F");
    tree->Branch("fPolZ",           &g_ads.fPolZ,           "fPolZ/F");
    tree->Branch("fZetaX",          &g_ads.fZetaX,          "fZetaX/F");
    tree->Branch("fZetaY",          &g_ads.fZetaY,          "fZetaY/F");
    tree->Branch("fZetaZ",          &g_ads.fZetaZ,          "fZetaZ/F");
    tree->Branch("fXf",             &g_ads.fXf,             "fXf/F");
    tree->Branch("fPolXi",          &g_ads.fPolXi,          "fPolXi/F");
    tree->Branch("fPolYi",          &g_ads.fPolYi,          "fPolYi/F");
    tree->Branch("fPolZi",          &g_ads.fPolZi,          "fPolZi/F");
    tree->Branch("fZetaXi",         &g_ads.fZetaXi,         "fZetaXi/F");
    tree->Branch("fZetaYi",         &g_ads.fZetaYi,         "fZetaYi/F");
    tree->Branch("fZetaZi",         &g_ads.fZetaZi,         "fZetaZi/F");
    tree->Branch("fXfi",            &g_ads.fXfi,            "fXfi/F");

    tree->Branch("fMinTrackDist",   &g_ads.fMTD,            "fMinTrackDist/F" );
    tree->Branch("fVertDistX",      &g_ads.fVertDistX,      "fVertDistX/F");
    tree->Branch("fPVA",            &g_ads.fPVA,            "fPVA/F");
    tree->Branch("fPVA_miss",       &g_ads.fPVA_miss,       "fPVA_miss/F");

    tree->Branch("fPrimVertexX",    &g_ads.fPrimVertexX,    "fPrimVertexX/F");
    tree->Branch("fPrimVertexY",    &g_ads.fPrimVertexY,    "fPrimVertexY/F");
    tree->Branch("fPrimVertexZ",    &g_ads.fPrimVertexZ,    "fPrimVertexZ/F");
    tree->Branch("fPrimVertexR",    &g_ads.fPrimVertexR,    "fPrimVertexR/F");

    tree->Branch("fFitVertexX",     &g_ads.fFitVertexX,     "fFitVertexX/F");
    tree->Branch("fFitVertexY",     &g_ads.fFitVertexY,     "fFitVertexY/F");
    tree->Branch("fFitVertexZ",     &g_ads.fFitVertexZ,     "fFitVertexZ/F");

    tree->Branch("fDecayVertexX",   &g_ads.fDecayVertexX,   "fDecayVertexX/F");
    tree->Branch("fDecayVertexY",   &g_ads.fDecayVertexY,   "fDecayVertexY/F");
    tree->Branch("fDecayVertexZ",   &g_ads.fDecayVertexZ,   "fDecayVertexZ/F");
    tree->Branch("fDecayVertexR",   &g_ads.fDecayVertexR,   "fDecayVertexR/F");

    tree->Branch("fEventVertexX",   &g_ads.fEventVertexX,   "fEventVertexX/F");
    tree->Branch("fEventVertexY",   &g_ads.fEventVertexY,   "fEventVertexY/F");
    tree->Branch("fEventVertexZ",   &g_ads.fEventVertexZ,   "fEventVertexZ/F");

    tree->Branch("fDVres",          &g_ads.fDVres,          "fDVres/F");

    tree->Branch("fMassA",          &g_ads.fMassA,          "fMassA/F");
    tree->Branch("fMassB",          &g_ads.fMassB,          "fMassB/F" );
    tree->Branch("fMomA",           &g_ads.fMomA,           "fMomA/F");
    tree->Branch("fMomB",           &g_ads.fMomB,           "fMomB/F" );
//    tree->Branch("fMomAx",        &g_ads.fMomAx,          "fMomAx/F");
//    tree->Branch("fMomAy",        &g_ads.fMomAy,          "fMomAy/F");
//    tree->Branch("fMomAz",        &g_ads.fMomAz,          "fMomAz/F");
//    tree->Branch("fMomBx",        &g_ads.fMomBx,          "fMomBx/F" );
//    tree->Branch("fMomBy",        &g_ads.fMomBy,          "fMomBy/F" );
//    tree->Branch("fMomBz",        &g_ads.fMomBz,          "fMomBz/F" );
    tree->Branch("fMomA_cms",       &g_ads.fMomA_cms,       "fMomA_cms/F");
    tree->Branch("fMomB_cms",       &g_ads.fMomB_cms,       "fMomB_cms/F");
    tree->Branch("fAngleAB",        &g_ads.fAngleAB,        "fAngleAB/F");
    tree->Branch("fRelAngleA",      &g_ads.fRelAngleA,      "fRelAngleA/F");
    tree->Branch("fRelAngleB",      &g_ads.fRelAngleB,      "fRelAngleB/F");

    tree->Branch("fMDCdEdxA",       &g_ads.fMDCdEdxA,       "fMDCdEdxA/F");
    tree->Branch("fMDCdEdxB",       &g_ads.fMDCdEdxB,       "fMDCdEdxB/F");

    tree->Branch("fChiA",           &g_ads.fChiA,           "fChiA/F");
    tree->Branch("fChiB",           &g_ads.fChiB,           "fChiB/F");
    tree->Branch("fMetaMatchQA",    &g_ads.fMetaMatchQA,    "fMetaMatchQA/F");
    tree->Branch("fMetaMatchQB",    &g_ads.fMetaMatchQB,    "fMetaMatchQB/F");
    tree->Branch("fVertDistA",      &g_ads.fVertDistA,      "fVertDistA/F");
    tree->Branch("fVertDistB",      &g_ads.fVertDistB,      "fVertDistB/F");
    tree->Branch("fHadesTracksNum", &g_ads.fHadesTracksNum, "fHadesTracksNum/I");
    tree->Branch("fFwDetTracksNum", &g_ads.fFwDetTracksNum, "fFwDetTracksNum/I");
    tree->Branch("fIsFwDetData",    &g_ads.fIsFwDetData,    "fIsFwDetData/I");

    tree->Branch("fMultA",          &g_ads.fMultA,          "fMultA/I");
    tree->Branch("fMultB",          &g_ads.fMultB,          "fMultB/I");

    tree->Branch("fAccA",           &g_ads.fAccA,           "fAccA/I");
    tree->Branch("fAccB",           &g_ads.fAccB,           "fAccB/I");

    tree->Branch("fBestTry",        &g_ads.fBestTry,        "fBestTry/I");
    tree->Branch("fEventLCounter",  &g_ads.fEventLCounter,  "fEventLCounter/I");

    if (analysisType == KT::Sim)
    {
        tree->Branch("fRealLambda", &g_ads.fRealLambda,     "fRealLambda/I");
        tree->Branch("fPrimLambda", &g_ads.fPrimLambda,     "fPrimLambda/I");
        tree->Branch("fGeantInfoNum",   &g_ads.fGeantInfoNum,   "fGeantInfoNum/I");
        tree->Branch("fGeantWeight",    &g_ads.fGeantWeight,"fGeantWeight/F");
    }

    tree->Branch("fWallT",          &g_ads.fWallT,          "fWallT/F");
    tree->Branch("fWallX",          &g_ads.fWallX,          "fWallX/F");
    tree->Branch("fWallY",          &g_ads.fWallY,          "fWallY/F");
    tree->Branch("fWallR",          &g_ads.fWallR,          "fWallR/F");
    tree->Branch("fWallCharge",     &g_ads.fWallCharge,     "fWallCharge/F");
    tree->Branch("fWallP",          &g_ads.fWallP,          "fWallP/F");
    tree->Branch("fWallPx",         &g_ads.fWallPx,         "fWallPx/F");
    tree->Branch("fWallPy",         &g_ads.fWallPy,         "fWallPy/F");
    tree->Branch("fWallPz",         &g_ads.fWallPz,         "fWallPz/F");
    tree->Branch("fWallBeta",       &g_ads.fWallBeta,       "fWallBeta/F");
    tree->Branch("fWallGamma",      &g_ads.fWallGamma,      "fWallGamma/F");
    tree->Branch("fWallClusterSize",&g_ads.fWallClusterSize,"fWallClusterSize/I");
    tree->Branch("fWallClustersNum",&g_ads.fWallClustersNum,"fWallClustersNum/I");

    tree->Branch("fGeaP",           &g_ads.fGeaP,           "fGeaP/F");
    tree->Branch("fGeaPx",          &g_ads.fGeaPx,          "fGeaPx/F");
    tree->Branch("fGeaPy",          &g_ads.fGeaPy,          "fGeaPy/F");
    tree->Branch("fGeaPz",          &g_ads.fGeaPz,          "fGeaPz/F");
    tree->Branch("fGeaTheta",       &g_ads.fGeaTheta,       "fGeaTheta/F");
    tree->Branch("fGeaPhi",         &g_ads.fGeaPhi,         "fGeaPhi/F");
    tree->Branch("fGeaPolX",        &g_ads.fGeaPolX,        "fGeaPolX/F");
    tree->Branch("fGeaPolY",        &g_ads.fGeaPolY,        "fGeaPolY/F");
    tree->Branch("fGeaPolZ",        &g_ads.fGeaPolZ,        "fGeaPolZ/F");
    tree->Branch("fGeaZetaX",       &g_ads.fGeaZetaX,       "fGeaZetaX/F");
    tree->Branch("fGeaZetaY",       &g_ads.fGeaZetaY,       "fGeaZetaY/F");
    tree->Branch("fGeaZetaZ",       &g_ads.fGeaZetaZ,       "fGeaZetaZ/F");
    tree->Branch("fGeaXf",          &g_ads.fGeaXf,          "fGeaXf/F");
    tree->Branch("fGeaAngleAB",     &g_ads.fGeaAngleAB,     "fGeaAngleAB/F");

    tree->Branch("fPVtype",         &g_ads.fPVtype,         "fPVtype/I");

#ifdef SHOWREJECTED
    tree->Branch("fRejected",       &g_ads.fRejected,       "fRejected/I");
#endif /*SHOWREJECTED*/

    tree->Branch("fSortOrderMtd",   &g_ads.fSortOrderMtd,   "fSortOrderMtd/I");
}

void pp45_Lambda::configureGraphicalCuts(KTrackInspector & trackInsp)
{
    const TString jsieben_pNb_cuts = "/scratch/e12f/knucl/jsieben/pNb/Cuts/";
    const TString aschmah_pp35_cuts = "/scratch/e12f/schmah/GraphicalCuts/pp35/";
    const TString jchen_pp35_cuts = "/home/gu27buz/hadesdst/pp35/";
    const TString jchen_pp35_cuts_sim = "/scratch/e12l/knucl/hades/jchen/pp35/GraphicalCuts/Sim/";
    const TString jchen_pp35_cuts_sim2 = "/scratch/e12f/knucl/rlalik/pp35/LambdaAnalysis/graph_cuts/";//"/scratch/e12l/knucl/hades/jchen/pp35/new_backup/GraphicalCuts/Sim/";

//    setChargePID(kTRUE);

    const TString rlalik_cuts = "/scratch/e12m/knucl/rlalik/pp35/LambdaAnalysis/Exp/";
    if (analysisType == KT::Sim)
    {
    // protons
//         trackInsp.registerCut(KT::MDC, KT::cut_p, jchen_pp35_cuts_sim2 + "Modified_PID_Cuts_Poly5_ChiiV5_Meth2.root", "Mdc_dEdx_P_cut_mod_ChiiV1_Sim_mod");
    // pions-
//         trackInsp.registerCut(KT::MDC, KT::cut_pim, jchen_pp35_cuts_sim2 + "Modified_PID_Cuts_Poly5_ChiiV5_Meth2.root", "Mdc_dEdx_PiP_cut_PID_mod_ChiiV2_Sim_mod_PiM");

//         trackInsp.registerdEdxPlot(KT::TOF);
//         trackInsp.registerdEdxPlot(KT::TOFINO);
    }
    else if (analysisType == KT::Exp)
    {
    // protons
//         trackInsp.registerCut(KT::MDC, KT::cut_p, jchen_pp35_cuts + "Mdc_dEdx_P_cut_mod_ChiiV1.root", "Mdc_dEdx_P_cut_mod_ChiiV1");
    // pions-
//         trackInsp.registerCut(KT::MDC, KT::cut_pim, jchen_pp35_cuts + "Mdc_dEdx_PiP_cut_PID_mod_ChiiV2.root", "Mdc_dEdx_PiP_cut_PID_mod_ChiiV2", kFALSE);    // new cat is already for pim
    }

    trackInsp.configureMetaSystem(KT::cut_p, KT::MDC);
    trackInsp.configureMetaSystem(KT::cut_pim, KT::MDC);
}

void pp45_Lambda::initAnalysis(KT::Experiment exp, KT::AnalysisType analysisType)
{
    KAbstractAnalysis::initAnalysis(exp, analysisType);

    refBeamVector = getTargetGeomVector();
    beamCal = new KBeamCalibration(getTargetGeometry());
//     beamCal->initBeamCorrArray(analysisType);
}

void pp45_Lambda::finalizeAnalysis()
{
    delete beamCal;
    KAbstractAnalysis::finalizeAnalysis();
}

int pp45_Lambda::Configure(int argc, char ** argv)
{
    int c;
    optind = 1;
    
    while (1) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            { "nosecvercuts",        no_argument,        &flag_nosecvertcuts,    1 },
            { "elosscorr",            no_argument,        &flag_elosscorr,        1 },
            { "no-sigmas",            no_argument,        &flag_nosigmas,            1 },
            { "use-event-vertex",    no_argument,        &flag_useeventvertex,    1 },
            { "use-wall",            no_argument,        &flag_usewall,    1 },
            /* These options don't set a flag.
            We distingu*ish them by their indices. */
//            {"events",     no_argument,       0, 'e'},
            { "Mtd",                required_argument,    0, 'm'},
            { "VertDistX",            required_argument,    0, 'x'},
            { "VertDistA",            required_argument,    0, 'p'},
            { "VertDistB",            required_argument,    0, 'q'},
            { "help",                no_argument,        0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "h", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
            case 0:
                if (long_options[option_index].flag != 0)
                    break;
                break;
            case 'm':
                par_Mtd = atol(optarg);
                break;
            case 'x':
                par_VertDistX = atol(optarg);
                break;
            case 'p':
                par_VertDistA = atol(optarg);
                break;
            case 'q':
                par_VertDistB = atol(optarg);
                break;
            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                Usage();
                break;

            default:
                abort();
                break;
        }
    }

    return optind;
}

void pp45_Lambda::Usage() const
{
    std::cout <<
    "Analysis options: \n" <<
    "      --nosecvercuts\t\t\t - disable secondary vertex cuts\n" <<
    "\n\n";
}
