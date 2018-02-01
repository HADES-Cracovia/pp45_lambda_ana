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

#include "KTools.h"
#include "KTrackInspector.h"
#include "KBeamCalibration.h"

#include "pp35_Lambda.h"

#include "HParticleCandSim.h"

#ifdef SHOWREJECTED
//int fRejected;
#endif /*SHOWREJECTED*/

static const float E_kin_beam       = 3500.0;
static const float E_kin_target     = 0.0;
static const float E_total_beam     = E_kin_beam + HPhysicsConstants::mass(14);
static const float E_total_target   = E_kin_target + HPhysicsConstants::mass(14);
static const float pz_beam          = sqrt(E_total_beam*E_total_beam-HPhysicsConstants::mass(14)*HPhysicsConstants::mass(14));

static const float LAMBDA_MASS = 1115.9;

static const float D2R = TMath::DegToRad();
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

pp35_Lambda::pp35_Lambda(const TString& analysisName, const TString& treeName) : KAbstractAnalysis(analysisName, treeName),
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

bool pp35_Lambda::analysis(HEvent * fEvent, Int_t event_num, HCategory * pcand, Int_t cand_size)
{
//     printf("%d  %d\n", cand_size, vect_size);
    if ( (cand_size) < 2 )
        return false;

    int A_PID = KT::p;    // proton
    int B_PID = KT::pim;    // pi-

    g_ads.fHadesTracksNum = cand_size;
    g_ads.fFwDetTracksNum = 0;
    g_ads.fIsFwDetData = 0;

    g_ads.fGeantWeight = 0;

    g_ads.fEventVertexX = fEvent->getHeader()->getVertex().getX();
    g_ads.fEventVertexY = fEvent->getHeader()->getVertex().getY();
    g_ads.fEventVertexZ = fEvent->getHeader()->getVertex().getZ();

    g_ads.fEventLCounter = 0;

    g_ads.fMultA = 0;
    g_ads.fMultB = 0;

    g_ads.fIsFwDetData = 0;

    std::vector<AnaDataSet> ads_arr;
    ads_arr.reserve(10000);

    size_t combo_cnt = 0;
    std::vector<float> mtd_list;

    for(int i = 0; i < /*fMultA*/g_ads.fHadesTracksNum; ++i)
    {
        for(int j = 0; j < /*fMultB*/g_ads.fHadesTracksNum; ++j)
        {
            if (!hades_tracks[i].pid[A_PID][KT::Charge])
                break;

            if (i == j)
                continue;

            if (!hades_tracks[j].pid[B_PID][KT::Charge])
                continue;

            AnaDataSet ads_ret = singlePairAnalysis(fEvent, event_num, A_PID, B_PID, pcand, i, j);
            ads_ret.fIsFwDetData = 0;

            if (ads_ret.ret == ERR_MASS_OUT_OF_RANGE)
                continue;

            if (ads_ret.ret == ERR_NO_LAMBDA)
                continue;

            if (ads_ret.ret <= NULL_DATA)
                continue;

            hades_tracks[i].is_used = true;
            hades_tracks[j].is_used = true;

            mtd_list.push_back(ads_ret.fMTD);

            ads_arr.push_back(ads_ret);
            ++combo_cnt;
        }
    }

    std::sort(mtd_list.begin(), mtd_list.end());

    for (size_t i = 0; i < combo_cnt; ++i)
    {
        for (size_t j = 0; j < combo_cnt; ++j)
        {
            if (ads_arr[i].fMTD == mtd_list[j])
                ads_arr[i].fSortOrderMtd = j;
        }
    }

    for (size_t i = 0; i < combo_cnt; ++i)
    {
        g_ads = ads_arr[i];

        track_lambda_cms = ads_arr[i].tr_lambda_cms;
        track_lambda_a = ads_arr[i].tr_lambda_a;
        track_lambda_b = ads_arr[i].tr_lambda_b;
        vertex_primary = ads_arr[i].vx_primary;
        vertex_lambda = ads_arr[i].vx_lambda;
        vector_pol = ads_arr[i].vr_pol;
        vector_zeta = ads_arr[i].vr_zeta;
        vector_poli = ads_arr[i].vr_poli;
        vector_zetai = ads_arr[i].vr_zetai;
        trec_lambda = ads_arr[i].trec_lambda;
        trec_lambdaF = ads_arr[i].trec_lambdaF;

        track_lambda_cms.fill();
        track_lambda_a.fill();
        track_lambda_b.fill();
        vertex_primary.fill();
        vertex_lambda.fill();
        vector_pol.fill();
        vector_zeta.fill();
        vector_poli.fill();
        vector_zetai.fill();
        trec_lambda.fill();
        trec_lambdaF.fill();

        getTree()->Fill();
    }

    return true;
}

AnaDataSet pp35_Lambda::singlePairAnalysis(HEvent * fEvent, int /*event_num*/, UInt_t pid_a, UInt_t pid_b, HCategory * pcand, int trackA_num, int trackB_num, bool quick_run)
{
    HGeomVector beamVector; // FIXME
//     if (analysisType == KT::Exp)
//     {
//         beamVector = beamCal->calculateBeamOffset(event->getRunId());
//     } else {
        beamVector = refBeamVector;
//     }

    TVector3 p_beam_vec(beamVector.X(), beamVector.Y(), 3000.0);
//     p_beam_vec.SetMag(pz_beam); FIXME

    const TLorentzVector Vec_pp35_beam      = TLorentzVector(p_beam_vec.X(), p_beam_vec.Y(), p_beam_vec.Z(), E_total_beam);
    const TLorentzVector Vec_pp35_target    = TLorentzVector(0.0, 0.0, 0.0, E_total_target);
    const TLorentzVector Vec_pp35_sum       = Vec_pp35_beam + Vec_pp35_target;
    const float cmrap                       = Vec_pp35_sum.Rapidity();

    TLorentzVector vec_beam_cms             = Vec_pp35_beam;
    vec_beam_cms.Boost(-Vec_pp35_sum.BoostVector());


    AnaDataSet ads = g_ads;
    ads.init();
//     ads.fGeantWeight = pcand->getGeantGenweight(); FIXME

    HGeomVector dirMother, PrimVertexMother;

    TObject * o_a = pcand->getObject(trackA_num);
    TObject * o_b = pcand->getObject(trackB_num);

    HParticleCand trackA = *(HPidTrackCand*)o_a;
    HParticleCand trackB = *(HPidTrackCand*)o_b;

    trackA.calc4vectorProperties(HPhysicsConstants::mass(pid_a));
    trackB.calc4vectorProperties(HPhysicsConstants::mass(pid_b));

    float momentum_A_corr = 0;
    float momentum_B_corr = 0;

    Float_t fMomA = trackA.getMomentum();
    Float_t fMomB = trackB.getMomentum();

    if (flag_elosscorr and analysisType == KT::Exp)
    {
        // with corr
        momentum_A_corr = eLossCorr->getCorrMom(pid_a, fMomA, trackA.getTheta());
        momentum_B_corr = eLossCorr->getCorrMom(pid_b, fMomB, trackB.getTheta());
    }
    else
    {
        // no corr
        momentum_A_corr = fMomA;
        momentum_B_corr = fMomB;
        momentum_A_corr = eLossCorr->getCorrMom(pid_a, fMomA, trackA.getTheta());
        momentum_B_corr = eLossCorr->getCorrMom(pid_b, fMomB, trackB.getTheta());

    }

    trackA.setMomentum(momentum_A_corr);
    trackB.setMomentum(momentum_B_corr);

    trackA.calc4vectorProperties(HPhysicsConstants::mass(pid_a));
    trackB.calc4vectorProperties(HPhysicsConstants::mass(pid_b));

    ads.fMomAx = trackA.Px();
    ads.fMomAy = trackA.Py();
    ads.fMomAz = trackA.Pz();

    ads.fMomBx = trackB.Px();
    ads.fMomBy = trackB.Py();
    ads.fMomBz = trackB.Pz();

    ads.tr_lambda_a = trackA;
    ads.tr_lambda_b = trackB;
    ads.trec_lambda.reconstruct(trackA, trackB);

//    printf("L: px=%f  py=%f  pz=%f  e=%f\n", lambda.Px(), lambda.Py(), lambda. Pz(), ads.fE);
//    printf("beta=%f p_p=%f p_l=%f\n", beta_vec.Mag(), trackA.P(), l_lambda.P());
//    PR(ads.fPol);

    // I do not need so many data!
    if (ads.trec_lambda.M() > 1200)
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
        ads.ret = (LAMBDA_MASS - ads.trec_lambda.M());
        return ads;// * fMinTrackDist;
    }

    ads.fMTD = ads.trec_lambda.getMTD();        // minimum distance between the two tracks
//     if (quick_run) printf(" mtd: %f\n", fMinTrackDist);

    float GeantxVertexA    = 0;
    float GeantyVertexA    = 0;
    float GeantzVertexA    = 0;
    float GeantxVertexB    = 0;
    float GeantyVertexB    = 0;
    float GeantzVertexB    = 0;

    // extra checks for the simulation analysis
    HPidTrackCandSim * tcs_a = dynamic_cast<HPidTrackCandSim*>(o_a);
    HPidTrackCandSim * tcs_b = dynamic_cast<HPidTrackCandSim*>(o_b);

    if (analysisType == KT::Sim && tcs_a && tcs_b)
    {
        HParticleCandSim trackA = *tcs_a;
        HParticleCandSim trackB = *tcs_b;

        TLorentzVector geaA; geaA.SetXYZM(trackA.getGeantxMom(), trackA.getGeantyMom(), trackA.getGeantzMom(), HPhysicsConstants::mass(pid_a));
        TLorentzVector geaB; geaB.SetXYZM(trackB.getGeantxMom(), trackB.getGeantyMom(), trackB.getGeantzMom(), HPhysicsConstants::mass(pid_b));
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

        int GeantPIDA = trackA.getGeantPID();
        int GeantPIDB = trackB.getGeantPID();

        int GeantPIDAparent = trackA.getGeantParentPID();
        int GeantPIDBparent = trackB.getGeantParentPID();

        int GeantPIDAGparent = trackA.getGeantGrandParentPID();
        int GeantPIDBGparent = trackB.getGeantGrandParentPID();

        // the simulated vertex of particleA and particleB has to be the same
        ads.fRealLambda = ( (GeantPIDA == 14 and GeantPIDAparent == 18 and GeantPIDB == 9 and GeantPIDBparent == 18) or (GeantPIDA == 14 and GeantPIDAparent == 18 and GeantPIDB == 6 and GeantPIDBparent == 9 and GeantPIDBGparent == 18));

//         if (ads.fRealLambda)
//         {
//             printf("*************** TRACK A ***************\n");
//             trackA.print(1<<4 | 1<<2);
//             printf("*************** TRACK B ***************\n");
//             trackB.print(1<<4 | 1<<2);
//         }

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

    dirMother.setXYZ(ads.trec_lambda.X(), ads.trec_lambda.Y(), ads.trec_lambda.Z());    // direction vector of the mother particle
    ads.vx_lambda = ads.trec_lambda.getDecayVertex();

    if (analysisType == KT::Sim)
    {
        ads.fDVres = TMath::Sqrt(
                TMath::Power(ads.vx_lambda.X() - GeantxVertexA, 2) +
                TMath::Power(ads.vx_lambda.Y() - GeantyVertexA, 2) +
                TMath::Power(ads.vx_lambda.Z() - GeantzVertexA, 2)
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
    PrimVertexMother.setXYZ(ads.fEventVertexX, ads.fEventVertexY, ads.fEventVertexZ);

    ads.vx_primary = PrimVertexMother;
// printf("vertex z: prim=%f   sec=%f\n", ads.fPrimVertexZ, ads.fDecayVertexZ);
    if ((ads.vx_lambda.Z() - ads.vx_primary.Z()) < 0)
    {
#ifdef SHOWREJECTED
        ads.fRejected = ERR_VERTEX_Z_MISSMATCH;
#else
//        ads.ret = ERR_VERTEX_Z_MISSMATCH;
//        return ads;
#endif /*SHOWREJECTED*/
    }

    ads.fFitVertexX = (float)vvectors.second.getX();
    ads.fFitVertexY = (float)vvectors.second.getY();
    ads.fFitVertexZ = (float)vvectors.second.getZ();

    HGeomVector v1 = ads.vx_lambda - PrimVertexMother;
//    HGeomVector v2 = DecayVertex - trackAB;
    //vvectors.second;

    ads.fVertDistX = (v1).length();

    TVector3 _v1; _v1.SetXYZ(v1.X(), v1.Y(), v1.Z());
    TVector3 _v2; _v2.SetXYZ(ads.trec_lambda.Px(), ads.trec_lambda.Py(), ads.trec_lambda.Pz());

    ads.fPVA = calcAngleVar(v1, ads.trec_lambda);

    ads.fVertDistA = ads.trec_lambda.getMTDa();
    ads.fVertDistB = ads.trec_lambda.getMTDb();

    // HERE we assign PV-SV vactor for Lambda

//    TVector3 vLpsv = _v1;
//    vLpsv.SetMag(ads.trec_lambda.P());

//    TParticleTrack Lpsv;
//    Lpsv.setMomentum(vLpsv.Mag());
//    Lpsv.setTheta(vLpsv.Theta() * R2D);
//    Lpsv.setPhi(vLpsv.Phi() * R2D);
//    Lpsv.calc4vectorProperties(HPhysicsConstants::mass(18));

    TLorentzVector lambdaAB = ads.trec_lambda;
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
        ads.vr_pol.setX( p_proton.Dot(n_x) / (p_proton.Mag() * n_x.Mag()) );
        ads.vr_pol.setY( p_proton.Dot(n_y) / (p_proton.Mag() * n_y.Mag()) );
        ads.vr_pol.setZ( p_proton.Dot(n_z) / (p_proton.Mag() * n_z.Mag()) );

        ads.vr_zeta.setX( TMath::ACos(ads.vr_pol.getX()) * R2D );
        ads.vr_zeta.setY( TMath::ACos(ads.vr_pol.getY()) * R2D );
        ads.vr_zeta.setZ( TMath::ACos(ads.vr_pol.getZ()) * R2D );
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
        ads.vr_poli.setX( p_proton.Dot(n_x) / (p_proton.Mag() * n_x.Mag()) );
        ads.vr_poli.setY( p_proton.Dot(n_y) / (p_proton.Mag() * n_y.Mag()) );
        ads.vr_poli.setZ( p_proton.Dot(n_z) / (p_proton.Mag() * n_z.Mag()) );

        ads.vr_zetai.setX( TMath::ACos(ads.vr_poli.getX()) * R2D );
        ads.vr_zetai.setY( TMath::ACos(ads.vr_poli.getY()) * R2D );
        ads.vr_zetai.setZ( TMath::ACos(ads.vr_poli.getZ()) * R2D );
    }

//     double dist2 = pow(PrimVertexMother.getX() - beamVector.getX(), 2.0) +
//         pow(PrimVertexMother.getY() - beamVector.getY(), 2.0);

//     if ( !(dist2 < 100.0 and PrimVertexMother.getZ() < 0.0 and PrimVertexMother.getZ() > -90.0) )
//         return ERR_NOT_IN_VERTEX;
    if ( !(ads.vx_primary.getR() < 10.0 and ads.vx_primary.getZ() < 0.0 and ads.vx_primary.getZ() > -90.0) )
    {
#ifdef SHOWREJECTED
        ads.fRejected = ERR_NOT_IN_VERTEX;
#else
        ads.ret = ERR_NOT_IN_VERTEX;
        return ads;
#endif /*SHOWREJECTED*/
    }

    ads.fAngleAB    = trackA.Angle(trackB.Vect());
    ads.fRelAngleA  = ads.trec_lambda.Angle(trackA.Vect());
    ads.fRelAngleB  = ads.trec_lambda.Angle(trackB.Vect());

    //Boost in CMS: ***********************************************
    TLorentzVector Vec_beam_target = Vec_pp35_beam+Vec_pp35_target;
    double CMS_Beta     = Vec_beam_target.Beta();
    double CMS_Beta_x   = CMS_Beta*sin(Vec_beam_target.Theta())*cos(Vec_beam_target.Phi());     // x component of BetaVector
    double CMS_Beta_y   = CMS_Beta*sin(Vec_beam_target.Theta())*sin(Vec_beam_target.Phi());     // y component of BetaVector
    double CMS_Beta_z   = CMS_Beta*cos(Vec_beam_target.Theta());                  // z component of BetaVector

    TLorentzVector trackAB_CMS = ads.trec_lambda;

    trackAB_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.tr_lambda_cms = trackAB_CMS;

    TLorentzVector trackA_CMS = trackA;
    trackA_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.fMomA_cms            = trackA_CMS.P();

    TLorentzVector trackB_CMS = trackB;
    trackB_CMS.Boost(-CMS_Beta_x, -CMS_Beta_y, -CMS_Beta_z);
    ads.fMomB_cms            = trackB_CMS.P();

    ads.fXf = fabs(trackAB_CMS.Pz()/vec_beam_cms.Pz());
    ads.fXfi = ads.fXf;

    //*************************************************************

    ads.fMt                = ads.trec_lambda.Mt();               // Transverse mass

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

//     A_PID = KT::pip;
//     trackA = *(HParticleCandSim*)pcand->getObject(trackA_num);
//     trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));
//     KTifini::CalcSegVector(trackA.getZ(), trackA.getR(), trackA.getPhi(), trackA.getTheta(), baseA, dirA);
//     ads.fMomA = trackA.getMomentum();
// 
//     if (flag_elosscorr)
//     {
//         // with corr
//         momentum_A_corr = eLossCorr->getCorrMom(A_PID, ads.fMomA, trackA.getTheta());
//     }
//     else
//     {
//         // no corr
//         momentum_A_corr = ads.fMomA;
//     }
// 
//     trackA.setMomentum(momentum_A_corr);
// 
//     trackA.calc4vectorProperties(HPhysicsConstants::mass(A_PID));
// 
//     TLorentzVector trackAB_miss = trackA + trackB;
// 
//     ads.fM_miss = trackAB_miss.M();
//     ads.fPVA_miss = calcAngleVar(v1, trackAB_miss);
// 
//     ads.ret = (LAMBDA_MASS - ads.tr_lambda.M());
    return ads;
}

void pp35_Lambda::configureTree(TTree * tree)
{
    trec_lambda.setTree(tree, "Lambda_");
    trec_lambdaF.setTree(tree, "LambdaF_");
    track_lambda_cms.setTree(tree, "Lambda_cms_", KTrack::bCosTheta | KTrack::bP | KTrack::bY);
    track_lambda_a.setTree(tree, "partA_", 0x1ff);
    track_lambda_b.setTree(tree, "partB_", 0x1ff);
    vertex_primary.setTree(tree, "PrimaryVertex", KVertex::bXYZ | KVertex::bR);
    vertex_lambda.setTree(tree, "LambdaDecay", KVertex::bXYZ | KVertex::bR);
    vector_pol.setTree(tree, "Pol");
    vector_zeta.setTree(tree, "Zeta");
    vector_poli.setTree(tree, "PolI");
    vector_zetai.setTree(tree, "ZetaI");

    tree->Branch("fM_miss",         &g_ads.fM_miss,         "fM_miss/F");
    tree->Branch("fMt",             &g_ads.fMt,             "fMt/F");
    tree->Branch("fXf",             &g_ads.fXf,             "fXf/F");
    tree->Branch("fXfi",            &g_ads.fXfi,            "fXfi/F");

    tree->Branch("fMinTrackDist",   &g_ads.fMTD,            "fMinTrackDist/F" );
    tree->Branch("fVertDistX",      &g_ads.fVertDistX,      "fVertDistX/F");
    tree->Branch("fPVA",            &g_ads.fPVA,            "fPVA/F");
    tree->Branch("fPVA_miss",       &g_ads.fPVA_miss,       "fPVA_miss/F");

    tree->Branch("fFitVertexX",     &g_ads.fFitVertexX,     "fFitVertexX/F");
    tree->Branch("fFitVertexY",     &g_ads.fFitVertexY,     "fFitVertexY/F");
    tree->Branch("fFitVertexZ",     &g_ads.fFitVertexZ,     "fFitVertexZ/F");

    tree->Branch("fEventVertexX",   &g_ads.fEventVertexX,   "fEventVertexX/F");
    tree->Branch("fEventVertexY",   &g_ads.fEventVertexY,   "fEventVertexY/F");
    tree->Branch("fEventVertexZ",   &g_ads.fEventVertexZ,   "fEventVertexZ/F");

    tree->Branch("fDVres",          &g_ads.fDVres,          "fDVres/F");

    tree->Branch("fMomA_cms",       &g_ads.fMomA_cms,       "fMomA_cms/F");
    tree->Branch("fMomB_cms",       &g_ads.fMomB_cms,       "fMomB_cms/F");
    tree->Branch("fAngleAB",        &g_ads.fAngleAB,        "fAngleAB/F");
    tree->Branch("fRelAngleA",      &g_ads.fRelAngleA,      "fRelAngleA/F");
    tree->Branch("fRelAngleB",      &g_ads.fRelAngleB,      "fRelAngleB/F");

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

void pp35_Lambda::configureGraphicalCuts(KTrackInspector & trackInsp)
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

void pp35_Lambda::initAnalysis(KT::Experiment exp, KT::AnalysisType analysisType)
{
    KAbstractAnalysis::initAnalysis(exp, analysisType);

    refBeamVector = getTargetGeomVector();
    beamCal = new KBeamCalibration(getTargetGeometry());
//     beamCal->initBeamCorrArray(analysisType);
}

void pp35_Lambda::finalizeAnalysis()
{
    delete beamCal;
    KAbstractAnalysis::finalizeAnalysis();
}

int pp35_Lambda::Configure(int argc, char ** argv)
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

void pp35_Lambda::Usage() const
{
    std::cout <<
    "Analysis options: \n" <<
    "      --nosecvercuts\t\t\t - disable secondary vertex cuts\n" <<
    "\n\n";
}
