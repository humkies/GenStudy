#include "NTupleReader.h"

#include <sstream>
#include <iostream>
#include <fstream>

#include "samples.h"
#include "customize.h"

#include "TStopwatch.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"

#include "tdrstyle.h"
#include "TPrincipal.h"
#include "run/pca_TTbar.C"

const bool debug = false;

const double mW_ = 80.385, mTop_ = 173.5;

std::vector<std::string> keyStringCachedVec;
std::vector<double> scaleMCCachedVec;
std::vector<int> colorCachedVec;

double cnt_passLeptVeto_WeightedScaledMC = 0, cnt_passLeptVeto_WeightedErrorScaledMC = 0;
double cnt_passnJets_WeightedScaledMC = 0, cnt_passnJets_WeightedErrorScaledMC = 0;
double cnt_passdPhis_WeightedScaledMC = 0, cnt_passdPhis_WeightedErrorScaledMC = 0;
double cnt_passBJets_WeightedScaledMC = 0, cnt_passBJets_WeightedErrorScaledMC = 0;
double cnt_passMET_WeightedScaledMC = 0, cnt_passMET_WeightedErrorScaledMC = 0;
double cnt_passTagger_WeightedScaledMC = 0, cnt_passTagger_WeightedErrorScaledMC = 0;
double cnt_passBaseline_WeightedScaledMC = 0, cnt_passBaseline_WeightedErrorScaledMC = 0;
double cnt_passMore_WeightedScaledMC = 0, cnt_passMore_WeightedErrorScaledMC = 0;

// onelepton includes: W->e nu, W->mu nu, W->tau->e nu, W->tau->mu nu;
// dileptons includes combinations of both W's to e/mu or W->tau->e/mu
// diWtauhads includes combinations of both W's to W->tau->had.
// dileptons_inc_Wtauhad includes combination of one W to e/mu or W->tau->e/mu, the other W to W->tau-had
// allHad includes both W's decay hadronically
double cnt_onelepton_WeightedScaledMC = 0, cnt_onelepton_WeightedErrorScaledMC = 0;
double cnt_Wtauhad_WeightedScaledMC = 0, cnt_Wtauhad_WeightedErrorScaledMC = 0;
double cnt_dileptons_WeightedScaledMC = 0, cnt_dileptons_WeightedErrorScaledMC = 0;
double cnt_diWtauhads_WeightedScaledMC = 0, cnt_diWtauhads_WeightedErrorScaledMC = 0;
double cnt_dileptons_inc_Wtauhad_WeightedScaledMC = 0, cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC = 0;
double cnt_allHad_WeightedScaledMC = 0, cnt_allHad_WeightedErrorScaledMC = 0;
                               
const double nTops_SR_lo[] = { 0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3 };
const double nTops_SR_hi[] = { 0,  0,  0,  1,  1,  1,  2,  2,  2, -1, -1, -1 };

const double nbJets_SR_lo[]  = { 1,  2,  3,  1,  2,  3,  1,  2,  3,  1,  2,  3 };
const double nbJets_SR_hi[]  = { 1,  2, -1,  1,  2, -1,  1,  2, -1,  1,  2, -1 };
/*
// nTopsLE3 --> >=200
const double met_SR_lo[] = { 200 };
const double met_SR_hi[] = {  -1 };

const std::string keyStr_met_SR[] = { "LE200" };
const std::string disStr_met_SR[] = { "MET>=200" };
*/
// nTopsEQ1 --> 200-350; 350-400; 400-450; 450-500; >=500;
const double met_SR_lo[] = { 200,  350,  400,  450,  500};
const double met_SR_hi[] = { 350,  400,  450,  500,   -1};

const std::string keyStr_met_SR[] = {     "200to350",     "350to400",     "400to450",     "450to500",    "LE500" };
const std::string disStr_met_SR[] = { "200<=MET<350", "350<=MET<400", "400<=MET<450", "450<=MET<500", "MET>=500" };

/*
// nTopsEQ2 --> 200-350; 350-400; 400-475; >=475;
const double met_SR_lo[] = { 200,  350,  400,  475 };
const double met_SR_hi[] = { 350,  400,  475,   -1 };

const std::string keyStr_met_SR[] = {     "200to350",     "350to400",     "400to475",        "LE475" };
const std::string disStr_met_SR[] = { "200<=MET<350", "350<=MET<400", "400<=MET<475",     "MET>=475" };
*/
/*
const double met_SR_lo[] = { 200,  350,  375,  400,  425,  450,  475,  525,  550};
const double met_SR_hi[] = { 350,  375,  400,  425,  450,  475,  500,  550,   -1};

const std::string keyStr_met_SR[] = {     "200to350",     "350to375",     "375to400",     "400to425",     "425to450",     "450to475",     "475to500",     "500to525",     "525to550",    "LE550" };
const std::string disStr_met_SR[] = { "200<=MET<350", "350<=MET<375", "375<=MET<400", "400<=MET<425", "425<=MET<450", "450<=MET<475", "475<=MET<500", "500<=MET<525", "525<=MET<550", "MET>=550" };
*/

const std::string keyStr_nTops_SR[] = { "EQ0", "EQ0", "EQ0", "EQ1", "EQ1", "EQ1", "EQ2", "EQ2", "EQ2", "LE3", "LE3", "LE3" };
const std::string disStr_nTops_SR[] = {  "=0",  "=0",  "=0",  "=1",  "=1",  "=1",  "=2",  "=2",  "=2", ">=3", ">=3", ">=3" };

const std::string keyStr_nbJets_SR[] = { "EQ1", "EQ2", "LE3", "EQ1", "EQ2", "LE3", "EQ1", "EQ2", "LE3", "EQ1", "EQ2", "LE3" };
const std::string disStr_nbJets_SR[] = {  "=1",  "=2", ">=3",  "=1",  "=2", ">=3",  "=1",  "=2", ">=3",  "=1",  "=2", ">=3" };

const int nSR = sizeof(nTops_SR_lo)/sizeof(nTops_SR_lo[0]);
const int nSR_met = sizeof(met_SR_lo)/sizeof(met_SR_lo[0]);

bool doIsoTrksVeto = true;
bool doMT2mTcombCuts = true;

std::vector<TH2D*> h2_evtCnt_nbJets_vs_nTopsVec;
TH2D * h2_evtCnt_sumSM_nbJets_vs_nTops = 0;

std::vector<TPrincipal*> pcaVec;
std::vector<TString> keyWordVec;

std::vector<TH1D*> h1_cutFlowVec, h1_cutFlow_auxVec;
std::vector<TH1D*> h1_dR_gen_reco_topVec;
std::vector<TH1D*> h1_dPhi_lept_metVec, h1_dPhi_nu_metVec, h1_dPhi_sum_lept_nu_metVec;
std::vector<TH1D*> h1_dR_b_topVec, h1_dPhi_b_topVec;
std::vector<TH2D*> h2_mt_b_vs_dPhi_b_topVec;
std::vector<TH1D*> h1_dPhi_rJet_metVec, h1_dPhi_rJet_bVec;
std::vector<TH2D*> h2_dPhi_rJet_b_vs_metVec;

std::vector<TH2D*> h2_dalitzVec, h2_gen_m23overm123vsarctanm13overm12Vec, h2_m23overm123vsarctanm13overm12Vec;

std::vector<TH1D*> h1_dPhi_top_metVec, h1_minDphi_tJet_metVec;
std::vector<TH1D*> h1_MT_tJetVec, h1_MT_tDiJetsVec;

std::vector<TH1D*> h1_pt_genWdau1_nrJets_EQ1Vec, h1_eta_genWdau1_nrJets_EQ1Vec;
std::vector<TH1D*> h1_pt_genWdau2_nrJets_EQ1Vec, h1_eta_genWdau2_nrJets_EQ1Vec;
std::vector<TH1D*> h1_mass_rbJet_mtch0_nrJets_EQ2Vec, h1_mass_rbJet_mtch1_nrJets_EQ2Vec;
std::vector<TH1D*> h1_dR_rbJet_mtch0_nrJets_EQ2Vec, h1_dR_rbJet_mtch1_nrJets_EQ2Vec;
std::vector<TH1D*> h1_dPhi_rJet_met_mtch0_nrJets_EQ2Vec, h1_dPhi_rJet_met_mtch1_nrJets_EQ2Vec;

std::vector<TH1D*> h1_Wjets_bTagCatsVec, h1_Topjets_bTagCatsVec;
std::vector<TH1D*> h1_bTagged_Wjets_CSVVec;
std::vector<TH2D*> h2_bTagged_Wjets_minDR_genb_vs_minDR_genWVec, h2_bTagged_Topjets_minDR_genb_vs_minDR_genTopVec;
std::vector<TH2D*> h2_nobTagged_Wjets_minDR_genb_vs_minDR_genWVec, h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTopVec;
std::vector<TH2D*> h2_bTagged_Wjets_mass_vs_minDR_genbVec;
std::vector<TH1D*> h1_bTagged_Wjets_massVec, h1_nobTagged_Wjets_massVec, h1_bTagged_Wjets_minDR_otherbJetsVec;
std::vector<TH1D*> h1_bTagged_Topjets_massVec, h1_nobTagged_Topjets_massVec;
std::vector<TH1D*> h1_bTagged_Wjets_deltaR1b_genDausVec, h1_bTagged_Wjets_deltaR2b_genDausVec, h1_bTagged_Wjets_deltaR12_genDausVec;
std::vector<TH1D*> h1_bTagged_rndmComb_Wjets_ptVec, h1_bTagged_sameTop_Wjets_ptVec;

std::vector<TH1D*> h1_pt_gentbVec, h1_eta_gentbVec, h1_pt_genrbVec, h1_eta_genrbVec;
std::vector<TH1D*> h1_pt_genrb_match0Vec, h1_eta_genrb_match0Vec;
std::vector<TH1D*> h1_csvs_fakebVec;

std::vector<TH1D*> h1_gen1b_massVec, h1_gen2b_massVec;
std::vector<TH1D*> h1_gen1b_deltaRVec, h1_gen2b_deltaRVec, h1_gen12_deltaRVec;
std::vector<TH1D*> h1_gen1MET_deltaPhiVec, h1_gen2MET_deltaPhiVec;

std::vector<TH1D*> h1_dPhi_genTopsVec, h1_dPhi_genb_genTopVec, h1_dPhi_genbsVec;
std::vector<TH1D*> h1_nComb_lept_brJetVec, h1_nComb_had_brJetVec;
std::vector<TH1D*> h1_MT_lept_brJetVec, h1_MT_had_brJetVec;
std::vector<TH1D*> h1_MT_bJetVec;
std::vector<TH1D*> h1_mass_rTopVec;
std::vector<TH1D*> h1_ori_MT_lept_brJetVec, h1_ori_MT_had_brJetVec;
std::vector<TH1D*> h1_MT2_lept_brJetVec, h1_MT2_had_brJetVec;
std::vector<TH1D*> h1_ori_MT2_lept_brJetVec, h1_ori_MT2_had_brJetVec;
std::vector<TH2D*> h2_nComb_lept_vs_had_brJetVec;
std::vector<TH1D*> h1_dPhi_sum_METVec, h1_MT_sumVec;

std::vector<TH1D*> h1_nComb_top_brJetsVec;

std::vector<TH2D*> h2_MT_vs_MT2_lept_brJetVec, h2_MT_vs_MT2_had_brJetVec, h2_ori_MT2_vs_MT2_had_brJetVec, h2_mTcomb_vs_MT2_lept_brJetVec, h2_mTcomb_vs_MT2_had_brJetVec;
std::vector<TH2D*> h2_MT_lept_vs_MT_bJetVec, h2_MT_bJet_vs_mTcomb_had_brJetVec;
std::vector<TH2D*> h2_mTcomb_vs_MT_lept_brJetVec, h2_mTcomb_vs_MT_had_brJetVec;
std::vector<TH2D*> h2_MT_lept_vs_top_brJetVec, h2_MT_had_vs_top_brJetVec, h2_aft_PCA_MT_had_vs_top_brJetVec;
std::vector<TH2D*> h2_MT2_lept_vs_had_brJetVec, h2_MT_lept_vs_had_brJetVec, h2_MT_lept_vs_MT2_had_brJetVec, h2_MT2_lept_vs_MT_had_brJetVec, h2_mTcomb_lept_vs_had_brJetVec;
std::vector<TH2D*> h2_MT_lept_vs_mTcomb_had_brJetVec, h2_mTcomb_lept_vs_MT_had_brJetVec, h2_MT_lept_vs_aft_PCA_had_brJetVec;
std::vector<TH2D*> h2_MT2_lept_vs_had_aft_MT_cuts_brJetVec;
std::vector<TH1D*> h1_MT_same_lept_had_brJetVec, h1_MT2_same_lept_had_brJetVec, h1_mTcomb_same_lept_had_brJetVec, h1_MT2_same_lept_had_aft_MT_cuts_brJetsVec;

std::vector<TH1D*> h1_adj_MT_lept_brJetVec, h1_adj_MT_had_brJetVec;
std::vector<TH1D*> h1_adj_MT2_lept_brJetVec, h1_adj_MT2_had_brJetVec;
std::vector<TH2D*> h2_adj_MT_vs_MT2_lept_brJetVec, h2_adj_MT_vs_MT2_had_brJetVec, h2_adj_mTcomb_vs_MT2_lept_brJetVec, h2_adj_mTcomb_vs_MT2_had_brJetVec;
std::vector<TH2D*> h2_adj_mTcomb_vs_MT_lept_brJetVec, h2_adj_mTcomb_vs_MT_had_brJetVec;
std::vector<TH2D*> h2_adj_MT_lept_vs_top_brJetVec, h2_adj_MT_had_vs_top_brJetVec;
std::vector<TH2D*> h2_adj_MT2_lept_vs_had_brJetVec, h2_adj_MT_lept_vs_had_brJetVec, h2_adj_MT_lept_vs_MT2_had_brJetVec, h2_adj_MT2_lept_vs_MT_had_brJetVec, h2_adj_mTcomb_lept_vs_had_brJetVec;
std::vector<TH2D*> h2_adj_MT_lept_vs_mTcomb_had_brJetVec, h2_adj_mTcomb_lept_vs_MT_had_brJetVec;
std::vector<TH2D*> h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJetVec;
std::vector<TH1D*> h1_adj_MT_same_lept_had_brJetVec, h1_adj_MT2_same_lept_had_brJetVec, h1_adj_mTcomb_same_lept_had_brJetVec, h1_adj_MT2_same_lept_had_aft_MT_cuts_brJetsVec;


std::vector<std::vector<TH1D*> > h1_nJetsVec, h1_metVec, h1_MT2Vec, h1_mTcombVec, h1_HTVec, h1_nJetsRsysVec;
std::vector<std::string> declaredSampleStrVec;
std::vector<std::vector<TH2D*> > h2_met_vs_nJetsVec, h2_MT2_vs_metVec, h2_coarse_bin_MT2_vs_metVec;

// 1st idx: SR of (nbJets, nTops); 2nd idx: SR of met; 3rd idx: sample
std::vector<std::vector<std::vector<double> > > cnt_passSRmet_WeightedScaledMCVec, cnt_passSRmet_WeightedErrorScaledMCVec;
// 1st idx: SR of (nbJets, nTops); 2nd idx: SR of met
std::vector<std::vector<double> > cnt_passSRmet_sumSM_WeightedScaledMCVec, cnt_passSRmet_sumSM_WeightedErrorScaledMCVec;

std::vector<TH1D*> h1_genLeptORprong_PtVec, h1_genLeptORprong_EtaVec, h1_genLeptORprong_PhiVec;
std::vector<TH1D*> h1_genLeptORprong_minDR_tripletsVec, h1_genLeptORprong_minDphi_tripletsVec, h1_genLeptORprong_minDR_RsysVec, h1_genLeptORprong_minDphi_RsysVec;
std::vector<TH1D*> h1_genLeptORprong_dPhi_metVec;

std::vector<TH2D*> h2_genLeptORprong_minDR_Rsys_vs_tripletsVec, h2_genLeptORprong_minDphi_Rsys_vs_tripletsVec;

std::vector<TH1D*> h1_minMTj_lepJetVec, h1_minDphi_met_lepJetVec;
std::vector<TH2D*> h2_minMTj_pickedIdx_vs_minDR_pickedIdxVec, h2_minDphi_met_pickedIdx_vs_minDR_pickedIdxVec;

std::vector<TH1D*> h1_relPt_genLeptORprongOVERjetPt_tripletsVec, h1_relPt_genLeptORprongOVERjetPt_RsysVec;
std::vector<TH2D*> h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPtVec, h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPtVec;

const int nTopCandToPlot = 3;
std::vector<std::vector<TH1D*> > h1_topCand_MVec, h1_topCand_PtVec, h1_topCand_EtaVec, h1_topCand_nJetsVec;
std::vector<std::vector<TH1D*> > h1_W_MVec, h1_W_PtVec, h1_W_EtaVec;
std::vector<std::vector<TH1D*> > h1_W_dausDRVec;
std::vector<std::vector<TH2D*> > h2_W_dausDR_versus_W_MVec;
std::vector<std::vector<TH1D*> > h1_b_MVec, h1_b_PtVec, h1_b_EtaVec;
std::vector<std::vector<TH1D*> > h1_ratio_mW_over_mTopVec;
std::vector<std::vector<TH2D*> > h2_mW_versus_mTopVec, h2_ratio_mW_over_mTop_versus_mTopVec;
std::vector<std::vector<TH1D*> > h1_mTtopVec, h1_mTWVec, h1_mTbVec;
std::vector<std::vector<TH1D*> > h1_deltaR12Vec, h1_deltaR1bVec, h1_deltaR2bVec;
std::vector<std::vector<TH2D*> > h2_deltaR12_vs_deltaR1bVec, h2_deltaR12_vs_deltaR2bVec, h2_deltaR1b_vs_deltaR2bVec;
std::vector<std::vector<TH1D*> > h1_relPt12Vec, h1_relPt1bVec, h1_relPt2bVec;
std::vector<std::vector<TH2D*> > h2_relPt12_vs_relPt1bVec, h2_relPt12_vs_relPt2bVec, h2_relPt1b_vs_relPt2bVec;

std::vector<TH1D*> h1_MT2_nTopsEQ2Vec, h1_MT2_TopAndbLept_nTopsEQ2Vec, h1_minDphiLeptMET_nTopsEQ2Vec;
std::vector<TH2D*> h2_mTW1_vs_mTtop1_nTopsEQ2Vec, h2_mTb1_vs_mTtop1_nTopsEQ2Vec, h2_mTb1_vs_mTW1_nTopsEQ2Vec;
std::vector<TH2D*> h2_mTW2_vs_mTtop2_nTopsEQ2Vec, h2_mTb2_vs_mTtop2_nTopsEQ2Vec, h2_mTb2_vs_mTW2_nTopsEQ2Vec;
std::vector<TH2D*> h2_mTtop2_vs_mTtop1_nTopsEQ2Vec, h2_mTW2_vs_mTtop1_nTopsEQ2Vec, h2_mTb2_vs_mTtop1_nTopsEQ2Vec;
std::vector<TH2D*> h2_mTtop2_vs_mTW1_nTopsEQ2Vec, h2_mTW2_vs_mTW1_nTopsEQ2Vec, h2_mTb2_vs_mTW1_nTopsEQ2Vec;
std::vector<TH2D*> h2_mTtop2_vs_mTb1_nTopsEQ2Vec, h2_mTW2_vs_mTb1_nTopsEQ2Vec, h2_mTb2_vs_mTb1_nTopsEQ2Vec;

char names[200], dispt[200];

TStopwatch timer;

void draw1DallINone(TCanvas *cs, const int lastPadIdx, const std::vector<TH1D*> &h1_inputVec, const int nRebin =1, const TString optDrawStr ="");
void draw2DallINone(TCanvas *cs, const int lastPadIdx, const std::vector<TH2D*> &h2_inputVec, const TString optDrawStr = "");

bool find_mother(int momIdx, int dauIdx, const std::vector<int> &genDecayIdxVec, const std::vector<int> &genDecayMomIdxVec);
int find_idx(int genIdx, const std::vector<int> &genDecayIdxVec);
double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec);

void drawOverFlowBin(TH1 *histToAdjust){
   int nbins = histToAdjust->GetXaxis()->GetNbins();
   double overflow = histToAdjust->GetBinContent(nbins+1);
   double lastCont = histToAdjust->GetBinContent(nbins);
   histToAdjust->SetBinContent(nbins, overflow+lastCont);
}

void passBaselineFunc(NTupleReader &tr){
   bool passBaseline = true;

// Form TLorentzVector of MET
   TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);

// Calculate number of leptons
   int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
   int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
   int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);

// Calculate number of jets and b-tagged jets
   int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
   int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
   int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
   int cntNJetsPt30      = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);

// Calculate deltaPhi
   std::vector<double> * dPhiVec = new std::vector<double>();
   (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);

// Prepare jets and b-tag working points for top tagger
   std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
   AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
   if( debug ) std::cout<<"\njetsLVec_forTagger->size : "<<jetsLVec_forTagger->size()<<"  recoJetsBtag_forTagger->size : "<<recoJetsBtag_forTagger->size()<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass lepton veto?
   bool passLeptVeto = true;
   if( nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passLeptVeto = false; }
   if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passLeptVeto = false; }
// Isolated track veto is disabled for now
   if( doIsoTrksVeto && nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passLeptVeto = false; }
   if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass number of jets?
   bool passnJets = true;
   if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passnJets = false; }
   if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passnJets = false; }
   if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass deltaPhi?
   bool passdPhis = true;
   if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passdPhis = false; }
   if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass number of b-tagged jets?
   bool passBJets = true;
   if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
   if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass the baseline MET requirement?
   bool passMET = true;
   if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passMET = false; }
   if( debug ) std::cout<<"met : "<<tr.getVar<double>("met")<<"  defaultMETcut : "<<AnaConsts::defaultMETcut<<"  passBaseline : "<<passBaseline<<std::endl;

// Calculate top tagger related variables. 
// Note that to save speed, only do the calculation after previous base line requirements.
   int bestTopJetIdx = -1;
   bool remainPassCSVS = false;
   int pickedRemainingCombfatJetIdx = -1;
   double bestTopJetMass = -1, bestTopJetEta = 0;
   int nTopCandSortedCnt = 0;
   double MT2 = -1;
   double mTcomb = -1;

   if( passBaseline && cntNJetsPt30 >= AnaConsts::nJetsSel ){
      type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
      bestTopJetIdx = type3Ptr->bestTopJetIdx;
      remainPassCSVS = type3Ptr->remainPassCSVS;
      pickedRemainingCombfatJetIdx = type3Ptr->pickedRemainingCombfatJetIdx;
      if( bestTopJetIdx != -1 ){ bestTopJetMass = type3Ptr->bestTopJetLVec.M(); bestTopJetEta = type3Ptr->bestTopJetLVec.Eta(); }

      nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
      MT2 = type3Ptr->MT2;
      mTcomb = type3Ptr->mTbJet + 0.5*type3Ptr->mTbestTopJet;
   }

// Pass top tagger requirement?
   bool passTagger = true;

// bestTopJetIdx != -1 means at least 1 top candidate!
   if( nTopCandSortedCnt ==1 ){
      if( bestTopJetIdx == -1 ){ passBaseline = false; passTagger = false; }
      if( ! remainPassCSVS ){ passBaseline = false; passTagger = false; }
      if( pickedRemainingCombfatJetIdx == -1 && jetsLVec_forTagger->size()>=6 ){ passBaseline = false; passTagger = false; }
      if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ && bestTopJetEta > AnaConsts::lowTopEtaCut_ && bestTopJetEta < AnaConsts::highTopEtaCut_) ){ passBaseline = false; passTagger = false; }
   }

   if( debug ) std::cout<<"bestTopJetidx : "<<bestTopJetIdx<<"  remainPassCSVS : "<<remainPassCSVS<<"  pickedRemainingCombfatJetIdx : "<<pickedRemainingCombfatJetIdx<<"  bestTopJetMass : "<<bestTopJetMass<<"  passBaseline : "<<passBaseline<<std::endl;

// Register all the calculated variables
   tr.registerDerivedVar("nMuons_CUT", nMuons);
   tr.registerDerivedVar("nElectrons_CUT", nElectrons);
   tr.registerDerivedVar("nIsoTrks_CUT", nIsoTrks);

   tr.registerDerivedVar("cntNJetsPt50Eta24", cntNJetsPt50Eta24);
   tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);

   tr.registerDerivedVec("dPhiVec", dPhiVec);

   tr.registerDerivedVar("cntCSVS", cntCSVS);

   tr.registerDerivedVec("jetsLVec_forTagger", jetsLVec_forTagger);
   tr.registerDerivedVec("recoJetsBtag_forTagger", recoJetsBtag_forTagger);

   tr.registerDerivedVar("cntNJetsPt30", cntNJetsPt30);

   tr.registerDerivedVar("bestTopJetIdx", bestTopJetIdx);
   tr.registerDerivedVar("remainPassCSVS", remainPassCSVS);
   tr.registerDerivedVar("pickedRemainingCombfatJetIdx", pickedRemainingCombfatJetIdx);
   tr.registerDerivedVar("bestTopJetMass", bestTopJetMass);
   tr.registerDerivedVar("bestTopJetEta", bestTopJetEta);

   tr.registerDerivedVar("passLeptVeto", passLeptVeto);
   tr.registerDerivedVar("passnJets", passnJets);
   tr.registerDerivedVar("passdPhis", passdPhis);
   tr.registerDerivedVar("passBJets", passBJets);
   tr.registerDerivedVar("passMET", passMET);
   tr.registerDerivedVar("passTagger", passTagger);
   tr.registerDerivedVar("passBaseline", passBaseline);

   if( debug ) std::cout<<"nTopCandSortedCnt : "<<nTopCandSortedCnt<<"  MT2 : "<<MT2<<"  mTcomb : "<<mTcomb<<"  passBaseline : "<<passBaseline<<std::endl;

   tr.registerDerivedVar("nTopCandSortedCnt", nTopCandSortedCnt);
   tr.registerDerivedVar("MT2", MT2);
   tr.registerDerivedVar("mTcomb", mTcomb);

   double HT = AnaFunctions::calcHT(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
   tr.registerDerivedVar("HT", HT);

   if( debug ) std::cout<<"passBaseline : "<<passBaseline<<"  passBaseline : "<<passBaseline<<std::endl;
}

void anaFunc(NTupleReader *tr, std::vector<TTree *> treeVec, const std::vector<std::string> &subSampleKeysVec, const std::string sampleKeyString="ttbar", int verbose=0){
  TString sampleKeyStringT(sampleKeyString);

  keyStringCachedVec.push_back(sampleKeyString);
  double sampleScaleMC = 1.0; int sampleColor = 1;
  for(int ib=0; ib<nMC; ib++){
     TString permcStrT(mcStr[ib]);
     if( permcStrT.Contains(sampleKeyString) ) { sampleColor = colors[ib]; }
  }
  scaleMCCachedVec.push_back( sampleScaleMC );
  colorCachedVec.push_back( sampleColor );

  cnt_passLeptVeto_WeightedScaledMC = 0; cnt_passLeptVeto_WeightedErrorScaledMC = 0;
  cnt_passnJets_WeightedScaledMC = 0; cnt_passnJets_WeightedErrorScaledMC = 0;
  cnt_passdPhis_WeightedScaledMC = 0; cnt_passdPhis_WeightedErrorScaledMC = 0;
  cnt_passBJets_WeightedScaledMC = 0; cnt_passBJets_WeightedErrorScaledMC = 0;
  cnt_passMET_WeightedScaledMC = 0; cnt_passMET_WeightedErrorScaledMC = 0;
  cnt_passTagger_WeightedScaledMC = 0; cnt_passTagger_WeightedErrorScaledMC = 0;
  cnt_passBaseline_WeightedScaledMC = 0; cnt_passBaseline_WeightedErrorScaledMC = 0;
  cnt_passMore_WeightedScaledMC = 0; cnt_passMore_WeightedErrorScaledMC = 0;

  cnt_onelepton_WeightedScaledMC = 0; cnt_onelepton_WeightedErrorScaledMC = 0;
  cnt_Wtauhad_WeightedScaledMC = 0; cnt_Wtauhad_WeightedErrorScaledMC = 0;
  cnt_dileptons_WeightedScaledMC = 0; cnt_dileptons_WeightedErrorScaledMC = 0;
  cnt_diWtauhads_WeightedScaledMC = 0; cnt_diWtauhads_WeightedErrorScaledMC = 0;
  cnt_dileptons_inc_Wtauhad_WeightedScaledMC = 0; cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC = 0;
  cnt_allHad_WeightedScaledMC = 0; cnt_allHad_WeightedErrorScaledMC = 0;

  TPrincipal * pca = new TPrincipal(2, "ND"); pcaVec.push_back(pca);
  keyWordVec.push_back(sampleKeyStringT);

  TH1D * h1_cutFlow = new TH1D(sampleKeyStringT+"_h1_cutFlow", sampleKeyStringT+": cut flow table", 20, 0, 20); h1_cutFlow->SetBit(TH1::kCanRebin); h1_cutFlowVec.push_back((TH1D*)h1_cutFlow->Clone());
  TH1D * h1_cutFlow_aux = new TH1D(sampleKeyStringT+"_h1_cutFlow_aux", sampleKeyStringT+": more cut flow table", 20, 0, 20); h1_cutFlow_aux->SetBit(TH1::kCanRebin); h1_cutFlow_auxVec.push_back((TH1D*)h1_cutFlow_aux->Clone());

  TH1D * h1_pt_genWdau1_nrJets_EQ1 = new TH1D(sampleKeyStringT+"_h1_pt_genWdau1_nrJets_EQ1", sampleKeyStringT+": pt of first W dau for nrJets = 1; P_{T} (GeV)", 100, 0, 200); h1_pt_genWdau1_nrJets_EQ1Vec.push_back((TH1D*)h1_pt_genWdau1_nrJets_EQ1->Clone());
  TH1D * h1_eta_genWdau1_nrJets_EQ1 = new TH1D(sampleKeyStringT+"_h1_eta_genWdau1_nrJets_EQ1", sampleKeyStringT+": #eta of first W dau for nrJets = 1; #eta", 100, -5, 5); h1_eta_genWdau1_nrJets_EQ1Vec.push_back((TH1D*)h1_eta_genWdau1_nrJets_EQ1->Clone());

  TH1D * h1_pt_genWdau2_nrJets_EQ1 = new TH1D(sampleKeyStringT+"_h1_pt_genWdau2_nrJets_EQ1", sampleKeyStringT+": pt of second W dau for nrJets = 1; P_{T} (GeV)", 100, 0, 200); h1_pt_genWdau2_nrJets_EQ1Vec.push_back((TH1D*)h1_pt_genWdau2_nrJets_EQ1->Clone());
  TH1D * h1_eta_genWdau2_nrJets_EQ1 = new TH1D(sampleKeyStringT+"_h1_eta_genWdau2_nrJets_EQ1", sampleKeyStringT+": #eta of second W dau for nrJets =1; #eta", 100, -5, 5); h1_eta_genWdau2_nrJets_EQ1Vec.push_back((TH1D*)h1_eta_genWdau2_nrJets_EQ1->Clone());

  TH1D * h1_mass_rbJet_mtch0_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_mass_rbJet_mtch0_nrJets_EQ2", sampleKeyStringT+" match0 : mass of rbJet for nrJets =2; M (GeV)", 100, 0, 300); h1_mass_rbJet_mtch0_nrJets_EQ2Vec.push_back((TH1D*)h1_mass_rbJet_mtch0_nrJets_EQ2->Clone());
  TH1D * h1_mass_rbJet_mtch1_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_mass_rbJet_mtch1_nrJets_EQ2", sampleKeyStringT+" match1 : mass of rbJet for nrJets =2; M (GeV)", 100, 0, 300); h1_mass_rbJet_mtch1_nrJets_EQ2Vec.push_back((TH1D*)h1_mass_rbJet_mtch1_nrJets_EQ2->Clone());

  TH1D * h1_dR_rbJet_mtch0_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_dR_rbJet_mtch0_nrJets_EQ2", sampleKeyStringT+" match0 : dR of rbJet for nrJets =2; #DeltaR", 100, 0, 5); h1_dR_rbJet_mtch0_nrJets_EQ2Vec.push_back((TH1D*)h1_dR_rbJet_mtch0_nrJets_EQ2->Clone());
  TH1D * h1_dR_rbJet_mtch1_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_dR_rbJet_mtch1_nrJets_EQ2", sampleKeyStringT+" match1 : dR of rbJet for nrJets =2; #DeltaR", 100, 0, 5); h1_dR_rbJet_mtch1_nrJets_EQ2Vec.push_back((TH1D*)h1_dR_rbJet_mtch1_nrJets_EQ2->Clone());

  TH1D * h1_dPhi_rJet_met_mtch0_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_dPhi_rJet_met_mtch0_nrJets_EQ2", sampleKeyStringT+" match0 : dPhi of rJet and met for nrJets =2; #Delta#phi", 100, -3.2, 3.2); h1_dPhi_rJet_met_mtch0_nrJets_EQ2Vec.push_back((TH1D*)h1_dPhi_rJet_met_mtch0_nrJets_EQ2->Clone());
  TH1D * h1_dPhi_rJet_met_mtch1_nrJets_EQ2 = new TH1D(sampleKeyStringT+"_h1_dPhi_rJet_met_mtch1_nrJets_EQ2", sampleKeyStringT+" match1 : dPhi of rJet and met for nrJets =2; #Delta#phi", 100, -3.2, 3.2); h1_dPhi_rJet_met_mtch1_nrJets_EQ2Vec.push_back((TH1D*)h1_dPhi_rJet_met_mtch1_nrJets_EQ2->Clone());

  TH1D * h1_dPhi_rJet_met = new TH1D(sampleKeyStringT+"_h1_dPhi_rJet_met", sampleKeyStringT+": #Delta#phi(rJet, met); #Delta#phi(rJet, met)", 100, -3.2, 3.2); h1_dPhi_rJet_metVec.push_back((TH1D*)h1_dPhi_rJet_met->Clone());
  TH1D * h1_dPhi_rJet_b = new TH1D(sampleKeyStringT+"_h1_dPhi_rJet_b", sampleKeyStringT+": #Delta#phi(rJet, b); #Delta#phi(rJet, b)", 100, -3.2, 3.2); h1_dPhi_rJet_bVec.push_back((TH1D*)h1_dPhi_rJet_b->Clone());

  TH1D * h1_dPhi_top_met = new TH1D(sampleKeyStringT+"_h1_dPhi_top_met", sampleKeyStringT+": #Delta#phi(top, met); #Delta#phi(top, met)", 100, -3.2, 3.2); h1_dPhi_top_metVec.push_back((TH1D*)h1_dPhi_top_met->Clone());
  TH1D * h1_minDphi_tJet_met = new TH1D(sampleKeyStringT+"_h1_minDphi_tJet_met", sampleKeyStringT+": #Delta#phi(tJet, met); #Delta#phi(tJet, met)", 100, -3.2, 3.2); h1_minDphi_tJet_metVec.push_back((TH1D*)h1_minDphi_tJet_met->Clone());

  TH1D * h1_MT_tJet = new TH1D(sampleKeyStringT+"_h1_MT_tJet", sampleKeyStringT+": MT(tJet, met); MT(tJet, met) (GeV)", 100, 0, 1000); h1_MT_tJetVec.push_back((TH1D*)h1_MT_tJet->Clone());
  TH1D * h1_MT_tDiJets = new TH1D(sampleKeyStringT+"_h1_MT_tDiJetsVec", sampleKeyStringT+": MT(tDiJetsVec, met); MT(tDiJetsVec, met) (GeV)", 100, 0, 1000); h1_MT_tDiJetsVec.push_back((TH1D*)h1_MT_tDiJets->Clone());

  TH2D * h2_dalitz = new TH2D(sampleKeyStringT+"_h2_dalitz", sampleKeyStringT+": mb1^{2} versus mb2^{2}; mb2^{2} (GeV^{2}); mb1^{2} (GeV^{2}", 100, 0, 200, 100, 0, 200); h2_dalitzVec.push_back((TH2D*)h2_dalitz->Clone());
  TH2D * h2_gen_m23overm123vsarctanm13overm12 = new TH2D(sampleKeyStringT+"_h2_gen_m23overm123vsarctanm13overm12", sampleKeyStringT+"  gen level : m23/m123 versus arctan m13/m12; arctan m13/m12; m23/m123", 100, 0, 1.5, 100, 0, 1.0); h2_gen_m23overm123vsarctanm13overm12Vec.push_back((TH2D*)h2_gen_m23overm123vsarctanm13overm12->Clone());
  TH2D * h2_m23overm123vsarctanm13overm12 = new TH2D(sampleKeyStringT+"_h2_m23overm123vsarctanm13overm12", sampleKeyStringT+": m23/m123 versus arctan m13/m12; arctan m13/m12; m23/m123", 100, 0, 1.5, 100, 0, 1.0); h2_m23overm123vsarctanm13overm12Vec.push_back((TH2D*)h2_m23overm123vsarctanm13overm12->Clone());

  TH2D * h2_dPhi_rJet_b_vs_met = new TH2D(sampleKeyStringT+"_h2_dPhi_rJet_b_vs_met", sampleKeyStringT+": #Delta#phi(rJet, b) vs. #Delta#phi(rJet, met); #Delta#phi(rJet, met); #Delta#phi(rJet, b)", 100, -3.2, 3.2, 100, -3.2, 3.2); h2_dPhi_rJet_b_vs_metVec.push_back((TH2D*)h2_dPhi_rJet_b_vs_met->Clone());

  TH1D * h1_dR_gen_reco_top = new TH1D(sampleKeyStringT+"_h1_dR_gen_reco_top", sampleKeyStringT+": #DeltaR(gen, reco) of top; #DeltaR(gen, reco)", 100, 0, 5.0); h1_dR_gen_reco_topVec.push_back((TH1D*)h1_dR_gen_reco_top->Clone());

  TH1D * h1_dPhi_sum_lept_nu_met = new TH1D(sampleKeyStringT+"_h1_dPhi_sum_lept_nu_met", sampleKeyStringT+": #Delta#phi(lept+nu, met); #Delta#phi(lept+nu, met)", 100, -1.5, 1.5); h1_dPhi_sum_lept_nu_metVec.push_back((TH1D*)h1_dPhi_sum_lept_nu_met->Clone());
  TH1D * h1_dPhi_lept_met = new TH1D(sampleKeyStringT+"_h1_dPhi_lept_met", sampleKeyStringT+": #Delta#phi(lept, met); #Delta#phi(lept, met)", 100, -3.2, 3.2); h1_dPhi_lept_metVec.push_back((TH1D*)h1_dPhi_lept_met->Clone());
  TH1D * h1_dPhi_nu_met = new TH1D(sampleKeyStringT+"_h1_dPhi_nu_met", sampleKeyStringT+": #Delta#phi(nu, met); #Delta#phi(nu, met)", 100, -1.5, 1.5); h1_dPhi_nu_metVec.push_back((TH1D*)h1_dPhi_nu_met->Clone());

  TH1D * h1_dR_b_top = new TH1D(sampleKeyStringT+"_h1_dR_b_top", sampleKeyStringT+": #DeltaR(b, top); #DeltaR(b, top)", 100, 0, 5.0); h1_dR_b_topVec.push_back((TH1D*)h1_dR_b_top->Clone());
  TH1D * h1_dPhi_b_top = new TH1D(sampleKeyStringT+"_h1_dPhi_b_top", sampleKeyStringT+": #Delta#phi(b, top); #Delta#phi(b, top)", 100, -3.2, 3.2); h1_dPhi_b_topVec.push_back((TH1D*)h1_dPhi_b_top->Clone());

  TH2D * h2_mt_b_vs_dPhi_b_top = new TH2D(sampleKeyStringT+"_h2_mt_b_vs_dPhi_b_top", sampleKeyStringT+": M_{T}^{b} vs. #Delta#phi(b, top); #Delta#phi(b, top); M_{T}^{b} (GeV)", 100, -3.2, 3.2, 100, 0, 500); h2_mt_b_vs_dPhi_b_topVec.push_back((TH2D*)h2_mt_b_vs_dPhi_b_top->Clone());
                               
  TH2D * h2_evtCnt_nbJets_vs_nTops = new TH2D(sampleKeyStringT+"_h2_evtCnt_nbJets_vs_nTops", sampleKeyStringT+": event counts nbJets versus nTops; nTops; nbJets", 4, 0, 4, 3, 1, 4); h2_evtCnt_nbJets_vs_nTopsVec.push_back((TH2D*) h2_evtCnt_nbJets_vs_nTops->Clone());

  TH1D * h1_nComb_lept_brJet = new TH1D(sampleKeyStringT+"_h1_nComb_lept_brJet", sampleKeyStringT+": lept-like  number of combination of b and rjet; nComb", 5, 0, 5); h1_nComb_lept_brJetVec.push_back((TH1D*) h1_nComb_lept_brJet->Clone());
  TH1D * h1_nComb_had_brJet = new TH1D(sampleKeyStringT+"_h1_nComb_had_brJet", sampleKeyStringT+": had-like  number of combination of b and rjet; nComb", 5, 0, 5); h1_nComb_had_brJetVec.push_back((TH1D*) h1_nComb_had_brJet->Clone());

  TH1D * h1_nComb_top_brJets = new TH1D(sampleKeyStringT+"_h1_nComb_top_brJets", sampleKeyStringT+": top jet in the remaining system; nComb", 5, 0, 5); h1_nComb_top_brJetsVec.push_back((TH1D*) h1_nComb_top_brJets->Clone());

  TH2D * h2_nComb_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_nComb_lept_vs_had_brJet", sampleKeyStringT+": lept vs. had. of combination of b and rjet; nComb_{had}; nComb_{lept}", 5, 0, 5, 5, 0, 5); h2_nComb_lept_vs_had_brJetVec.push_back((TH2D*)h2_nComb_lept_vs_had_brJet->Clone());

  TH1D * h1_dPhi_genTops = new TH1D(sampleKeyStringT+"_h1_dPhi_genTops", sampleKeyStringT+": #Delta#phi(genTop1, genTop2); #Delta#phi", 100, -3.2, 3.2); h1_dPhi_genTopsVec.push_back((TH1D*) h1_dPhi_genTops->Clone());
  TH1D * h1_dPhi_genbs = new TH1D(sampleKeyStringT+"_h1_dPhi_genbs", sampleKeyStringT+": #Delta#phi(genb1, genb2); #Delta#phi", 100, -3.2, 3.2); h1_dPhi_genbsVec.push_back((TH1D*) h1_dPhi_genbs->Clone());
  TH1D * h1_dPhi_genb_genTop = new TH1D(sampleKeyStringT+"_h1_dPhi_genb_genTop", sampleKeyStringT+": #Delta#phi(genTop1, genb2); #Delta#phi", 100, -3.2, 3.2); h1_dPhi_genb_genTopVec.push_back((TH1D*) h1_dPhi_genb_genTop->Clone());

  TH1D * h1_MT_lept_brJet = new TH1D(sampleKeyStringT+"_h1_MT_lept_brJet", sampleKeyStringT+": lept-like  MT(brJet, met); MT(brJet, met)", 100, 0, 1000); h1_MT_lept_brJetVec.push_back((TH1D*) h1_MT_lept_brJet->Clone());
  TH1D * h1_MT_had_brJet = new TH1D(sampleKeyStringT+"_h1_MT_had_brJet", sampleKeyStringT+": had-like  MT(brJet, met); MT(brJet, met)", 100, 0, 1000); h1_MT_had_brJetVec.push_back((TH1D*) h1_MT_had_brJet->Clone());

  TH1D * h1_MT_bJet = new TH1D(sampleKeyStringT+"_h1_MT_bJet", sampleKeyStringT+": MT(bJet, met); MT(bJet, met)", 100, 0, 1000); h1_MT_bJetVec.push_back((TH1D*) h1_MT_bJet->Clone());

  TH1D * h1_mass_rTop = new TH1D(sampleKeyStringT+"_h1_mass_rTop", sampleKeyStringT+": mass of rTop; M(rTop) (GeV)", 100, 0, 300); h1_mass_rTopVec.push_back((TH1D*) h1_mass_rTop->Clone());

  TH1D * h1_ori_MT_lept_brJet = new TH1D(sampleKeyStringT+"_h1_ori_MT_lept_brJet", sampleKeyStringT+": original lept-like  MT(brJet, met); MT(brJet, met)", 100, 0, 1000); h1_ori_MT_lept_brJetVec.push_back((TH1D*) h1_ori_MT_lept_brJet->Clone());
  TH1D * h1_ori_MT_had_brJet = new TH1D(sampleKeyStringT+"_h1_ori_MT_had_brJet", sampleKeyStringT+": original had-like  MT(brJet, met); MT(brJet, met)", 100, 0, 1000); h1_ori_MT_had_brJetVec.push_back((TH1D*) h1_ori_MT_had_brJet->Clone());

  TH1D * h1_ori_MT2_lept_brJet = new TH1D(sampleKeyStringT+"_h1_ori_MT2_lept_brJet", sampleKeyStringT+": original lept-like  MT2; MT2", 100, 0, 1000); h1_ori_MT2_lept_brJetVec.push_back((TH1D*) h1_ori_MT2_lept_brJet->Clone());
  TH1D * h1_ori_MT2_had_brJet = new TH1D(sampleKeyStringT+"_h1_ori_MT2_had_brJet", sampleKeyStringT+": original had-like  MT2; MT2", 100, 0, 1000); h1_ori_MT2_had_brJetVec.push_back((TH1D*) h1_ori_MT2_had_brJet->Clone());

  TH1D * h1_dPhi_sum_MET = new TH1D(sampleKeyStringT+"_h1_dPhi_sum_MET", sampleKeyStringT+": #Delta#phi(-sum, met); #Delta#phi", 100, -3.2, 3.2); h1_dPhi_sum_METVec.push_back((TH1D*)h1_dPhi_sum_MET->Clone());
  TH1D * h1_MT_sum = new TH1D(sampleKeyStringT+"_h1_MT_sum", sampleKeyStringT+": MT(-sum, met); MT", 100, 0, 1000); h1_MT_sumVec.push_back((TH1D*)h1_MT_sum->Clone());

  TH1D * h1_MT2_lept_brJet = new TH1D(sampleKeyStringT+"_h1_MT2_lept_brJet", sampleKeyStringT+": lept-like  MT2; MT2", 100, 0, 1000); h1_MT2_lept_brJetVec.push_back((TH1D*) h1_MT2_lept_brJet->Clone());
  TH1D * h1_MT2_had_brJet = new TH1D(sampleKeyStringT+"_h1_MT2_had_brJet", sampleKeyStringT+": had-like  MT2; MT2", 100, 0, 1000); h1_MT2_had_brJetVec.push_back((TH1D*) h1_MT2_had_brJet->Clone());

  TH2D * h2_MT2_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT2_lept_vs_had_brJet", sampleKeyStringT+": MT2 lept-like versus had-like; MT2(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT2_lept_vs_had_brJetVec.push_back((TH2D*)h2_MT2_lept_vs_had_brJet->Clone());
  TH1D * h1_MT2_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_MT2_same_lept_had_brJet", sampleKeyStringT+": MT2 same lept-like and had-like; MT2(same)", 100, 0, 1000); h1_MT2_same_lept_had_brJetVec.push_back((TH1D*)h1_MT2_same_lept_had_brJet->Clone());
  TH1D * h1_MT_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_MT_same_lept_had_brJet", sampleKeyStringT+": MT same lept-like and had-like; MT(same)", 100, 0, 1000); h1_MT_same_lept_had_brJetVec.push_back((TH1D*)h1_MT_same_lept_had_brJet->Clone());
  TH1D * h1_mTcomb_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_mTcomb_same_lept_had_brJet", sampleKeyStringT+": mTcomb same lept-like and had-like; mTcomb(same)", 100, 0, 1500); h1_mTcomb_same_lept_had_brJetVec.push_back((TH1D*)h1_mTcomb_same_lept_had_brJet->Clone());

  TH2D * h2_MT2_lept_vs_had_aft_MT_cuts_brJet = new TH2D(sampleKeyStringT+"_h2_MT2_lept_vs_had_aft_MT_cuts_brJet", sampleKeyStringT+" aft MT cuts : MT2 lept-like versus had-like; MT2(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT2_lept_vs_had_aft_MT_cuts_brJetVec.push_back((TH2D*)h2_MT2_lept_vs_had_aft_MT_cuts_brJet->Clone());
  TH1D * h1_MT2_same_lept_had_aft_MT_cuts_brJets = new TH1D(sampleKeyStringT+"_h1_MT2_same_lept_had_aft_MT_cuts_brJets", sampleKeyStringT+" aft MT cuts : MT2 same lept-like and had-like; MT2(same)", 100, 0, 1000); h1_MT2_same_lept_had_aft_MT_cuts_brJetsVec.push_back((TH1D*)h1_MT2_same_lept_had_aft_MT_cuts_brJets->Clone());

  TH2D * h2_MT_lept_vs_MT_bJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_MT_bJet", sampleKeyStringT+": lept-like MT versus bJet MT; MT(bJet, met); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_lept_vs_MT_bJetVec.push_back((TH2D*)h2_MT_lept_vs_MT_bJet->Clone());
  TH2D * h2_MT_bJet_vs_mTcomb_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_bJet_vs_mTcomb_had_brJet", sampleKeyStringT+": MT(bJet, met) versus had-like mTcomb; mTcomb(had); MT(bJet, met)", 100, 0, 1000, 100, 0, 1000); h2_MT_bJet_vs_mTcomb_had_brJetVec.push_back((TH2D*)h2_MT_bJet_vs_mTcomb_had_brJet->Clone());

  TH2D * h2_MT_lept_vs_top_brJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_top_brJet", sampleKeyStringT+": lept-like MT versus top MT; MT(top); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_lept_vs_top_brJetVec.push_back((TH2D*)h2_MT_lept_vs_top_brJet->Clone());
  TH2D * h2_mTcomb_vs_MT_lept_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_vs_MT_lept_brJet", sampleKeyStringT+": lept-like mTcomb versus MT; MT(lept); mTcomb(lept)", 100, 0, 1000, 100, 0, 1500); h2_mTcomb_vs_MT_lept_brJetVec.push_back((TH2D*)h2_mTcomb_vs_MT_lept_brJet->Clone());
  TH2D * h2_MT_vs_MT2_lept_brJet = new TH2D(sampleKeyStringT+"_h2_MT_vs_MT2_lept_brJet", sampleKeyStringT+": lept-like MT versus MT2; MT2(lept); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_vs_MT2_lept_brJetVec.push_back((TH2D*)h2_MT_vs_MT2_lept_brJet->Clone());
  TH2D * h2_mTcomb_vs_MT2_lept_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_vs_MT2_lept_brJet", sampleKeyStringT+": lept-like mTcomb versus MT2; MT2(lept); mTcomb(lept)", 100, 0, 1000, 100, 0, 1500); h2_mTcomb_vs_MT2_lept_brJetVec.push_back((TH2D*)h2_mTcomb_vs_MT2_lept_brJet->Clone());
  TH2D * h2_MT_had_vs_top_brJet = new TH2D(sampleKeyStringT+"_h2_MT_had_vs_top_brJet", sampleKeyStringT+": had-like MT versus top MT; MT(top); MT(had)", 100, 0, 1000, 100, 0, 1000); h2_MT_had_vs_top_brJetVec.push_back((TH2D*)h2_MT_had_vs_top_brJet->Clone());
  TH2D * h2_aft_PCA_MT_had_vs_top_brJet = new TH2D(sampleKeyStringT+"_h2_aft_PCA_MT_had_vs_top_brJet", sampleKeyStringT+" aft PCA : had-like MT versus top MT; MT(top); MT(had)", 100, -5, 5, 100, -5, 5); h2_aft_PCA_MT_had_vs_top_brJetVec.push_back((TH2D*)h2_aft_PCA_MT_had_vs_top_brJet->Clone());
  TH2D * h2_mTcomb_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_vs_MT_had_brJet", sampleKeyStringT+": had-like mTcomb versus MT; MT(had); mTcomb(had)", 100, 0, 1000, 100, 0, 1500); h2_mTcomb_vs_MT_had_brJetVec.push_back((TH2D*)h2_mTcomb_vs_MT_had_brJet->Clone());
  TH2D * h2_MT_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_vs_MT2_had_brJet", sampleKeyStringT+": had-like MT versus MT2; MT2(had); MT(had)", 100, 0, 1000, 100, 0, 1000); h2_MT_vs_MT2_had_brJetVec.push_back((TH2D*)h2_MT_vs_MT2_had_brJet->Clone());
  TH2D * h2_ori_MT2_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_ori_MT2_vs_MT2_had_brJet", sampleKeyStringT+": had-like MT2 versus ori MT2; MT2(had); MT2(ori)", 100, 0, 1000, 100, 0, 1000); h2_ori_MT2_vs_MT2_had_brJetVec.push_back((TH2D*)h2_ori_MT2_vs_MT2_had_brJet->Clone());
  TH2D * h2_mTcomb_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_vs_MT2_had_brJet", sampleKeyStringT+": had-like mTcomb versus MT2; MT2(had); mTcomb(had)", 100, 0, 1000, 100, 0, 1500); h2_mTcomb_vs_MT2_had_brJetVec.push_back((TH2D*)h2_mTcomb_vs_MT2_had_brJet->Clone());
  TH2D * h2_MT_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_had_brJet", sampleKeyStringT+": MT lept-like versus had-like; MT(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_lept_vs_had_brJetVec.push_back((TH2D*)h2_MT_lept_vs_had_brJet->Clone());
  TH2D * h2_MT_lept_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_MT2_had_brJet", sampleKeyStringT+": MT lept-like versus MT2 had-like; MT2(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_lept_vs_MT2_had_brJetVec.push_back((TH2D*)h2_MT_lept_vs_MT2_had_brJet->Clone());
  TH2D * h2_MT2_lept_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT2_lept_vs_MT_had_brJet", sampleKeyStringT+": MT2 lept-like versus MT had-like; MT(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT2_lept_vs_MT_had_brJetVec.push_back((TH2D*)h2_MT2_lept_vs_MT_had_brJet->Clone());

  TH2D * h2_MT_lept_vs_mTcomb_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_mTcomb_had_brJet", sampleKeyStringT+": MT lept-like versus mTcomb had-like; mTcomb(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_MT_lept_vs_mTcomb_had_brJetVec.push_back((TH2D*)h2_MT_lept_vs_mTcomb_had_brJet->Clone());
  TH2D * h2_MT_lept_vs_aft_PCA_had_brJet = new TH2D(sampleKeyStringT+"_h2_MT_lept_vs_aft_PCA_had_brJet", sampleKeyStringT+": MT lept-like versus aft_PCA had-like; aft_PCA(had); MT(lept)", 100, -5, 5, 100, 0, 1000); h2_MT_lept_vs_aft_PCA_had_brJetVec.push_back((TH2D*)h2_MT_lept_vs_aft_PCA_had_brJet->Clone());
  TH2D * h2_mTcomb_lept_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_lept_vs_MT_had_brJet", sampleKeyStringT+": mTcomb lept-like versus MT had-like; MT(had); mTcomb(lept)", 100, 0, 1000, 100, 0, 1000); h2_mTcomb_lept_vs_MT_had_brJetVec.push_back((TH2D*)h2_mTcomb_lept_vs_MT_had_brJet->Clone());

  TH2D * h2_mTcomb_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_mTcomb_lept_vs_had_brJet", sampleKeyStringT+": mTcomb lept-like versus had-like; mTcomb(had); mTcomb(lept)", 100, 0, 1500, 100, 0, 1500); h2_mTcomb_lept_vs_had_brJetVec.push_back((TH2D*)h2_mTcomb_lept_vs_had_brJet->Clone());

  TH1D * h1_csvs_fakeb = new TH1D(sampleKeyStringT+"_h1_csvs_fakeb", sampleKeyStringT+": csvs of the fake b; csvs", 100, 0, 1); h1_csvs_fakebVec.push_back((TH1D*)h1_csvs_fakeb->Clone());

  TH1D * h1_pt_gentb = new TH1D(sampleKeyStringT+"_h1_pt_gentb", sampleKeyStringT+": pt of gen b from top; p_{T}(GeV)", 100, 0, 300); h1_pt_gentbVec.push_back((TH1D*)h1_pt_gentb->Clone());
  TH1D * h1_eta_gentb = new TH1D(sampleKeyStringT+"_h1_eta_gentb", sampleKeyStringT+": eta of gen b from top; #eta", 100, -5, 5); h1_eta_gentbVec.push_back((TH1D*)h1_eta_gentb->Clone());
  TH1D * h1_pt_genrb = new TH1D(sampleKeyStringT+"_h1_pt_genrb", sampleKeyStringT+": pt of gen b from rtop; p_{T}(GeV)", 100, 0, 300); h1_pt_genrbVec.push_back((TH1D*)h1_pt_genrb->Clone());
  TH1D * h1_eta_genrb = new TH1D(sampleKeyStringT+"_h1_eta_genrb", sampleKeyStringT+": eta of gen b from rtop; #eta", 100, -5, 5); h1_eta_genrbVec.push_back((TH1D*)h1_eta_genrb->Clone());
  TH1D * h1_pt_genrb_match0 = new TH1D(sampleKeyStringT+"_h1_pt_genrb_match0", sampleKeyStringT+"  match0 : pt of gen b from rtop; p_{T}(GeV)", 100, 0, 300); h1_pt_genrb_match0Vec.push_back((TH1D*)h1_pt_genrb_match0->Clone());
  TH1D * h1_eta_genrb_match0 = new TH1D(sampleKeyStringT+"_h1_eta_genrb_match0", sampleKeyStringT+"  match0 : eta of gen b from rtop; #eta", 100, -5, 5); h1_eta_genrb_match0Vec.push_back((TH1D*)h1_eta_genrb_match0->Clone());

// begin adj
  TH1D * h1_adj_MT_lept_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT_lept_brJet", sampleKeyStringT+" adj : lept-like  MT(lept); MT", 100, 0, 1000); h1_adj_MT_lept_brJetVec.push_back((TH1D*) h1_adj_MT_lept_brJet->Clone());
  TH1D * h1_adj_MT_had_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT_had_brJet", sampleKeyStringT+" adj : had-like  MT(had); MT", 100, 0, 1000); h1_adj_MT_had_brJetVec.push_back((TH1D*) h1_adj_MT_had_brJet->Clone());
  TH1D * h1_adj_MT2_lept_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT2_lept_brJet", sampleKeyStringT+" adj : lept-like  MT2; MT2", 100, 0, 1000); h1_adj_MT2_lept_brJetVec.push_back((TH1D*) h1_adj_MT2_lept_brJet->Clone());
  TH1D * h1_adj_MT2_had_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT2_had_brJet", sampleKeyStringT+" adj : had-like  MT2; MT2", 100, 0, 1000); h1_adj_MT2_had_brJetVec.push_back((TH1D*) h1_adj_MT2_had_brJet->Clone());

  TH2D * h2_adj_MT2_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT2_lept_vs_had_brJet", sampleKeyStringT+" adj :  MT2 lept-like versus had-like; MT2(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT2_lept_vs_had_brJetVec.push_back((TH2D*)h2_adj_MT2_lept_vs_had_brJet->Clone());
  TH1D * h1_adj_MT2_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT2_same_lept_had_brJet", sampleKeyStringT+" adj :  MT2 same lept-like and had-like; MT2(same)", 100, 0, 1000); h1_adj_MT2_same_lept_had_brJetVec.push_back((TH1D*)h1_adj_MT2_same_lept_had_brJet->Clone());
  TH1D * h1_adj_MT_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_adj_MT_same_lept_had_brJet", sampleKeyStringT+" adj :  MT same lept-like and had-like; MT(same)", 100, 0, 1000); h1_adj_MT_same_lept_had_brJetVec.push_back((TH1D*)h1_adj_MT_same_lept_had_brJet->Clone());
  TH1D * h1_adj_mTcomb_same_lept_had_brJet = new TH1D(sampleKeyStringT+"_h1_adj_mTcomb_same_lept_had_brJet", sampleKeyStringT+" adj :  mTcomb same lept-like and had-like; mTcomb(same)", 100, 0, 1500); h1_adj_mTcomb_same_lept_had_brJetVec.push_back((TH1D*)h1_adj_mTcomb_same_lept_had_brJet->Clone());

  TH2D * h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJet", sampleKeyStringT+" aft MT cuts : MT2 lept-like versus had-like; MT2(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJetVec.push_back((TH2D*)h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJet->Clone());
  TH1D * h1_adj_MT2_same_lept_had_aft_MT_cuts_brJets = new TH1D(sampleKeyStringT+"_h1_adj_MT2_same_lept_had_aft_MT_cuts_brJets", sampleKeyStringT+" aft MT cuts : MT2 same lept-like and had-like; MT2(same)", 100, 0, 1000); h1_adj_MT2_same_lept_had_aft_MT_cuts_brJetsVec.push_back((TH1D*)h1_adj_MT2_same_lept_had_aft_MT_cuts_brJets->Clone());

  TH2D * h2_adj_MT_lept_vs_top_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_lept_vs_top_brJet", sampleKeyStringT+" adj :  lept-like MT versus top MT; MT(top); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_lept_vs_top_brJetVec.push_back((TH2D*)h2_adj_MT_lept_vs_top_brJet->Clone());
  TH2D * h2_adj_mTcomb_vs_MT_lept_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_vs_MT_lept_brJet", sampleKeyStringT+" adj :  lept-like mTcomb versus MT; MT(lept); mTcomb(lept)", 100, 0, 1000, 100, 0, 1500); h2_adj_mTcomb_vs_MT_lept_brJetVec.push_back((TH2D*)h2_adj_mTcomb_vs_MT_lept_brJet->Clone());
  TH2D * h2_adj_MT_vs_MT2_lept_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_vs_MT2_lept_brJet", sampleKeyStringT+" adj :  lept-like MT versus MT2; MT2(lept); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_vs_MT2_lept_brJetVec.push_back((TH2D*)h2_adj_MT_vs_MT2_lept_brJet->Clone());
  TH2D * h2_adj_mTcomb_vs_MT2_lept_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_vs_MT2_lept_brJet", sampleKeyStringT+" adj :  lept-like mTcomb versus MT2; MT2(lept); mTcomb(lept)", 100, 0, 1000, 100, 0, 1500); h2_adj_mTcomb_vs_MT2_lept_brJetVec.push_back((TH2D*)h2_adj_mTcomb_vs_MT2_lept_brJet->Clone());
  TH2D * h2_adj_MT_had_vs_top_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_had_vs_top_brJet", sampleKeyStringT+" adj :  had-like MT versus top MT; MT(top); MT(had)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_had_vs_top_brJetVec.push_back((TH2D*)h2_adj_MT_had_vs_top_brJet->Clone());
  TH2D * h2_adj_mTcomb_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_vs_MT_had_brJet", sampleKeyStringT+" adj :  had-like mTcomb versus MT; MT(had); mTcomb(had)", 100, 0, 1000, 100, 0, 1500); h2_adj_mTcomb_vs_MT_had_brJetVec.push_back((TH2D*)h2_adj_mTcomb_vs_MT_had_brJet->Clone());
  TH2D * h2_adj_MT_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_vs_MT2_had_brJet", sampleKeyStringT+" adj :  had-like MT versus MT2; MT2(had); MT(had)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_vs_MT2_had_brJetVec.push_back((TH2D*)h2_adj_MT_vs_MT2_had_brJet->Clone());
  TH2D * h2_adj_mTcomb_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_vs_MT2_had_brJet", sampleKeyStringT+" adj :  had-like mTcomb versus MT2; MT2(had); mTcomb(had)", 100, 0, 1000, 100, 0, 1500); h2_adj_mTcomb_vs_MT2_had_brJetVec.push_back((TH2D*)h2_adj_mTcomb_vs_MT2_had_brJet->Clone());
  TH2D * h2_adj_MT_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_lept_vs_had_brJet", sampleKeyStringT+" adj :  MT lept-like versus had-like; MT(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_lept_vs_had_brJetVec.push_back((TH2D*)h2_adj_MT_lept_vs_had_brJet->Clone());
  TH2D * h2_adj_MT_lept_vs_MT2_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_lept_vs_MT2_had_brJet", sampleKeyStringT+" adj :  MT lept-like versus MT2 had-like; MT2(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_lept_vs_MT2_had_brJetVec.push_back((TH2D*)h2_adj_MT_lept_vs_MT2_had_brJet->Clone());
  TH2D * h2_adj_MT2_lept_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT2_lept_vs_MT_had_brJet", sampleKeyStringT+" adj :  MT2 lept-like versus MT had-like; MT(had); MT2(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT2_lept_vs_MT_had_brJetVec.push_back((TH2D*)h2_adj_MT2_lept_vs_MT_had_brJet->Clone());

  TH2D * h2_adj_MT_lept_vs_mTcomb_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_MT_lept_vs_mTcomb_had_brJet", sampleKeyStringT+" adj :  MT lept-like versus mTcomb had-like; mTcomb(had); MT(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_MT_lept_vs_mTcomb_had_brJetVec.push_back((TH2D*)h2_adj_MT_lept_vs_mTcomb_had_brJet->Clone());
  TH2D * h2_adj_mTcomb_lept_vs_MT_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_lept_vs_MT_had_brJet", sampleKeyStringT+" adj :  mTcomb lept-like versus MT had-like; MT(had); mTcomb(lept)", 100, 0, 1000, 100, 0, 1000); h2_adj_mTcomb_lept_vs_MT_had_brJetVec.push_back((TH2D*)h2_adj_mTcomb_lept_vs_MT_had_brJet->Clone());

  TH2D * h2_adj_mTcomb_lept_vs_had_brJet = new TH2D(sampleKeyStringT+"_h2_adj_mTcomb_lept_vs_had_brJet", sampleKeyStringT+" adj :  mTcomb lept-like versus had-like; mTcomb(had); mTcomb(lept)", 100, 0, 1500, 100, 0, 1500); h2_adj_mTcomb_lept_vs_had_brJetVec.push_back((TH2D*)h2_adj_mTcomb_lept_vs_had_brJet->Clone());

// end adj

  TH1D * h1_gen1b_mass = new TH1D(sampleKeyStringT+"_h1_gen1b_mass", sampleKeyStringT+": mass of b quark and closest quark from W; m_{1b} (GeV)", 100, 0, 200); h1_gen1b_massVec.push_back((TH1D*) h1_gen1b_mass->Clone());
  TH1D * h1_gen2b_mass = new TH1D(sampleKeyStringT+"_h1_gen2b_mass", sampleKeyStringT+": mass of b quark and second closest quark from W; m_{2b} (GeV)", 100, 0, 200); h1_gen2b_massVec.push_back((TH1D*) h1_gen2b_mass->Clone());

  TH1D * h1_gen1b_deltaR = new TH1D(sampleKeyStringT+"_h1_gen1b_deltaR", sampleKeyStringT+": deltaR of b quark and closest quark from W; #DeltaR_{1b}", 100, 0, 5); h1_gen1b_deltaRVec.push_back((TH1D*) h1_gen1b_deltaR->Clone());
  TH1D * h1_gen2b_deltaR = new TH1D(sampleKeyStringT+"_h1_gen2b_deltaR", sampleKeyStringT+": deltaR of b quark and second closest quark from W; #DeltaR_{2b}", 100, 0, 5); h1_gen2b_deltaRVec.push_back((TH1D*) h1_gen2b_deltaR->Clone());
  TH1D * h1_gen12_deltaR = new TH1D(sampleKeyStringT+"_h1_gen12_deltaR", sampleKeyStringT+": deltaR of 1st and 2nd quarks from W; #DeltaR_{12}", 100, 0, 5); h1_gen12_deltaRVec.push_back((TH1D*) h1_gen12_deltaR->Clone());

  TH1D * h1_gen1MET_deltaPhi = new TH1D(sampleKeyStringT+"_h1_gen1MET_deltaPhi", sampleKeyStringT+": deltaPhi of 1st quark from W to MET; #Delta#phi_{1MET}", 100, -3.2, 3.2); h1_gen1MET_deltaPhiVec.push_back((TH1D*)h1_gen1MET_deltaPhi->Clone());
  TH1D * h1_gen2MET_deltaPhi = new TH1D(sampleKeyStringT+"_h1_gen2MET_deltaPhi", sampleKeyStringT+": deltaPhi of 2nd quark from W to MET; #Delta#phi_{2MET}", 100, -3.2, 3.2); h1_gen2MET_deltaPhiVec.push_back((TH1D*)h1_gen2MET_deltaPhi->Clone());

  TH1D * h1_Wjets_bTagCats = new TH1D(sampleKeyStringT+"_h1_Wjets_bTagCats", sampleKeyStringT+": b tagged or not for W jets; bTagCats", 2, 0, 2); h1_Wjets_bTagCatsVec.push_back((TH1D*)h1_Wjets_bTagCats->Clone());
  TH1D * h1_Topjets_bTagCats = new TH1D(sampleKeyStringT+"_h1_Topjets_bTagCats", sampleKeyStringT+": b tagged or not for Top jets; bTagCats", 2, 0, 2); h1_Topjets_bTagCatsVec.push_back((TH1D*)h1_Topjets_bTagCats->Clone());

  TH1D * h1_bTagged_Wjets_CSV = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_CSV", sampleKeyStringT+": bTagged W jets : CSV; CSV", 25, AnaConsts::cutCSVS, 1.0); h1_bTagged_Wjets_CSVVec.push_back((TH1D*)h1_bTagged_Wjets_CSV->Clone());

  TH2D * h2_bTagged_Wjets_minDR_genb_vs_minDR_genW = new TH2D(sampleKeyStringT+"_h2_bTagged_Wjets_minDR_genb_vs_minDR_genW", sampleKeyStringT+": bTagged W jets :  minDR(Wjet, genb) versus minDR(Wjet, genW); #DeltaR(Wjet, genW); #DeltaR(Wjet, genb)", 25, 0, 1.5, 25, 0, 1.5); h2_bTagged_Wjets_minDR_genb_vs_minDR_genWVec.push_back((TH2D*)h2_bTagged_Wjets_minDR_genb_vs_minDR_genW->Clone());
  TH2D * h2_bTagged_Topjets_minDR_genb_vs_minDR_genTop = new TH2D(sampleKeyStringT+"_h2_bTagged_Topjets_minDR_genb_vs_minDR_genTop", sampleKeyStringT+": bTagged Top jets :  minDR(Topjet, genb) versus minDR(Topjet, genTop); #DeltaR(Topjet, genTop); #DeltaR(Topjet, genb)", 25, 0, 1.5, 25, 0, 1.5); h2_bTagged_Topjets_minDR_genb_vs_minDR_genTopVec.push_back((TH2D*)h2_bTagged_Topjets_minDR_genb_vs_minDR_genTop->Clone());

  TH2D * h2_bTagged_Wjets_mass_vs_minDR_genb = new TH2D(sampleKeyStringT+"_h2_bTagged_Wjets_mass_vs_minDR_genb", sampleKeyStringT+": bTagged W jets :  Wjet mass versus minDR(Wjet, genb); #DeltaR(Wjet, genb); M_{Wjet} (GeV)", 25, 0, 1.5, 40, 70, 110); h2_bTagged_Wjets_mass_vs_minDR_genbVec.push_back((TH2D*)h2_bTagged_Wjets_mass_vs_minDR_genb->Clone());

  TH1D * h1_bTagged_Wjets_mass = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_mass", sampleKeyStringT+": bTagged W jets : Wjet mass; M_{Wjet} (GeV)", 40, 70, 110); h1_bTagged_Wjets_massVec.push_back((TH1D*) h1_bTagged_Wjets_mass->Clone());
  TH1D * h1_nobTagged_Wjets_mass = new TH1D(sampleKeyStringT+"_h1_nobTagged_Wjets_mass", sampleKeyStringT+": nobTagged W jets : Wjet mass; M_{Wjet} (GeV)", 40, 70, 110); h1_nobTagged_Wjets_massVec.push_back((TH1D*) h1_nobTagged_Wjets_mass->Clone());

  TH1D * h1_bTagged_rndmComb_Wjets_pt = new TH1D(sampleKeyStringT+"_h1_bTagged_rndmComb_Wjets_pt", sampleKeyStringT+": bTagged W jets : rndm combination Wjet pt; P_{T}^{Wjet} (GeV)", 100, 0, 1000); h1_bTagged_rndmComb_Wjets_ptVec.push_back((TH1D*) h1_bTagged_rndmComb_Wjets_pt->Clone());
  TH1D * h1_bTagged_sameTop_Wjets_pt = new TH1D(sampleKeyStringT+"_h1_bTagged_sameTop_Wjets_pt", sampleKeyStringT+": bTagged W jets : same top combination Wjet pt; P_{T}^{Wjet} (GeV)", 100, 0, 1000); h1_bTagged_sameTop_Wjets_ptVec.push_back((TH1D*) h1_bTagged_sameTop_Wjets_pt->Clone());

  TH1D * h1_bTagged_Wjets_deltaR1b_genDaus = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_deltaR1b_genDaus", sampleKeyStringT+": bTagged W jets : #DeltaR(1b) gen top daus; #DeltaR(1b)", 25, 0, 5); h1_bTagged_Wjets_deltaR1b_genDausVec.push_back((TH1D*)h1_bTagged_Wjets_deltaR1b_genDaus->Clone());
  TH1D * h1_bTagged_Wjets_deltaR2b_genDaus = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_deltaR2b_genDaus", sampleKeyStringT+": bTagged W jets : #DeltaR(2b) gen top daus; #DeltaR(2b)", 25, 0, 5); h1_bTagged_Wjets_deltaR2b_genDausVec.push_back((TH1D*)h1_bTagged_Wjets_deltaR2b_genDaus->Clone());
  TH1D * h1_bTagged_Wjets_deltaR12_genDaus = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_deltaR12_genDaus", sampleKeyStringT+": bTagged W jets : #DeltaR(12) gen top daus; #DeltaR(12)", 25, 0, 5); h1_bTagged_Wjets_deltaR12_genDausVec.push_back((TH1D*)h1_bTagged_Wjets_deltaR12_genDaus->Clone());

  TH1D * h1_bTagged_Wjets_minDR_otherbJets = new TH1D(sampleKeyStringT+"_h1_bTagged_Wjets_minDR_otherbJets", sampleKeyStringT+": bTagged W jets : min#DeltaR from another b jet to it; #DeltaR(Wjet, bJet)", 25, 0, 1.5); h1_bTagged_Wjets_minDR_otherbJetsVec.push_back((TH1D*)h1_bTagged_Wjets_minDR_otherbJets->Clone());

  TH2D * h2_nobTagged_Wjets_minDR_genb_vs_minDR_genW = new TH2D(sampleKeyStringT+"_h2_nobTagged_Wjets_minDR_genb_vs_minDR_genW", sampleKeyStringT+": nobTagged W jets :  minDR(Wjet, genb) versus minDR(Wjet, genW); #DeltaR(Wjet, genW); #DeltaR(Wjet, genb)", 25, 0, 1.5, 25, 0, 1.5); h2_nobTagged_Wjets_minDR_genb_vs_minDR_genWVec.push_back((TH2D*)h2_nobTagged_Wjets_minDR_genb_vs_minDR_genW->Clone());
  TH2D * h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTop = new TH2D(sampleKeyStringT+"_h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTop", sampleKeyStringT+": nobTagged Top jets :  minDR(Topjet, genb) versus minDR(Topjet, genTop); #DeltaR(Topjet, genTop); #DeltaR(Topjet, genb)", 25, 0, 1.5, 25, 0, 1.5); h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTopVec.push_back((TH2D*)h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTop->Clone());

  TH1D * h1_bTagged_Topjets_mass = new TH1D(sampleKeyStringT+"_h1_bTagged_Topjets_mass", sampleKeyStringT+": bTagged Top jets : Topjet mass; M_{Topjet} (GeV)", 45, 130, 220); h1_bTagged_Topjets_massVec.push_back((TH1D*) h1_bTagged_Topjets_mass->Clone());
  TH1D * h1_nobTagged_Topjets_mass = new TH1D(sampleKeyStringT+"_h1_nobTagged_Topjets_mass", sampleKeyStringT+": nobTagged Top jets : Topjet mass; M_{Topjet} (GeV)", 45, 130, 220); h1_nobTagged_Topjets_massVec.push_back((TH1D*) h1_nobTagged_Topjets_mass->Clone());

  TH1D * h1_genLeptORprong_Pt = new TH1D(sampleKeyStringT+"_h1_genLeptORprong_Pt", sampleKeyStringT+": gen lepton or prong; Pt (GeV)", 100, 0, 500); h1_genLeptORprong_PtVec.push_back((TH1D*) h1_genLeptORprong_Pt->Clone());
  TH1D * h1_genLeptORprong_Eta = new TH1D(sampleKeyStringT+"_h1_genLeptORprong_Eta", sampleKeyStringT+": gen lepton or prong; #eta", 100, -5, 5); h1_genLeptORprong_EtaVec.push_back((TH1D*) h1_genLeptORprong_Eta->Clone());
  TH1D * h1_genLeptORprong_Phi = new TH1D(sampleKeyStringT+"_h1_genLeptORprong_Phi", sampleKeyStringT+": gen lepton or prong; #phi", 100, -3.2, 3.2); h1_genLeptORprong_PhiVec.push_back((TH1D*) h1_genLeptORprong_Phi->Clone());

  TH1D * h1_genLeptORprong_minDR_triplets = new TH1D(sampleKeyStringT+"_ h1_genLeptORprong_minDR_triplets", sampleKeyStringT+": triplets: min #DeltaR between gen lepton or prong and jets; #DeltaR", 100, 0, 5); h1_genLeptORprong_minDR_tripletsVec.push_back((TH1D*) h1_genLeptORprong_minDR_triplets->Clone());
  TH1D * h1_genLeptORprong_minDphi_triplets = new TH1D(sampleKeyStringT+"_ h1_genLeptORprong_minDphi_triplets", sampleKeyStringT+": triplets: min #Delta#phi between gen lepton or prong and jets; #Delta#phi", 100, 0, 5); h1_genLeptORprong_minDphi_tripletsVec.push_back((TH1D*) h1_genLeptORprong_minDphi_triplets->Clone());
  TH1D * h1_genLeptORprong_minDR_Rsys = new TH1D(sampleKeyStringT+"_ h1_genLeptORprong_minDR_Rsys", sampleKeyStringT+": Rsys: min #DeltaR between gen lepton or prong and jets; #DeltaR", 100, 0, 5); h1_genLeptORprong_minDR_RsysVec.push_back((TH1D*) h1_genLeptORprong_minDR_Rsys->Clone());
  TH1D * h1_genLeptORprong_minDphi_Rsys = new TH1D(sampleKeyStringT+"_ h1_genLeptORprong_minDphi_Rsys", sampleKeyStringT+": Rsys: min #Delta#phi between gen lepton or prong and jets; #Delta#phi", 100, 0, 5); h1_genLeptORprong_minDphi_RsysVec.push_back((TH1D*) h1_genLeptORprong_minDphi_Rsys->Clone());
  TH1D * h1_genLeptORprong_dPhi_met = new TH1D(sampleKeyStringT+"_ h1_genLeptORprong_dPhi_met", sampleKeyStringT+": min #Delta#phi between gen lepton or prong and met; #Delta#phi", 100, 0, 5); h1_genLeptORprong_dPhi_metVec.push_back((TH1D*) h1_genLeptORprong_dPhi_met->Clone());

  TH2D * h2_genLeptORprong_minDR_Rsys_vs_triplets = new TH2D(sampleKeyStringT+"_h2_genLeptORprong_minDR_Rsys_vs_triplets", sampleKeyStringT+": min #DeltaR Rsys versus triplets; triplets; Rsys", 100, 0, 5, 100, 0, 5); h2_genLeptORprong_minDR_Rsys_vs_tripletsVec.push_back((TH2D*)h2_genLeptORprong_minDR_Rsys_vs_triplets->Clone());
  TH2D * h2_genLeptORprong_minDphi_Rsys_vs_triplets = new TH2D(sampleKeyStringT+"_h2_genLeptORprong_minDphi_Rsys_vs_triplets", sampleKeyStringT+": min #Delta#phi Rsys versus triplets; triplets; Rsys", 100, 0, 5, 100, 0, 5); h2_genLeptORprong_minDphi_Rsys_vs_tripletsVec.push_back((TH2D*)h2_genLeptORprong_minDphi_Rsys_vs_triplets->Clone());

  TH1D * h1_minMTj_lepJet = new TH1D(sampleKeyStringT+"_h1_minMTj_lepJet", sampleKeyStringT+": min MT between jets and MET (jets other than triplet and b); MTj (GeV)", 100, 0, 300); h1_minMTj_lepJetVec.push_back((TH1D*) h1_minMTj_lepJet->Clone());
  TH2D * h2_minMTj_pickedIdx_vs_minDR_pickedIdx = new TH2D(sampleKeyStringT+"_h2_minMTj_pickedIdx_vs_minDR_pickedIdx", sampleKeyStringT+": pickedIdx from minMTj versus pickedIdx from minDR; idx_{minDR}; idx_{minMTj}", 20, 0, 20, 20, 0, 20); h2_minMTj_pickedIdx_vs_minDR_pickedIdxVec.push_back((TH2D*)h2_minMTj_pickedIdx_vs_minDR_pickedIdx->Clone());

  TH1D * h1_minDphi_met_lepJet = new TH1D(sampleKeyStringT+"_h1_minDphi_met_lepJet", sampleKeyStringT+": min #Delta#phi between jets and MET (jets other than triplet and b); #Delta#phi", 100, 0, 5); h1_minDphi_met_lepJetVec.push_back((TH1D*) h1_minDphi_met_lepJet->Clone());
  TH2D * h2_minDphi_met_pickedIdx_vs_minDR_pickedIdx = new TH2D(sampleKeyStringT+"_h2_minDphi_met_pickedIdx_vs_minDR_pickedIdx", sampleKeyStringT+": pickedIdx from minDphi_met versus pickedIdx from minDR; idx_{minDR}; idx_{minDphi_met}", 20, 0, 20, 20, 0, 20); h2_minDphi_met_pickedIdx_vs_minDR_pickedIdxVec.push_back((TH2D*)h2_minDphi_met_pickedIdx_vs_minDR_pickedIdx->Clone());

  TH1D * h1_relPt_genLeptORprongOVERjetPt_triplets = new TH1D(sampleKeyStringT+"_h1_relPt_genLeptORprongOVERjetPt_triplets", sampleKeyStringT+": Pt(genLeptORprong)/Pt(jet in triplets); Pt(genLept)/Pt(tripletJet)", 100, 0, 1.5); h1_relPt_genLeptORprongOVERjetPt_tripletsVec.push_back((TH1D*)h1_relPt_genLeptORprongOVERjetPt_triplets->Clone());
  TH1D * h1_relPt_genLeptORprongOVERjetPt_Rsys = new TH1D(sampleKeyStringT+"_h1_relPt_genLeptORprongOVERjetPt_Rsys", sampleKeyStringT+": Pt(genLeptORprong)/Pt(jet in Rsys); Pt(genLept)/Pt(RsysJet)", 100, 0, 1.5); h1_relPt_genLeptORprongOVERjetPt_RsysVec.push_back((TH1D*)h1_relPt_genLeptORprongOVERjetPt_Rsys->Clone());

  TH2D * h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPt = new TH2D(sampleKeyStringT+"_h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPt", sampleKeyStringT+": Pt(genLeptORprong)/Pt(jet in triplets) versus genLeptORprong Pt; Pt(genLept) (GeV);  Pt(genLept)/Pt(tripletJet)", 100, 0, 500, 100, 0, 1.5); h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPtVec.push_back((TH2D*)h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPt->Clone()); 
  TH2D * h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPt = new TH2D(sampleKeyStringT+"_h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPt", sampleKeyStringT+": Pt(genLeptORprong)/Pt(jet in Rsys) versus genLeptORprong Pt; Pt(genLept) (GeV);  Pt(genLept)/Pt(RsysJet)", 100, 0, 500, 100, 0, 1.5); h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPtVec.push_back((TH2D*)h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPt->Clone());

  TH1D * h1_MT2_nTopsEQ2 = new TH1D(sampleKeyStringT+"_h1_MT2_nTopsEQ2", sampleKeyStringT+": MT2 for nTops = 2; M_{T2} (GeV)", 100, 0, 1000); h1_MT2_nTopsEQ2Vec.push_back((TH1D*)h1_MT2_nTopsEQ2->Clone());
  TH1D * h1_MT2_TopAndbLept_nTopsEQ2 = new TH1D(sampleKeyStringT+"_h1_MT2_TopAndbLept_nTopsEQ2", sampleKeyStringT+": MT2_TopAndbLept for nTops = 2; M_{T2} (GeV)", 100, 0, 1000); h1_MT2_TopAndbLept_nTopsEQ2Vec.push_back((TH1D*)h1_MT2_TopAndbLept_nTopsEQ2->Clone());
  TH1D * h1_minDphiLeptMET_nTopsEQ2 = new TH1D(sampleKeyStringT+"_h1_minDphiLeptMET_nTopsEQ2", sampleKeyStringT+": min#Delta#phi(lept, met) for nTops = 2; min#Delta#phi", 100, 0, 5); h1_minDphiLeptMET_nTopsEQ2Vec.push_back((TH1D*)h1_minDphiLeptMET_nTopsEQ2->Clone());
  TH2D * h2_mTW1_vs_mTtop1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTW1_vs_mTtop1_nTopsEQ2", sampleKeyStringT+": mTW1 versus mTtop1 for nTops = 2; M_{Ttop1} (GeV); M_{W1} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTW1_vs_mTtop1_nTopsEQ2Vec.push_back((TH2D*)h2_mTW1_vs_mTtop1_nTopsEQ2->Clone());
  TH2D * h2_mTb1_vs_mTtop1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb1_vs_mTtop1_nTopsEQ2", sampleKeyStringT+": mTb1 versus mTtop1 for nTops = 2; M_{Ttop1} (GeV); M_{b1} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb1_vs_mTtop1_nTopsEQ2Vec.push_back((TH2D*)h2_mTb1_vs_mTtop1_nTopsEQ2->Clone());
  TH2D * h2_mTb1_vs_mTW1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb1_vs_mTW1_nTopsEQ2", sampleKeyStringT+": mTb1 versus mTW1 for nTops = 2; M_{TW1} (GeV); M_{b1} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb1_vs_mTW1_nTopsEQ2Vec.push_back((TH2D*)h2_mTb1_vs_mTW1_nTopsEQ2->Clone());

  TH2D * h2_mTW2_vs_mTtop2_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTW2_vs_mTtop2_nTopsEQ2", sampleKeyStringT+": mTW2 versus mTtop2 for nTops = 2; M_{Ttop2} (GeV); M_{W2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTW2_vs_mTtop2_nTopsEQ2Vec.push_back((TH2D*)h2_mTW2_vs_mTtop2_nTopsEQ2->Clone());
  TH2D * h2_mTb2_vs_mTtop2_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb2_vs_mTtop2_nTopsEQ2", sampleKeyStringT+": mTb2 versus mTtop2 for nTops = 2; M_{Ttop2} (GeV); M_{b2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb2_vs_mTtop2_nTopsEQ2Vec.push_back((TH2D*)h2_mTb2_vs_mTtop2_nTopsEQ2->Clone());
  TH2D * h2_mTb2_vs_mTW2_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb2_vs_mTW2_nTopsEQ2", sampleKeyStringT+": mTb2 versus mTW2 for nTops = 2; M_{TW2} (GeV); M_{b2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb2_vs_mTW2_nTopsEQ2Vec.push_back((TH2D*)h2_mTb2_vs_mTW2_nTopsEQ2->Clone());

  TH2D * h2_mTtop2_vs_mTtop1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTtop2_vs_mTtop1_nTopsEQ2", sampleKeyStringT+": mTtop2 versus mTtop1 for nTops = 2; M_{Ttop1} (GeV); M_{top2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTtop2_vs_mTtop1_nTopsEQ2Vec.push_back((TH2D*)h2_mTtop2_vs_mTtop1_nTopsEQ2->Clone());
  TH2D * h2_mTW2_vs_mTtop1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTW2_vs_mTtop1_nTopsEQ2", sampleKeyStringT+": mTW2 versus mTtop1 for nTops = 2; M_{Ttop1} (GeV); M_{W2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTW2_vs_mTtop1_nTopsEQ2Vec.push_back((TH2D*)h2_mTW2_vs_mTtop1_nTopsEQ2->Clone());
  TH2D * h2_mTb2_vs_mTtop1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb2_vs_mTtop1_nTopsEQ2", sampleKeyStringT+": mTb2 versus mTtop1 for nTops = 2; M_{Ttop1} (GeV); M_{b2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb2_vs_mTtop1_nTopsEQ2Vec.push_back((TH2D*)h2_mTb2_vs_mTtop1_nTopsEQ2->Clone());

  TH2D * h2_mTtop2_vs_mTW1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTtop2_vs_mTW1_nTopsEQ2", sampleKeyStringT+": mTtop2 versus mTW1 for nTops = 2; M_{TW1} (GeV); M_{top2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTtop2_vs_mTW1_nTopsEQ2Vec.push_back((TH2D*)h2_mTtop2_vs_mTW1_nTopsEQ2->Clone());
  TH2D * h2_mTW2_vs_mTW1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTW2_vs_mTW1_nTopsEQ2", sampleKeyStringT+": mTW2 versus mTW1 for nTops = 2; M_{TW1} (GeV); M_{W2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTW2_vs_mTW1_nTopsEQ2Vec.push_back((TH2D*)h2_mTW2_vs_mTW1_nTopsEQ2->Clone());
  TH2D * h2_mTb2_vs_mTW1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb2_vs_mTW1_nTopsEQ2", sampleKeyStringT+": mTb2 versus mTW1 for nTops = 2; M_{TW1} (GeV); M_{b2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb2_vs_mTW1_nTopsEQ2Vec.push_back((TH2D*)h2_mTb2_vs_mTW1_nTopsEQ2->Clone());

  TH2D * h2_mTtop2_vs_mTb1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTtop2_vs_mTb1_nTopsEQ2", sampleKeyStringT+": mTtop2 versus mTb1 for nTops = 2; M_{Tb1} (GeV); M_{top2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTtop2_vs_mTb1_nTopsEQ2Vec.push_back((TH2D*)h2_mTtop2_vs_mTb1_nTopsEQ2->Clone());
  TH2D * h2_mTW2_vs_mTb1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTW2_vs_mTb1_nTopsEQ2", sampleKeyStringT+": mTW2 versus mTb1 for nTops = 2; M_{Tb1} (GeV); M_{W2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTW2_vs_mTb1_nTopsEQ2Vec.push_back((TH2D*)h2_mTW2_vs_mTb1_nTopsEQ2->Clone());
  TH2D * h2_mTb2_vs_mTb1_nTopsEQ2 = new TH2D(sampleKeyStringT+"_h2_mTb2_vs_mTb1_nTopsEQ2", sampleKeyStringT+": mTb2 versus mTb1 for nTops = 2; M_{Tb1} (GeV); M_{b2} (GeV)", 100, 0, 1000, 100, 0, 1000); h2_mTb2_vs_mTb1_nTopsEQ2Vec.push_back((TH2D*)h2_mTb2_vs_mTb1_nTopsEQ2->Clone());

  for(unsigned int it=0; it<nTopCandToPlot; it++){
     std::vector<TH1D*> h1_PS_topCand_MVec, h1_PS_topCand_PtVec, h1_PS_topCand_EtaVec, h1_PS_topCand_nJetsVec;
     sprintf(names, "%s_h1_topCand_M_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: topCand mass for the %dth nTop; M_{top} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_topCand_M = new TH1D(names, dispt, 100, 0, 300);
     if( h1_topCand_MVec.size() < it+1 ){ h1_PS_topCand_MVec.push_back((TH1D*)h1_topCand_M->Clone()); h1_topCand_MVec.push_back(h1_PS_topCand_MVec); }
     else{ h1_topCand_MVec[it].push_back((TH1D*)h1_topCand_M->Clone()); }

     sprintf(names, "%s_h1_topCand_Pt_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: topCand Pt for the %dth nTop; Pt_{top} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_topCand_Pt = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_topCand_PtVec.size() < it+1 ){ h1_PS_topCand_PtVec.push_back((TH1D*)h1_topCand_Pt->Clone()); h1_topCand_PtVec.push_back(h1_PS_topCand_PtVec); }
     else{ h1_topCand_PtVec[it].push_back((TH1D*)h1_topCand_Pt->Clone()); }

     sprintf(names, "%s_h1_topCand_Eta_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: topCand #eta for the %dth nTop; Eta_{top}", sampleKeyStringT.Data(), it);
     TH1D * h1_topCand_Eta = new TH1D(names, dispt, 100, -5, -5);
     if( h1_topCand_EtaVec.size() < it+1 ){ h1_PS_topCand_EtaVec.push_back((TH1D*)h1_topCand_Eta->Clone()); h1_topCand_EtaVec.push_back(h1_PS_topCand_EtaVec); }
     else{ h1_topCand_EtaVec[it].push_back((TH1D*)h1_topCand_Eta->Clone()); }

     sprintf(names, "%s_h1_topCand_nJets_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: nJets in topCands for the %dth nTop; nJets_{top}", sampleKeyStringT.Data(), it);
     TH1D * h1_topCand_nJets = new TH1D(names, dispt, 5, 0, 5);
     if( h1_topCand_nJetsVec.size() < it+1 ){ h1_PS_topCand_nJetsVec.push_back((TH1D*)h1_topCand_nJets->Clone()); h1_topCand_nJetsVec.push_back(h1_PS_topCand_nJetsVec); }
     else{ h1_topCand_nJetsVec[it].push_back((TH1D*)h1_topCand_nJets->Clone()); }

     std::vector<TH1D*> h1_PS_W_MVec, h1_PS_W_PtVec, h1_PS_W_EtaVec;
     sprintf(names, "%s_h1_W_M_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: W mass for the %dth nTop; M_{W} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_W_M = new TH1D(names, dispt, 100, 0, 300);
     if( h1_W_MVec.size() < it+1 ){ h1_PS_W_MVec.push_back((TH1D*)h1_W_M->Clone()); h1_W_MVec.push_back(h1_PS_W_MVec); }
     else{ h1_W_MVec[it].push_back((TH1D*)h1_W_M->Clone()); }

     sprintf(names, "%s_h1_W_Pt_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: W Pt for the %dth nTop; Pt_{W} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_W_Pt = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_W_PtVec.size() < it+1 ){ h1_PS_W_PtVec.push_back((TH1D*)h1_W_Pt->Clone()); h1_W_PtVec.push_back(h1_PS_W_PtVec); }
     else{ h1_W_PtVec[it].push_back((TH1D*)h1_W_Pt->Clone()); }

     sprintf(names, "%s_h1_W_Eta_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: W #eta for the %dth nTop; Eta_{W}", sampleKeyStringT.Data(), it);
     TH1D * h1_W_Eta = new TH1D(names, dispt, 100, -5, -5);
     if( h1_W_EtaVec.size() < it+1 ){ h1_PS_W_EtaVec.push_back((TH1D*)h1_W_Eta->Clone()); h1_W_EtaVec.push_back(h1_PS_W_EtaVec); }
     else{ h1_W_EtaVec[it].push_back((TH1D*)h1_W_Eta->Clone()); }

     std::vector<TH1D*> h1_PS_W_dausDRVec;
     std::vector<TH2D*> h2_PS_W_dausDR_versus_W_MVec;

     sprintf(names, "%s_h1_W_dausDR_nTopIdx_%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR between W daus for the %dth nTop; #DeltaR", sampleKeyStringT.Data(), it);
     TH1D * h1_W_dausDR = new TH1D(names, dispt, 100, 0, 5);
     if( h1_W_dausDRVec.size() < it+1 ){ h1_PS_W_dausDRVec.push_back((TH1D*) h1_W_dausDR->Clone()); h1_W_dausDRVec.push_back(h1_PS_W_dausDRVec); }
     else{ h1_W_dausDRVec[it].push_back((TH1D*) h1_W_dausDR->Clone()); }

     sprintf(names, "%s_h1_W_dausDR_versus_W_M_nTopIdx_%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR between W daus versus W mass for the %dth nTop; M_{W} (GeV); #DeltaR", sampleKeyStringT.Data(), it);
     TH2D * h2_W_dausDR_versus_W_M = new TH2D(names, dispt, 100, 0, 300, 100, 0, 5);
     if( h2_W_dausDR_versus_W_MVec.size() < it+1 ){ h2_PS_W_dausDR_versus_W_MVec.push_back((TH2D*) h2_W_dausDR_versus_W_M->Clone()); h2_W_dausDR_versus_W_MVec.push_back(h2_PS_W_dausDR_versus_W_MVec); }
     else{ h2_W_dausDR_versus_W_MVec[it].push_back((TH2D*) h2_W_dausDR_versus_W_M->Clone()); }

     std::vector<TH1D*> h1_PS_b_MVec, h1_PS_b_PtVec, h1_PS_b_EtaVec;
     sprintf(names, "%s_h1_b_M_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: b mass for the %dth nTop; M_{b} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_b_M = new TH1D(names, dispt, 100, 0, 300);
     if( h1_b_MVec.size() < it+1 ){ h1_PS_b_MVec.push_back((TH1D*)h1_b_M->Clone()); h1_b_MVec.push_back(h1_PS_b_MVec); }
     else{ h1_b_MVec[it].push_back((TH1D*)h1_b_M->Clone()); }

     sprintf(names, "%s_h1_b_Pt_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: b Pt for the %dth nTop; Pt_{b} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_b_Pt = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_b_PtVec.size() < it+1 ){ h1_PS_b_PtVec.push_back((TH1D*)h1_b_Pt->Clone()); h1_b_PtVec.push_back(h1_PS_b_PtVec); }
     else{ h1_b_PtVec[it].push_back((TH1D*)h1_b_Pt->Clone()); }

     sprintf(names, "%s_h1_b_Eta_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: b #eta for the %dth nTop; Eta_{b}", sampleKeyStringT.Data(), it);
     TH1D * h1_b_Eta = new TH1D(names, dispt, 100, -5, -5);
     if( h1_b_EtaVec.size() < it+1 ){ h1_PS_b_EtaVec.push_back((TH1D*)h1_b_Eta->Clone()); h1_b_EtaVec.push_back(h1_PS_b_EtaVec); }
     else{ h1_b_EtaVec[it].push_back((TH1D*)h1_b_Eta->Clone()); }

     std::vector<TH1D*> h1_PS_ratio_mW_over_mTopVec;
     sprintf(names, "%s_h1_ratio_mW_over_mTop_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mW/mTop for the %dth nTop; M_{W}/M_{top}", sampleKeyStringT.Data(), it);
     TH1D * h1_ratio_mW_over_mTop = new TH1D(names, dispt, 100, 0, 1.5);
     if( h1_ratio_mW_over_mTopVec.size() < it+1 ){ h1_PS_ratio_mW_over_mTopVec.push_back((TH1D*) h1_ratio_mW_over_mTop->Clone()); h1_ratio_mW_over_mTopVec.push_back(h1_PS_ratio_mW_over_mTopVec); }
     else{ h1_ratio_mW_over_mTopVec[it].push_back((TH1D*)h1_ratio_mW_over_mTop->Clone()); }

     std::vector<TH2D*> h2_PS_mW_versus_mTopVec;
     sprintf(names, "%s_h2_mW_versus_mTop_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mW versus mTop for the %dth nTop; M_{top} (GeV); M_{W} (GeV)", sampleKeyStringT.Data(), it);
     TH2D * h2_mW_versus_mTop = new TH2D(names, dispt, 100, 0, 300, 100, 0, 300);
     if( h2_mW_versus_mTopVec.size() < it+1 ){ h2_PS_mW_versus_mTopVec.push_back((TH2D*) h2_mW_versus_mTop->Clone()); h2_mW_versus_mTopVec.push_back(h2_PS_mW_versus_mTopVec); }
     else{ h2_mW_versus_mTopVec[it].push_back((TH2D*)h2_mW_versus_mTop->Clone()); }

     std::vector<TH2D*> h2_PS_ratio_mW_over_mTop_versus_mTopVec;
     sprintf(names, "%s_h2_ratio_mW_over_mTop_versus_mTop_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mW/mTop versus mTop for the %dth nTop; M_{top} (GeV); M_{W}/M_{top}", sampleKeyStringT.Data(), it);
     TH2D * h2_ratio_mW_over_mTop_versus_mTop = new TH2D(names, dispt, 100, 0, 300, 100, 0, 1.5);
     if( h2_ratio_mW_over_mTop_versus_mTopVec.size() < it+1 ){ h2_PS_ratio_mW_over_mTop_versus_mTopVec.push_back((TH2D*) h2_ratio_mW_over_mTop_versus_mTop->Clone()); h2_ratio_mW_over_mTop_versus_mTopVec.push_back(h2_PS_ratio_mW_over_mTop_versus_mTopVec); }
     else{ h2_ratio_mW_over_mTop_versus_mTopVec[it].push_back((TH2D*)h2_ratio_mW_over_mTop_versus_mTop->Clone()); }

     std::vector<TH1D*> h1_PS_mTtopVec;
     sprintf(names, "%s_h1_mTtop_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mTtop for the %dth nTop; M_{Ttop} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_mTtop = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_mTtopVec.size() < it+1 ){ h1_PS_mTtopVec.push_back((TH1D*) h1_mTtop->Clone()); h1_mTtopVec.push_back(h1_PS_mTtopVec); }
     else{ h1_mTtopVec[it].push_back((TH1D*)h1_mTtop->Clone()); }

     std::vector<TH1D*> h1_PS_mTWVec;
     sprintf(names, "%s_h1_mTW_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mTW for the %dth nTop; M_{TW} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_mTW = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_mTWVec.size() < it+1 ){ h1_PS_mTWVec.push_back((TH1D*) h1_mTW->Clone()); h1_mTWVec.push_back(h1_PS_mTWVec); }
     else{ h1_mTWVec[it].push_back((TH1D*)h1_mTW->Clone()); }

     std::vector<TH1D*> h1_PS_mTbVec;
     sprintf(names, "%s_h1_mTb_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: mTb for the %dth nTop; M_{Tb} (GeV)", sampleKeyStringT.Data(), it);
     TH1D * h1_mTb = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_mTbVec.size() < it+1 ){ h1_PS_mTbVec.push_back((TH1D*) h1_mTb->Clone()); h1_mTbVec.push_back(h1_PS_mTbVec); }
     else{ h1_mTbVec[it].push_back((TH1D*)h1_mTb->Clone()); }

     std::vector<TH1D*> h1_PS_deltaR12Vec, h1_PS_deltaR1bVec, h1_PS_deltaR2bVec;
     std::vector<TH2D*> h2_PS_deltaR12_vs_deltaR1bVec, h2_PS_deltaR12_vs_deltaR2bVec, h2_PS_deltaR1b_vs_deltaR2bVec;

     sprintf(names, "%s_h1_deltaR12_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{12} for the %dth nTop; #DeltaR_{12}", sampleKeyStringT.Data(), it);
     TH1D * h1_deltaR12 = new TH1D(names, dispt, 100, 0, 5.0);
     if( h1_deltaR12Vec.size() < it+1 ){ h1_PS_deltaR12Vec.push_back((TH1D*)h1_deltaR12->Clone()); h1_deltaR12Vec.push_back(h1_PS_deltaR12Vec); }
     else{ h1_deltaR12Vec[it].push_back((TH1D*)h1_deltaR12->Clone()); }

     sprintf(names, "%s_h1_deltaR1b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{1b} for the %dth nTop; #DeltaR_{1b}", sampleKeyStringT.Data(), it);
     TH1D * h1_deltaR1b = new TH1D(names, dispt, 100, 0, 5.0);
     if( h1_deltaR1bVec.size() < it+1 ){ h1_PS_deltaR1bVec.push_back((TH1D*)h1_deltaR1b->Clone()); h1_deltaR1bVec.push_back(h1_PS_deltaR1bVec); }
     else{ h1_deltaR1bVec[it].push_back((TH1D*)h1_deltaR1b->Clone()); }

     sprintf(names, "%s_h1_deltaR2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{2b} for the %dth nTop; #DeltaR_{2b}", sampleKeyStringT.Data(), it);
     TH1D * h1_deltaR2b = new TH1D(names, dispt, 100, 0, 5.0);
     if( h1_deltaR2bVec.size() < it+1 ){ h1_PS_deltaR2bVec.push_back((TH1D*)h1_deltaR2b->Clone()); h1_deltaR2bVec.push_back(h1_PS_deltaR2bVec); }
     else{ h1_deltaR2bVec[it].push_back((TH1D*)h1_deltaR2b->Clone()); }

     sprintf(names, "%s_h2_deltaR12_vs_deltaR1b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{12} vs #DeltaR_{1b} for the %dth nTop; #DeltaR_{1b}; #DeltaR_{12}", sampleKeyStringT.Data(), it);
     TH2D * h2_deltaR12_vs_deltaR1b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 5.0);
     if( h2_deltaR12_vs_deltaR1bVec.size() < it+1 ){ h2_PS_deltaR12_vs_deltaR1bVec.push_back((TH2D*)h2_deltaR12_vs_deltaR1b->Clone()); h2_deltaR12_vs_deltaR1bVec.push_back(h2_PS_deltaR12_vs_deltaR1bVec); }
     else{ h2_deltaR12_vs_deltaR1bVec[it].push_back((TH2D*)h2_deltaR12_vs_deltaR1b->Clone()); }

     sprintf(names, "%s_h2_deltaR12_vs_deltaR2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{12} vs #DeltaR_{2b} for the %dth nTop; #DeltaR_{2b}; #DeltaR_{12}", sampleKeyStringT.Data(), it);
     TH2D * h2_deltaR12_vs_deltaR2b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 5.0);
     if( h2_deltaR12_vs_deltaR2bVec.size() < it+1 ){ h2_PS_deltaR12_vs_deltaR2bVec.push_back((TH2D*)h2_deltaR12_vs_deltaR2b->Clone()); h2_deltaR12_vs_deltaR2bVec.push_back(h2_PS_deltaR12_vs_deltaR2bVec); }
     else{ h2_deltaR12_vs_deltaR2bVec[it].push_back((TH2D*)h2_deltaR12_vs_deltaR2b->Clone()); }

     sprintf(names, "%s_h2_deltaR1b_vs_deltaR2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: #DeltaR_{1b} vs #DeltaR_{2b} for the %dth nTop; #DeltaR_{2b}; #DeltaR_{1b}", sampleKeyStringT.Data(), it);
     TH2D * h2_deltaR1b_vs_deltaR2b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 5.0);
     if( h2_deltaR1b_vs_deltaR2bVec.size() < it+1 ){ h2_PS_deltaR1b_vs_deltaR2bVec.push_back((TH2D*)h2_deltaR1b_vs_deltaR2b->Clone()); h2_deltaR1b_vs_deltaR2bVec.push_back(h2_PS_deltaR1b_vs_deltaR2bVec); }
     else{ h2_deltaR1b_vs_deltaR2bVec[it].push_back((TH2D*)h2_deltaR1b_vs_deltaR2b->Clone()); }

     std::vector<TH1D*> h1_PS_relPt12Vec, h1_PS_relPt1bVec, h1_PS_relPt2bVec;
     std::vector<TH2D*> h2_PS_relPt12_vs_relPt1bVec, h2_PS_relPt12_vs_relPt2bVec, h2_PS_relPt1b_vs_relPt2bVec;

     sprintf(names, "%s_h1_relPt12_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{12} for the %dth nTop; relPt_{12}", sampleKeyStringT.Data(), it);
     TH1D * h1_relPt12 = new TH1D(names, dispt, 100, 0, 1.5);
     if( h1_relPt12Vec.size() < it+1 ){ h1_PS_relPt12Vec.push_back((TH1D*)h1_relPt12->Clone()); h1_relPt12Vec.push_back(h1_PS_relPt12Vec); }
     else{ h1_relPt12Vec[it].push_back((TH1D*)h1_relPt12->Clone()); }

     sprintf(names, "%s_h1_relPt1b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{1b} for the %dth nTop; relPt_{1b}", sampleKeyStringT.Data(), it);
     TH1D * h1_relPt1b = new TH1D(names, dispt, 100, 0, 5.0);
     if( h1_relPt1bVec.size() < it+1 ){ h1_PS_relPt1bVec.push_back((TH1D*)h1_relPt1b->Clone()); h1_relPt1bVec.push_back(h1_PS_relPt1bVec); }
     else{ h1_relPt1bVec[it].push_back((TH1D*)h1_relPt1b->Clone()); }

     sprintf(names, "%s_h1_relPt2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{2b} for the %dth nTop; relPt_{2b}", sampleKeyStringT.Data(), it);
     TH1D * h1_relPt2b = new TH1D(names, dispt, 100, 0, 5.0);
     if( h1_relPt2bVec.size() < it+1 ){ h1_PS_relPt2bVec.push_back((TH1D*)h1_relPt2b->Clone()); h1_relPt2bVec.push_back(h1_PS_relPt2bVec); }
     else{ h1_relPt2bVec[it].push_back((TH1D*)h1_relPt2b->Clone()); }

     sprintf(names, "%s_h2_relPt12_vs_relPt1b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{12} vs relPt_{1b} for the %dth nTop; relPt_{1b}; relPt_{12}", sampleKeyStringT.Data(), it);
     TH2D * h2_relPt12_vs_relPt1b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 1.5);
     if( h2_relPt12_vs_relPt1bVec.size() < it+1 ){ h2_PS_relPt12_vs_relPt1bVec.push_back((TH2D*)h2_relPt12_vs_relPt1b->Clone()); h2_relPt12_vs_relPt1bVec.push_back(h2_PS_relPt12_vs_relPt1bVec); }
     else{ h2_relPt12_vs_relPt1bVec[it].push_back((TH2D*)h2_relPt12_vs_relPt1b->Clone()); }

     sprintf(names, "%s_h2_relPt12_vs_relPt2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{12} vs relPt_{2b} for the %dth nTop; relPt_{2b}; relPt_{12}", sampleKeyStringT.Data(), it);
     TH2D * h2_relPt12_vs_relPt2b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 1.5);
     if( h2_relPt12_vs_relPt2bVec.size() < it+1 ){ h2_PS_relPt12_vs_relPt2bVec.push_back((TH2D*)h2_relPt12_vs_relPt2b->Clone()); h2_relPt12_vs_relPt2bVec.push_back(h2_PS_relPt12_vs_relPt2bVec); }
     else{ h2_relPt12_vs_relPt2bVec[it].push_back((TH2D*)h2_relPt12_vs_relPt2b->Clone()); }

     sprintf(names, "%s_h2_relPt1b_vs_relPt2b_nTopIdx%d", sampleKeyStringT.Data(), it);
     sprintf(dispt, "%s: relPt_{1b} vs relPt_{2b} for the %dth nTop; relPt_{2b}; relPt_{1b}", sampleKeyStringT.Data(), it);
     TH2D * h2_relPt1b_vs_relPt2b = new TH2D(names, dispt, 100, 0, 5.0, 100, 0, 5.0);
     if( h2_relPt1b_vs_relPt2bVec.size() < it+1 ){ h2_PS_relPt1b_vs_relPt2bVec.push_back((TH2D*)h2_relPt1b_vs_relPt2b->Clone()); h2_relPt1b_vs_relPt2bVec.push_back(h2_PS_relPt1b_vs_relPt2bVec); }
     else{ h2_relPt1b_vs_relPt2bVec[it].push_back((TH2D*)h2_relPt1b_vs_relPt2b->Clone()); }

  }

  bool isDeclared = false;
  if( std::find(declaredSampleStrVec.begin(), declaredSampleStrVec.end(), sampleKeyStringT.Data()) != declaredSampleStrVec.end() ) isDeclared = true;
  else{ declaredSampleStrVec.push_back(sampleKeyString); }

  for(unsigned int iSR=0; iSR<nSR && !isDeclared; iSR++){
     std::vector<TH1D*> h1_PS_nJetsVec, h1_PS_metVec, h1_PS_MT2Vec, h1_PS_mTcombVec, h1_PS_HTVec, h1_PS_nJetsRsysVec;
     sprintf(names, "%s_h1_nJets_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: nJets for nbJets%s and nTops%s; nJets", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_nJets = new TH1D(names, dispt, 20, 0, 20); 
     if( h1_nJetsVec.size() < iSR+1 ){ h1_PS_nJetsVec.push_back((TH1D*)h1_nJets->Clone()); h1_nJetsVec.push_back(h1_PS_nJetsVec); }
     else{ h1_nJetsVec[iSR].push_back((TH1D*)h1_nJets->Clone()); }

     sprintf(names, "%s_h1_met_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: met for nbJets%s and nTops%s; met (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_met = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_metVec.size() < iSR+1 ){ h1_PS_metVec.push_back((TH1D*)h1_met->Clone()); h1_metVec.push_back(h1_PS_metVec); }
     else{ h1_metVec[iSR].push_back((TH1D*)h1_met->Clone()); }

     sprintf(names, "%s_h1_MT2_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: MT2 for nbJets%s and nTops%s; MT2 (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_MT2 = new TH1D(names, dispt, 100, 0, 1000);
     if( h1_MT2Vec.size() < iSR+1 ){ h1_PS_MT2Vec.push_back((TH1D*)h1_MT2->Clone()); h1_MT2Vec.push_back(h1_PS_MT2Vec); }
     else{ h1_MT2Vec[iSR].push_back((TH1D*)h1_MT2->Clone()); }

     sprintf(names, "%s_h1_mTcomb_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: mTcomb for nbJets%s and nTops%s; MTcomb (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_mTcomb = new TH1D(names, dispt, 150, 0, 1500);
     if( h1_mTcombVec.size() < iSR+1 ){ h1_PS_mTcombVec.push_back((TH1D*)h1_mTcomb->Clone()); h1_mTcombVec.push_back(h1_PS_mTcombVec); }
     else{ h1_mTcombVec[iSR].push_back((TH1D*)h1_mTcomb->Clone()); }

     sprintf(names, "%s_h1_HT_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: HT for nbJets%s and nTops%s; HT (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_HT = new TH1D(names, dispt, 150, 0, 1500);
     if( h1_HTVec.size() < iSR+1 ){ h1_PS_HTVec.push_back((TH1D*)h1_HT->Clone()); h1_HTVec.push_back(h1_PS_HTVec); }
     else{ h1_HTVec[iSR].push_back((TH1D*)h1_HT->Clone()); }

     sprintf(names, "%s_h1_nJetsRsys_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: remaining nJets for nbJets%s and nTops%s; nJetsRsys", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH1D * h1_nJetsRsys = new TH1D(names, dispt, 10, 0, 10); 
     if( h1_nJetsRsysVec.size() < iSR+1 ){ h1_PS_nJetsRsysVec.push_back((TH1D*)h1_nJetsRsys->Clone()); h1_nJetsRsysVec.push_back(h1_PS_nJetsRsysVec); }
     else{ h1_nJetsRsysVec[iSR].push_back((TH1D*)h1_nJetsRsys->Clone()); }

     std::vector<TH2D*> h2_PS_met_vs_nJetsVec;
     sprintf(names, "%s_h2_met_vs_nJets_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: met versus nJets for nbJets%s and nTops%s; nJets; met (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH2D * h2_met_vs_nJets = new TH2D(names, dispt, 10, 4, 14, 50, 200, 1000); 
     if( h2_met_vs_nJetsVec.size() < iSR+1 ) { h2_PS_met_vs_nJetsVec.push_back((TH2D*)h2_met_vs_nJets->Clone()); h2_met_vs_nJetsVec.push_back(h2_PS_met_vs_nJetsVec); }
     else{ h2_met_vs_nJetsVec[iSR].push_back((TH2D*)h2_met_vs_nJets->Clone()); }

     std::vector<TH2D*> h2_PS_MT2_vs_metVec;
     sprintf(names, "%s_h2_MT2_vs_met_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: MT2 versus met for nbJets%s and nTops%s; met (GeV); MT2 (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH2D * h2_MT2_vs_met = new TH2D(names, dispt, 50, 200, 1000, 100, 0, 1000);
     if( h2_MT2_vs_metVec.size() < iSR+1 ) { h2_PS_MT2_vs_metVec.push_back((TH2D*)h2_MT2_vs_met->Clone()); h2_MT2_vs_metVec.push_back(h2_PS_MT2_vs_metVec); }
     else{ h2_MT2_vs_metVec[iSR].push_back((TH2D*)h2_MT2_vs_met->Clone()); }

     std::vector<TH2D*> h2_PS_coarse_bin_MT2_vs_metVec;
     sprintf(names, "%s_h2_coarse_bin_MT2_vs_met_nbJets%s_nTops%s", sampleKeyStringT.Data(), keyStr_nbJets_SR[iSR].c_str(), keyStr_nTops_SR[iSR].c_str());
     sprintf(dispt, "%s: MT2 versus met for nbJets%s and nTops%s; met (GeV); MT2 (GeV)", sampleKeyStringT.Data(), disStr_nbJets_SR[iSR].c_str(), disStr_nTops_SR[iSR].c_str());
     TH2D * h2_coarse_bin_MT2_vs_met = new TH2D(names, dispt, 16, 200, 1000, 10, 0, 1000);
     if( h2_coarse_bin_MT2_vs_metVec.size() < iSR+1 ) { h2_PS_coarse_bin_MT2_vs_metVec.push_back((TH2D*)h2_coarse_bin_MT2_vs_met->Clone()); h2_coarse_bin_MT2_vs_metVec.push_back(h2_PS_coarse_bin_MT2_vs_metVec); }
     else{ h2_coarse_bin_MT2_vs_metVec[iSR].push_back((TH2D*)h2_coarse_bin_MT2_vs_met->Clone()); }
  }

  cnt_passSRmet_WeightedScaledMCVec.push_back(cnt_passSRmet_sumSM_WeightedScaledMCVec);
  cnt_passSRmet_WeightedErrorScaledMCVec.push_back(cnt_passSRmet_sumSM_WeightedErrorScaledMCVec);
  for(int iSR = 0; iSR < nSR; iSR++){
     for(int iSR_met =0; iSR_met < nSR_met; iSR_met++){
        cnt_passSRmet_WeightedScaledMCVec.back()[iSR][iSR_met] = 0;
        cnt_passSRmet_WeightedErrorScaledMCVec.back()[iSR][iSR_met] = 0;
     }
  }

  for(unsigned int ist=0; ist<subSampleKeysVec.size(); ist++){

     bool isData = false;

     std::string keyString = subSampleKeysVec[ist];

     double scaleMC = 1.0;
     for(int ib=0; ib<nMC; ib++){ if( mcStr[ib] == keyString ){ scaleMC = scalesVec[ib]; } }
     TString keyStringT(keyString);
     if( keyStringT.Contains("Data") ){ scaleMC = dataScale; isData = true; }

     if( tr ) delete tr;
     tr = new NTupleReader(treeVec[ist], AnaConsts::activatedBranchNames);
     tr->registerFunction(&passBaselineFunc);

     int entries = tr->getNEntries();
     std::cout<<"\n\n"<<keyString.c_str()<<"_entries : "<<entries<<"  scaleMC : "<<scaleMC<<std::endl;

     while(tr->getNextEvent()){

        if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (tr->getEvtNum()-1)%(entries/10) == 0 ) std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
/*
        if(tr->getEvtNum() == 1){
           tr->printTupleMembers();
           FILE * fout = fopen("NTupleTypes.txt", "w");
           tr->printTupleMembers(fout);
           fclose(fout);
        }
*/
        // Internal evtWeight in the sample: default is 1.0 execept for MC samples with intrinsic weight, e.g., QCD flat sample.
        double iniWeight = tr->getVar<double>("evtWeight");
        double puWeight = 1.0; // currently set to be 1.0
        
        double evtWeight = iniWeight * puWeight;

        if( !tr->getVar<bool>("passLeptVeto") ) continue;
        cnt_passLeptVeto_WeightedScaledMC += evtWeight * scaleMC; cnt_passLeptVeto_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passnJets") ) continue;
        cnt_passnJets_WeightedScaledMC += evtWeight * scaleMC; cnt_passnJets_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passdPhis") ) continue;
        cnt_passdPhis_WeightedScaledMC += evtWeight * scaleMC; cnt_passdPhis_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passBJets") ) continue;
        cnt_passBJets_WeightedScaledMC += evtWeight * scaleMC; cnt_passBJets_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passMET") ) continue;
        cnt_passMET_WeightedScaledMC += evtWeight * scaleMC; cnt_passMET_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passTagger") ) continue;
        cnt_passTagger_WeightedScaledMC += evtWeight * scaleMC; cnt_passTagger_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        if( !tr->getVar<bool>("passBaseline") ) continue;
        cnt_passBaseline_WeightedScaledMC += evtWeight * scaleMC; cnt_passBaseline_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC;

        const std::vector<TLorentzVector> & elesLVec = tr->getVec<TLorentzVector>("elesLVec"); const std::vector<double> & elesRelIso = tr->getVec<double>("elesRelIso");
        const std::vector<TLorentzVector> & muonsLVec = tr->getVec<TLorentzVector>("muonsLVec"); const std::vector<double> & muonsRelIso = tr->getVec<double>("muonsRelIso");
        const std::vector<TLorentzVector> & loose_isoTrksLVec = tr->getVec<TLorentzVector>("loose_isoTrksLVec"); const std::vector<double> & loose_isoTrks_iso = tr->getVec<double>("loose_isoTrks_iso");

        const std::vector<int> & W_emuVec = tr->getVec<int>("W_emuVec"), & W_tau_emuVec = tr->getVec<int>("W_tau_emuVec"), & W_tau_prongsVec = tr->getVec<int>("W_tau_prongsVec");
        const std::vector<int> & genDecayPdgIdVec = tr->getVec<int>("genDecayPdgIdVec"), & genDecayIdxVec = tr->getVec<int>("genDecayIdxVec"), & genDecayMomIdxVec = tr->getVec<int>("genDecayMomIdxVec");
        const std::vector<TLorentzVector> & genDecayLVec = tr->getVec<TLorentzVector>("genDecayLVec");
 
        const std::vector<TLorentzVector> & jetsLVec_forTagger = tr->getVec<TLorentzVector>("jetsLVec_forTagger");
        const std::vector<double> & recoJetsBtag_forTagger = tr->getVec<double>("recoJetsBtag_forTagger");
        TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr->getVar<double>("met"), 0, tr->getVar<double>("metphi"), 0);

        std::vector<int> genW_from_top_momIdxVec, genW_from_top_idxVec;
        for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
           const int pdgId = genDecayPdgIdVec.at(ig);
           if( std::abs(pdgId) == 24 ){
              int genW_mom_odx = genDecayMomIdxVec[ig];
              const int momIdx = find_idx(genW_mom_odx, genDecayIdxVec);
              if( std::abs(genDecayPdgIdVec[momIdx]) == 6 ){
                 genW_from_top_idxVec.push_back(ig); genW_from_top_momIdxVec.push_back(momIdx);
              }
           }
        }
        for(unsigned int iw=0; iw<genW_from_top_idxVec.size(); iw++){
           const int genW_idx = genW_from_top_idxVec.at(iw);
           const int genW_odx = genDecayIdxVec[genW_idx];
           int cntLept = 0;
           std::vector<int> genW_daus_from_top_idxVec;
           std::vector<int> genb_from_top_idxVec;
           for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
              int pdgId = genDecayPdgIdVec.at(ig);
              if( genDecayMomIdxVec[ig] == genW_odx ){
                 if( std::abs(pdgId) >= 11 ) cntLept ++;
                 genW_daus_from_top_idxVec.push_back(ig);
              }
              if( std::abs(pdgId) == 5 ){
                 int genb_mom_odx = genDecayMomIdxVec[ig];
                 const int mom_idx = find_idx(genb_mom_odx, genDecayIdxVec);
                 if( mom_idx == genW_from_top_momIdxVec[iw] ) genb_from_top_idxVec.push_back(ig); 
              }
           }
//           if( !cntLept ){
           if( cntLept ){
              if( genb_from_top_idxVec.size() > 1 ) std::cout<<"\nMore than 1 b from a top quark?!"<<std::endl;
              if( genW_daus_from_top_idxVec.size() >=2 && genb_from_top_idxVec.size() ==1 ){
                 TLorentzVector genbLVec = genDecayLVec.at(genb_from_top_idxVec.at(0));
                 TLorentzVector genWdau1LVec = genDecayLVec.at(genW_daus_from_top_idxVec.at(0));
                 TLorentzVector genWdau2LVec = genDecayLVec.at(genW_daus_from_top_idxVec.at(1));
                 double dR1b = genbLVec.DeltaR(genWdau1LVec);
                 double dR2b = genbLVec.DeltaR(genWdau2LVec);
                 double dR12 = genWdau1LVec.DeltaR(genWdau2LVec);

                 double dPhi1MET = metLVec.DeltaPhi(genWdau1LVec);
                 double dPhi2MET = metLVec.DeltaPhi(genWdau2LVec);

                 double mass1 = (genbLVec+genWdau1LVec).M();
                 double mass2 = (genbLVec+genWdau2LVec).M();
                 double corrdR1b = dR1b;
                 double corrdR2b = dR2b;
                 double corrdPhi1MET = dPhi1MET;
                 double corrdPhi2MET = dPhi2MET;

                 int pdgId1 = genDecayPdgIdVec.at(genW_daus_from_top_idxVec.at(0)), pdgId2 = genDecayPdgIdVec.at(genW_daus_from_top_idxVec.at(1));
                 if( cntLept && (std::abs(pdgId1) == 12 || std::abs(pdgId1) == 14 || std::abs(pdgId1) == 16) ){
                    mass1 = (genbLVec+genWdau2LVec).M();
                    mass2 = (genbLVec+genWdau1LVec).M();
                    corrdR1b = dR2b;
                    corrdR2b = dR1b;
                    corrdPhi1MET = dPhi2MET;
                    corrdPhi2MET = dPhi1MET;
                 }

                 if( !cntLept && genWdau1LVec.Pt() < genWdau2LVec.Pt() ){
                    mass1 = (genbLVec+genWdau2LVec).M();
                    mass2 = (genbLVec+genWdau1LVec).M();
                    corrdR1b = dR2b;
                    corrdR2b = dR1b;
                    corrdPhi1MET = dPhi2MET;
                    corrdPhi2MET = dPhi1MET;
                 }

                 h1_gen1b_massVec.back()->Fill(mass1, evtWeight * scaleMC);
                 h1_gen2b_massVec.back()->Fill(mass2, evtWeight * scaleMC);

                 h1_gen1b_deltaRVec.back()->Fill(corrdR1b, evtWeight * scaleMC);
                 h1_gen2b_deltaRVec.back()->Fill(corrdR2b, evtWeight * scaleMC);
                 h1_gen12_deltaRVec.back()->Fill(dR12, evtWeight * scaleMC);

                 h1_gen1MET_deltaPhiVec.back()->Fill(corrdPhi1MET, evtWeight * scaleMC);
                 h1_gen2MET_deltaPhiVec.back()->Fill(corrdPhi2MET, evtWeight * scaleMC);
              }
           }
        }

        int addbJets = 0;
        for(unsigned int ip=0; ip<type3Ptr->pickedTopCandSortedVec.size(); ip++){
           const int combIdx = type3Ptr->pickedTopCandSortedVec[ip];
           const std::vector<int> & topCombIdxVec = type3Ptr->finalCombfatJets[combIdx];
           int foundbJets = 0;
           for(unsigned int ic=0; ic<topCombIdxVec.size(); ic++){
              std::vector<TLorentzVector> dummyLVec; dummyLVec.push_back(jetsLVec_forTagger.at(topCombIdxVec.at(ic)));
              std::vector<double> dummyCSVS; dummyCSVS.push_back(recoJetsBtag_forTagger.at(topCombIdxVec.at(ic)));
              if( AnaFunctions::countCSVS(dummyLVec, dummyCSVS, AnaConsts::cutCSVS, AnaConsts::bTagArr) ) foundbJets++;
           }
           if( foundbJets >=2 ) std::cout<<"\nMore than 2 b jets in a top candidate?!"<<std::endl<<std::endl;
           if( !foundbJets ) addbJets++;
        }

        const int nbJets = tr->getVar<int>("cntCSVS"), nTops = tr->getVar<int>("nTopCandSortedCnt");
        const double MT2 = tr->getVar<double>("MT2"), mTcomb = tr->getVar<double>("mTcomb");

        double best_lept_brJet_MT = 9999., best_lept_brJet_MT2 = 9999.;
        double best_had_brJet_MT = 9999., best_had_brJet_MT2 = 9999.;
        double best_lept_brJet_mTcomb = 9999., best_had_brJet_mTcomb = 9999.;
        double best_bJet_MT = 9999.;
        double best_aft_PCA_had_brJet_MT = 9999., best_aft_PCA_top_MT = 9999.;

//        if( nTops == 1 && nbJets ==1 ){
//        if( nTops == 2 ){
        if( nTops ==1 || nTops == 2 ){
//        if( nTops == 1 ){

           h1_cutFlowVec.back()->Fill("original", evtWeight * scaleMC);

           const int combIdx = type3Ptr->ori_pickedTopCandSortedVec[0];
           const std::vector<int> & topCombIdxVec = type3Ptr->finalCombfatJets[combIdx];
           const std::vector<int> & passStatusVec = type3Ptr->finalCombfatJetsPassStatusVec[combIdx];
           const TLorentzVector topLVec = type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec);

           int passCnt = 0, pass1st =0, pass2nd =0, pass3rd =0;
           int pass1 =0, pass2 = 0, pass3 = 0, pass12 = 0, pass13 = 0, pass23 = 0;
           for(unsigned int ip=0; ip<passStatusVec.size(); ip++){
              if( passStatusVec[ip] ){
                 passCnt++; 
                 if( ip == 0 ) pass1st = 1;
                 if( ip == 1 ) pass2nd = 1;
                 if( ip == 2 ) pass3rd = 1;
              }
           }
           if( pass1st && !pass2nd && !pass3rd ) pass1 = 1;
           if( !pass1st && pass2nd && !pass3rd ) pass2 = 1;
           if( !pass1st && !pass2nd && pass3rd ) pass3 = 1;
           if( pass1st && pass2nd && !pass3rd ) pass12 = 1;
           if( pass1st && !pass2nd && pass3rd ) pass13 = 1;
           if( !pass1st && pass2nd && pass3rd ) pass23 = 1;

           if( topCombIdxVec.size() >=3 ){
              TLorentzVector reco_top_dau1LVec = jetsLVec_forTagger.at(topCombIdxVec[0]);
              TLorentzVector reco_top_dau2LVec = jetsLVec_forTagger.at(topCombIdxVec[1]);
              TLorentzVector reco_top_dau3LVec = jetsLVec_forTagger.at(topCombIdxVec[2]);
              double reco_top_m12 = (reco_top_dau1LVec + reco_top_dau2LVec).M();
              double reco_top_m13 = (reco_top_dau1LVec + reco_top_dau3LVec).M();
              double reco_top_m23 = (reco_top_dau2LVec + reco_top_dau3LVec).M();
              double reco_top_m123 = topLVec.M();
              h2_m23overm123vsarctanm13overm12Vec.back()->Fill(TMath::ATan(reco_top_m13/reco_top_m12), reco_top_m23/reco_top_m123, evtWeight * scaleMC);
           }
   
           std::vector<TLorentzVector> remainJetsLVec; std::vector<double> remainJetsCSVS;
           std::vector<TLorentzVector> remainbJetsLVec; std::vector<double> remainbJetsCSVS;
           std::vector<int> remainbJetsIdxVec;
           int isFaked_b = 0;
           double fake_maxCSVS_LT_Mpt = -1; int fake_pickedfakebJet_idx = -1;
           for(unsigned int ij=0; ij<jetsLVec_forTagger.size(); ij++){
              if( std::find(topCombIdxVec.begin(), topCombIdxVec.end(), ij) != topCombIdxVec.end() ) continue;
   
              if( jetsLVec_forTagger.at(ij).Pt() > AnaConsts::bTagArr[2] && std::abs(jetsLVec_forTagger.at(ij).Eta()) < AnaConsts::bTagArr[1] ){
                 if( recoJetsBtag_forTagger.at(ij) < AnaConsts::cutCSVS ){
                    if( fake_maxCSVS_LT_Mpt ==-1 || fake_maxCSVS_LT_Mpt < recoJetsBtag_forTagger.at(ij) ){
                       fake_pickedfakebJet_idx = ij; fake_maxCSVS_LT_Mpt = recoJetsBtag_forTagger.at(ij);
                    }
                 }
              }
              std::vector<TLorentzVector> dummyLVec; dummyLVec.push_back(jetsLVec_forTagger.at(ij));
              std::vector<double> dummyCSVS; dummyCSVS.push_back(recoJetsBtag_forTagger.at(ij));
              if( AnaFunctions::countCSVS(dummyLVec, dummyCSVS, AnaConsts::cutCSVS, AnaConsts::bTagArr) ){
                 remainbJetsLVec.push_back(jetsLVec_forTagger.at(ij));
                 remainbJetsCSVS.push_back(recoJetsBtag_forTagger.at(ij));
                 remainbJetsIdxVec.push_back(ij);
              }else{
                 remainJetsLVec.push_back(jetsLVec_forTagger.at(ij));
                 remainJetsCSVS.push_back(recoJetsBtag_forTagger.at(ij));
              }
           }
           if( fake_pickedfakebJet_idx == -1 && remainbJetsLVec.empty() ) std::cout<<"WARNING ... cannot fake b jet?!"<<std::endl;

// no b jets in the remaining system, then fake one
           if( remainbJetsLVec.empty() ){
              isFaked_b ++;
              remainJetsLVec.clear(); remainJetsCSVS.clear();
              for(unsigned int ij=0; ij<jetsLVec_forTagger.size(); ij++){
                 if( std::find(topCombIdxVec.begin(), topCombIdxVec.end(), ij) != topCombIdxVec.end() ) continue;
                 if( fake_pickedfakebJet_idx == ij ){
		   remainbJetsLVec.push_back(jetsLVec_forTagger.at(ij)); 
		   remainbJetsCSVS.push_back(recoJetsBtag_forTagger.at(ij)); 
		 }
                 else{
                    remainJetsLVec.push_back(jetsLVec_forTagger.at(ij));
                    remainJetsCSVS.push_back(recoJetsBtag_forTagger.at(ij));
                 }
              }
           }
           h1_cutFlowVec.back()->Fill("isFaked_b", (isFaked_b !=0) * evtWeight * scaleMC);
	   /************************************************************************************/
	   /************************************************************************************/


           std::vector<int> cachedgenWIdxVec;
           for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
              const int pdgId = genDecayPdgIdVec.at(ig);
              if( std::abs(pdgId) == 24 ) cachedgenWIdxVec.push_back(ig);
           }

	   /************************************************************************************/
	   /************************************************************************************/
           std::vector<std::vector<TLorentzVector> > topDausLVec;
           std::vector<std::vector<int> > topDausIdxVec;
           std::vector<TLorentzVector> genTopLVec;

	   // 0: hadronic  1: leptonic
           std::vector<int> topDecayTypeVec;
           bool badGenInfo = false;
           for(unsigned int iw=0; iw<cachedgenWIdxVec.size(); iw++) //1st cachedgenWIdxVec
	     {
	       int idxW = cachedgenWIdxVec[iw];
	       int momIdx_odx = genDecayMomIdxVec[idxW];
	       int W_odx = genDecayIdxVec[idxW];
	       bool isFromTop = false;
	       std::vector<int> perbIdxVec, perWIdxVec;
	       std::vector<TLorentzVector> perbLVec, perWLVec;
	       TLorentzVector pergenTopLVec;
	       intw momIdx = -1;
	       int decayType = 0;
	       for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++) // 2nd genDecayPdgIdVec
		 {
		   const int pdgId = genDecayPdgIdVec.at(ig);
		   if( genDecayIdxVec[ig] == momIdx_odx && std::abs(pdgId) == 6 ){ isFromTop = true; momIdx = ig; }
		   if( isFromTop )
		     {
		       if( genDecayMomIdxVec[ig] == momIdx_odx && std::abs(pdgId) == 5 )
			 {
			   perbIdxVec.push_back(ig); 
			   perbLVec.push_back(genDecayLVec.at(ig));
			 }
		       if( genDecayMomIdxVec[ig] == W_odx )
			 {
			   perWIdxVec.push_back(ig);
			   perWLVec.push_back(genDecayLVec.at(ig));
			   if( std::abs(pdgId) >= 11 && std::abs(pdgId) <=16 ) decayType =1;
			 }
		     } 
		 }// 2nd genDecayPdgIdVec ends
	       if( perbIdxVec.size() != 1 || perWIdxVec.size() != 2 || momIdx == -1 )
		 {
		   badGenInfo = true;
		 }
              std::vector<int> perIdxVec = perbIdxVec;
              std::vector<TLorentzVector> perLVec = perbLVec;
              for(unsigned int ip=0; ip<perWIdxVec.size(); ip++){
                 perIdxVec.push_back(perWIdxVec.at(ip));
                 perLVec.push_back(perWLVec.at(ip));
              }
              topDausLVec.push_back(perLVec);
              topDausIdxVec.push_back(perIdxVec);
              topDecayTypeVec.push_back(decayType);
              if( momIdx != -1 ) pergenTopLVec = genDecayLVec.at(momIdx); 
              genTopLVec.push_back(pergenTopLVec);
	     } //1st cachedgenWIdxVec ends
	   /************************************************************************************/
           /************************************************************************************/

           int cntHadDecay = 0, cntLeptDecay = 0;
           for(unsigned int id=0; id<topDecayTypeVec.size(); id++){
              if( topDecayTypeVec[id] == 0 ) cntHadDecay ++;
              else if( topDecayTypeVec[id] == 1 ) cntLeptDecay ++;
           }


	   /************************************************************************************/
           /************************************************************************************/
	   // Reco kinematics
           std::map<double, TLorentzVector> top_brJetsLVecMap;
           for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
              TLorentzVector bJet = remainbJetsLVec.at(ib);
              std::vector<TLorentzVector> rJets12Vec;
              for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                 TLorentzVector rJet1 = remainJetsLVec.at(ir);
                 TLorentzVector rbJet = bJet + rJet1;
                 if( bJet.M() > AnaConsts::lowWjetMass_ && bJet.M() < AnaConsts::highWjetMass_ ){
                    double deltaM = std::abs(rbJet.M() - mTop_);
                    if( rbJet.M() > AnaConsts::lowTopCut_ && rbJet.M() < AnaConsts::highTopCut_ ) 
                    top_brJetsLVecMap[deltaM] = rbJet;
                 }else if( rJet1.M() > AnaConsts::lowWjetMass_ && rJet1.M() < AnaConsts::highWCut_ ){
                    double deltaM = std::abs(rbJet.M() - mTop_); 
                    if( rbJet.M() > AnaConsts::lowTopCut_ && rbJet.M() < AnaConsts::highTopCut_ ) 
                    top_brJetsLVecMap[deltaM] = rbJet;
                 }
                 for(unsigned int jr=ir+1; jr<remainJetsLVec.size(); jr++){
                    TLorentzVector rJet2 = remainJetsLVec.at(jr);
                    TLorentzVector rJets12 = rJet1 + rJet2;
                    if( rJets12.M() > AnaConsts::lowWCut_ && rJets12.M() < AnaConsts::highWCut_ ){
                       rJets12Vec.push_back(rJets12);
                    }
                 }
              }
              if( bJet.M() > AnaConsts::lowTopjetMass_ && bJet.M() < AnaConsts::highTopjetMass_ ){
                 double deltaM = std::abs(bJet.M() - mTop_);
                 if( bJet.M() > AnaConsts::lowTopCut_ && bJet.M() < AnaConsts::highTopCut_ )
                 top_brJetsLVecMap[deltaM] = bJet;
              }
              for(unsigned int ic=0; ic<rJets12Vec.size(); ic++){
                 TLorentzVector rtopJet = bJet + rJets12Vec.at(ic);
                 double deltaM = std::abs(rtopJet.M() - mTop_);
                 if( rtopJet.M() > AnaConsts::lowTopCut_ && rtopJet.M() < AnaConsts::highTopCut_ ) 
                 top_brJetsLVecMap[deltaM] = rtopJet;
              }
	   }//remainbJetsLVec ends
	   /************************************************************************************/
           /************************************************************************************/

           h1_nComb_top_brJetsVec.back()->Fill(top_brJetsLVecMap.size(), evtWeight * scaleMC);

           std::vector<TLorentzVector> remainDiJetsLVec;
           for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
              for(unsigned int jr=ir+1; jr<remainJetsLVec.size(); jr++){
                 TLorentzVector perDiJetLVec = remainJetsLVec.at(ir) + remainJetsLVec.at(jr);
                 if( perDiJetLVec.M() > AnaConsts::lowWCut_ && perDiJetLVec.M() < AnaConsts::highWCut_ ){
                    remainDiJetsLVec.push_back(perDiJetLVec);
                 } 
              }
           }
           for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
              TLorentzVector perrJet = remainJetsLVec.at(ir);
              if( perrJet.M() > AnaConsts::lowWjetMass_ && perrJet.M() < AnaConsts::highWjetMass_ ) remainDiJetsLVec.push_back(perrJet);
           }
           double best_rTop_mass = -1; TLorentzVector best_rTopLVec;
           for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
              TLorentzVector bJet = remainbJetsLVec.at(ib);
              for(unsigned int id=0; id<remainDiJetsLVec.size(); id++){
                 TLorentzVector rTop = bJet + remainDiJetsLVec.at(id);
                 if( best_rTop_mass == -1 || std::abs(best_rTop_mass - mTop_) > std::abs(rTop.M() - mTop_ ) ){
                    best_rTop_mass = rTop.M(); best_rTopLVec = rTop;
                 }
              }
           }
           for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
              TLorentzVector perrJet = remainJetsLVec.at(ir);
              if( perrJet.M() > AnaConsts::lowTopjetMass_ && perrJet.M() < AnaConsts::highTopjetMass_ ){
                 if( best_rTop_mass == -1 || std::abs(best_rTop_mass - mTop_) > std::abs(perrJet.M() - mTop_ ) ){
                    best_rTop_mass = perrJet.M(); best_rTopLVec = perrJet;
                 }
              }
           }

//           if( !(best_rTop_mass > AnaConsts::lowTopCut_ && best_rTop_mass < AnaConsts::highTopCut_) ) best_rTop_mass = -1;
           if( best_rTop_mass != -1 ) h1_mass_rTopVec.back()->Fill(best_rTop_mass, evtWeight * scaleMC);
           h1_cutFlowVec.back()->Fill("found_rTop", (best_rTop_mass != -1) * evtWeight * scaleMC);
           h1_cutFlowVec.back()->Fill("found_rTop_better", (best_rTop_mass != -1 && best_rTop_mass > AnaConsts::lowTopCut_ && best_rTop_mass < AnaConsts::highTopCut_) * evtWeight * scaleMC);
	   /************************************************************************************/
           /************************************************************************************/


           std::map<double, std::vector<int> > lept_brJetIdxMap;
           std::map<double, std::vector<int> > had_brJetIdxMap;
           std::map<double, TLorentzVector> lept_brJetLVecMap, had_brJetLVecMap;
           std::map<double, TLorentzVector> ori_lept_brJetLVecMap, ori_had_brJetLVecMap;
           int cntbJet = 0;
           for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
              cntbJet++;
              TLorentzVector bJet = remainbJetsLVec.at(ib);
              for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                 TLorentzVector rJet = remainJetsLVec.at(ir);
                 TLorentzVector brJet = bJet+rJet;
                 std::vector<int> brIdxVec; brIdxVec.push_back(ib); brIdxVec.push_back(ir);
                 double dPhi_rJet_MET = rJet.DeltaPhi(metLVec);
                 double dPhi_rJet_bJet = rJet.DeltaPhi(bJet);
                 double dR_brJet = bJet.DeltaR(rJet);
                 double perMT = type3Ptr->calcMT(brJet, metLVec);
                 double perMT2 = type3Ptr->calcMT2(topLVec, brJet, metLVec);
                 double perMTr = type3Ptr->calcMT(rJet, metLVec);
                 double perMTb = type3Ptr->calcMT(bJet, metLVec);
                 if( brJet.M() < mTop_ ){
                    ori_lept_brJetLVecMap[dR_brJet] = brJet;
                    ori_had_brJetLVecMap[dR_brJet] = brJet;
                 }
//                 if( std::abs(dPhi_rJet_MET) < 1.0 && std::abs(dPhi_rJet_bJet) < 1.5 ){
                 if( std::abs(dPhi_rJet_MET) < 1.0 ){
                    if( brJet.M() < mTop_ ){
//                       lept_brJetIdxMap[std::abs(dPhi_rJet_MET)] = brIdxVec;
                       lept_brJetIdxMap[perMT] = brIdxVec;
                       lept_brJetLVecMap[perMT] = brJet;
                    }
//                 }else if( std::abs(dPhi_rJet_MET) < 1.0 && std::abs(dPhi_rJet_bJet) >= 1.5 ){
//                    lept_brJetIdxMap[perMTb] = brIdxVec;
//                    lept_brJetLVecMap[perMTb] = bJet;
                 }//else if( std::abs(dPhi_rJet_MET) >= 1.0 ){
                    if( brJet.M() > 50 && brJet.M() < mTop_ ){
//                       had_brJetIdxMap[dR_brJet] = brIdxVec;
//                       had_brJetIdxMap[perMT2] = brIdxVec;
//                       had_brJetLVecMap[perMT2] = brJet;

//                       had_brJetIdxMap[perMT] = brIdxVec;
//                       had_brJetLVecMap[perMT] = brJet;

                       had_brJetIdxMap[dR_brJet] = brIdxVec;
                       had_brJetLVecMap[dR_brJet] = brJet;
                    }
//                 }
              }
           }
/*
           if( lept_brJetLVecMap.empty() ){
              if( !had_brJetLVecMap.empty() ){
                 lept_brJetLVecMap = had_brJetLVecMap;
              }else{
                 double min_dPhi_b_MET = 999.0;
                 TLorentzVector pickedbLVec;
                 for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                    double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                    if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                       min_dPhi_b_MET = std::abs(dPhi_b_MET);
                       pickedbLVec = remainbJetsLVec.at(ib);
                    }
                 }
                 lept_brJetLVecMap[0] = pickedbLVec;
              }
           }
*/
           if( lept_brJetLVecMap.empty() ){
              double min_dPhi_b_MET = 999.0;
              TLorentzVector pickedbLVec;
              for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                 double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                 if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                    min_dPhi_b_MET = std::abs(dPhi_b_MET);
                    pickedbLVec = remainbJetsLVec.at(ib);
                 }
              }
              lept_brJetLVecMap[0] = pickedbLVec;
           }

           if( ori_lept_brJetLVecMap.empty() ){
              double min_dPhi_b_MET = 999.0;
              TLorentzVector pickedbLVec;
              for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                 double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                 if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                    min_dPhi_b_MET = std::abs(dPhi_b_MET);
                    pickedbLVec = remainbJetsLVec.at(ib);
                 }
              }
              ori_lept_brJetLVecMap[0] = pickedbLVec;
           }
           if( had_brJetLVecMap.empty() ){
              double min_dPhi_b_MET = 999.0;
              TLorentzVector pickedbLVec;
              for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                 double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                 if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                    min_dPhi_b_MET = std::abs(dPhi_b_MET);
                    pickedbLVec = remainbJetsLVec.at(ib);
                 }
              }
              had_brJetLVecMap[0] = pickedbLVec;
           }
           if( ori_had_brJetLVecMap.empty() ){
              double min_dPhi_b_MET = 999.0;
              TLorentzVector pickedbLVec;
              for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                 double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                 if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                    min_dPhi_b_MET = std::abs(dPhi_b_MET);
                    pickedbLVec = remainbJetsLVec.at(ib);
                 }
              }
              ori_had_brJetLVecMap[0] = pickedbLVec;
           }

           double top_MT = type3Ptr->calcMT(topLVec, metLVec);

           if( best_bJet_MT = 999. ){
              double min_dPhi_b_MET = 999.0;
              TLorentzVector pickedbLVec;
              for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                 double dPhi_b_MET = remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                 if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                    min_dPhi_b_MET = std::abs(dPhi_b_MET);
                    pickedbLVec = remainbJetsLVec.at(ib);
                 }
              }
              best_bJet_MT = type3Ptr->calcMT(pickedbLVec, metLVec);
           }

           if( !lept_brJetLVecMap.empty() && !isFaked_b ){
              TLorentzVector best_lept_brJet = lept_brJetLVecMap.begin()->second; 
              best_lept_brJet_MT = type3Ptr->calcMT(best_lept_brJet, metLVec);
              best_lept_brJet_MT2 = type3Ptr->calcMT2(topLVec, best_lept_brJet, metLVec);
//              if( !top_brJetsLVecMap.empty() ) best_lept_brJet_MT2 = type3Ptr->calcMT2(topLVec, top_brJetsLVecMap.begin()->second, metLVec);
//              if( best_rTop_mass != -1 ) best_lept_brJet_MT2 = type3Ptr->calcMT2(topLVec, best_rTopLVec, metLVec);
              best_lept_brJet_mTcomb = best_lept_brJet_MT + 0.6 * top_MT;
//              best_lept_brJet_mTcomb = best_lept_brJet_MT + top_MT;
           }
           if( !had_brJetLVecMap.empty() && !isFaked_b ){
              TLorentzVector best_had_brJet = had_brJetLVecMap.begin()->second;
              best_had_brJet_MT = type3Ptr->calcMT(best_had_brJet, metLVec);
              best_had_brJet_MT2 = type3Ptr->calcMT2(topLVec, best_had_brJet, metLVec);
//              if( !top_brJetsLVecMap.empty() ) best_had_brJet_MT2 = type3Ptr->calcMT2(topLVec, top_brJetsLVecMap.begin()->second, metLVec);
//              if( best_rTop_mass != -1 ){ best_had_brJet_MT2 = type3Ptr->calcMT2(topLVec, best_rTopLVec, metLVec); best_had_brJet_MT = type3Ptr->calcMT(best_rTopLVec, metLVec); }
//              best_had_brJet_mTcomb = 0.5*best_had_brJet_MT + 0.35 * top_MT;
              best_had_brJet_mTcomb = 0.5*best_had_brJet_MT + 0.5 * top_MT;
//              best_had_brJet_mTcomb = best_had_brJet_MT + top_MT;

              double data_for_pca[2];
              data_for_pca[0] = top_MT; data_for_pca[1] = best_had_brJet_MT;
              pcaVec.back()->AddRow(data_for_pca);
              double data_aft_pca[2];
              X2P(data_for_pca, data_aft_pca);
              best_aft_PCA_top_MT = data_aft_pca[0]; best_aft_PCA_had_brJet_MT = -1*data_aft_pca[1];
           }

// start adj
           std::vector<TLorentzVector> adj_jetsLVec_forTagger = jetsLVec_forTagger;
           std::vector<double> adj_recoJetsBtag_forTagger = recoJetsBtag_forTagger;

           int pickedfakebJet_idx = -1; double maxCSVS_LT_Mpt = -1;
           if( nbJets ==1 ){
              for(unsigned int ij=0; ij<adj_recoJetsBtag_forTagger.size(); ij++){
                 if( adj_jetsLVec_forTagger.at(ij).Pt() > AnaConsts::bTagArr[2] && std::abs(adj_jetsLVec_forTagger.at(ij).Eta()) < AnaConsts::bTagArr[1] ){
                    if( adj_recoJetsBtag_forTagger.at(ij) < AnaConsts::cutCSVS ){
                       if( maxCSVS_LT_Mpt ==-1 || maxCSVS_LT_Mpt < adj_recoJetsBtag_forTagger.at(ij) ){
                          pickedfakebJet_idx = ij; maxCSVS_LT_Mpt = adj_recoJetsBtag_forTagger.at(ij);
                       }
                    }
                 }
              }
//              if( pickedfakebJet_idx != -1 ){
//                 adj_recoJetsBtag_forTagger.at(pickedfakebJet_idx) = AnaConsts::cutCSVS + 1e-10;
//              }
           }

           type3Ptr2->prepareFindingBestTopCandidate(adj_jetsLVec_forTagger, adj_recoJetsBtag_forTagger);
           double bestTopMass2 = -1; int bestTopIdx2 = -1;
           for(unsigned int ic=0; ic<type3Ptr2->combSize; ic++){
              double fatJetm123 = type3Ptr2->fatJetMassVec[ic];
              if( fatJetm123 < AnaConsts::lowTopCut_ || fatJetm123 > AnaConsts::highTopCut_ ) continue;
         // Find a top fat jet passing at least one of the three criteria
              std::vector<int> fatJetPassStatusVec;
              int isTopJet = type3Ptr2->checkTopCriteria(type3Ptr2->finalCombfatJets[ic], adj_jetsLVec_forTagger, adj_recoJetsBtag_forTagger, type3Ptr2->fatJetSubMassVec[ic], fatJetm123, fatJetPassStatusVec);
              if( isTopJet ){
                 if( bestTopMass2 == -1 || std::abs(fatJetm123 - mTop_) < std::abs(bestTopMass2 - mTop_) ){
                    bestTopMass2 = fatJetm123; bestTopIdx2 = ic;
                 }
              }
           }
//           h1_cutFlowVec.back()->Fill("adj", evtWeight * scaleMC * (bestTopMass2 != -1 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );

           double adj_best_lept_brJet_MT = -1, adj_best_lept_brJet_MT2 = -1, adj_best_lept_brJet_mTcomb = -1;
           double adj_best_had_brJet_MT = -1, adj_best_had_brJet_MT2 = -1, adj_best_had_brJet_mTcomb = -1;

           if( bestTopMass2 != -1 ){
              const std::vector<int> & adj_topCombIdxVec = type3Ptr2->finalCombfatJets[bestTopIdx2];
              const TLorentzVector adj_topLVec = type3Ptr2->buildLVec(adj_jetsLVec_forTagger, adj_topCombIdxVec);

              int pickedfakebJet_idx2 = -1; double maxCSVS_LT_Mpt2 = -1;
              int adj_cntb = 0;

              for(unsigned int ij=0; ij<adj_jetsLVec_forTagger.size(); ij++){
                 if( std::find(adj_topCombIdxVec.begin(), adj_topCombIdxVec.end(), ij) != adj_topCombIdxVec.end() ) continue;

                 std::vector<TLorentzVector> dummyLVec; dummyLVec.push_back(adj_jetsLVec_forTagger.at(ij));
                 std::vector<double> dummyCSVS; dummyCSVS.push_back(adj_recoJetsBtag_forTagger.at(ij));
                 if( AnaFunctions::countCSVS(dummyLVec, dummyCSVS, AnaConsts::cutCSVS, AnaConsts::bTagArr) ) adj_cntb ++;

                 if( adj_jetsLVec_forTagger.at(ij).Pt() > AnaConsts::bTagArr[2] && std::abs(adj_jetsLVec_forTagger.at(ij).Eta()) < AnaConsts::bTagArr[1] ){
                    if( adj_recoJetsBtag_forTagger.at(ij) < AnaConsts::cutCSVS ){
                       if( maxCSVS_LT_Mpt2 ==-1 || maxCSVS_LT_Mpt2 < adj_recoJetsBtag_forTagger.at(ij) ){
                          pickedfakebJet_idx2 = ij; maxCSVS_LT_Mpt2 = adj_recoJetsBtag_forTagger.at(ij);
                       }
                    }
                 }
              }
              if( adj_cntb == 0 && pickedfakebJet_idx2 != -1 ){
                 adj_recoJetsBtag_forTagger.at(pickedfakebJet_idx2) = AnaConsts::cutCSVS + 1e-10;
              }

              std::vector<TLorentzVector> adj_remainJetsLVec; std::vector<double> adj_remainJetsCSVS;
              std::vector<TLorentzVector> adj_remainbJetsLVec; std::vector<double> adj_remainbJetsCSVS;
              for(unsigned int ij=0; ij<adj_jetsLVec_forTagger.size(); ij++){
                 if( std::find(adj_topCombIdxVec.begin(), adj_topCombIdxVec.end(), ij) != adj_topCombIdxVec.end() ) continue;

                 std::vector<TLorentzVector> dummyLVec; dummyLVec.push_back(adj_jetsLVec_forTagger.at(ij));
                 std::vector<double> dummyCSVS; dummyCSVS.push_back(adj_recoJetsBtag_forTagger.at(ij));
                 if( AnaFunctions::countCSVS(dummyLVec, dummyCSVS, AnaConsts::cutCSVS, AnaConsts::bTagArr) ){
                    adj_remainbJetsLVec.push_back(adj_jetsLVec_forTagger.at(ij));
                    adj_remainbJetsCSVS.push_back(adj_recoJetsBtag_forTagger.at(ij));
                 }else{
                    adj_remainJetsLVec.push_back(adj_jetsLVec_forTagger.at(ij));
                    adj_remainJetsCSVS.push_back(adj_recoJetsBtag_forTagger.at(ij));
                 }
              }
              std::map<double, std::vector<int> > adj_lept_brJetIdxMap;
              std::map<double, std::vector<int> > adj_had_brJetIdxMap;
              std::map<double, TLorentzVector> adj_lept_brJetLVecMap, adj_had_brJetLVecMap;
              for(unsigned int ib=0; ib<adj_remainbJetsLVec.size(); ib++){
                 TLorentzVector bJet = adj_remainbJetsLVec.at(ib);
                 for(unsigned int ir=0; ir<adj_remainJetsLVec.size(); ir++){
                    TLorentzVector rJet = adj_remainJetsLVec.at(ir);
                    TLorentzVector brJet = bJet+rJet;
                    std::vector<int> brIdxVec; brIdxVec.push_back(ib); brIdxVec.push_back(ir);
                    double dPhi_rJet_MET = rJet.DeltaPhi(metLVec);
                    double dPhi_rJet_bJet = rJet.DeltaPhi(bJet);
                    double dR_brJet = bJet.DeltaR(rJet);
                    double perMT = type3Ptr2->calcMT(brJet, metLVec);
                    double perMT2 = type3Ptr2->calcMT2(adj_topLVec, brJet, metLVec);
                    double perMTr = type3Ptr2->calcMT(rJet, metLVec);
                    double perMTb = type3Ptr2->calcMT(bJet, metLVec);
                    if( std::abs(dPhi_rJet_MET) < 1.0 ){
                       if( brJet.M() < mTop_ ){
                          adj_lept_brJetIdxMap[perMT] = brIdxVec;
                          adj_lept_brJetLVecMap[perMT] = brJet;
                       }
                    }
                    if( brJet.M() > 50 && brJet.M() < mTop_ ){
                       adj_had_brJetIdxMap[dR_brJet] = brIdxVec;
                       adj_had_brJetLVecMap[dR_brJet] = brJet;
                    }
                 }
              }

              if( adj_lept_brJetLVecMap.empty() ){
                 double min_dPhi_b_MET = 999.0;
                 TLorentzVector pickedbLVec;
                 for(unsigned int ib=0; ib<adj_remainbJetsLVec.size(); ib++){
                    double dPhi_b_MET = adj_remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                    if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                       min_dPhi_b_MET = std::abs(dPhi_b_MET);
                       pickedbLVec = adj_remainbJetsLVec.at(ib);
                    }
                 }
                 adj_lept_brJetLVecMap[0] = pickedbLVec;
              }

              if( adj_had_brJetLVecMap.empty() ){
                 double min_dPhi_b_MET = 999.0;
                 TLorentzVector pickedbLVec;
                 for(unsigned int ib=0; ib<adj_remainbJetsLVec.size(); ib++){
                    double dPhi_b_MET = adj_remainbJetsLVec.at(ib).DeltaPhi(metLVec);
                    if( min_dPhi_b_MET > std::abs(dPhi_b_MET) ){
                       min_dPhi_b_MET = std::abs(dPhi_b_MET);
                       pickedbLVec = adj_remainbJetsLVec.at(ib);
                    }
                 }
                 adj_had_brJetLVecMap[0] = pickedbLVec;
              }

              double adj_top_MT = type3Ptr2->calcMT(adj_topLVec, metLVec);

              TLorentzVector adj_best_lept_brJet = adj_lept_brJetLVecMap.begin()->second;
              adj_best_lept_brJet_MT = type3Ptr2->calcMT(adj_best_lept_brJet, metLVec);
              adj_best_lept_brJet_MT2 = type3Ptr2->calcMT2(adj_topLVec, adj_best_lept_brJet, metLVec);
              adj_best_lept_brJet_mTcomb = adj_best_lept_brJet_MT + 0.6 * adj_top_MT;

              TLorentzVector adj_best_had_brJet = adj_had_brJetLVecMap.begin()->second;
              adj_best_had_brJet_MT = type3Ptr2->calcMT(adj_best_had_brJet, metLVec);
              adj_best_had_brJet_MT2 = type3Ptr2->calcMT2(adj_topLVec, adj_best_had_brJet, metLVec);
              adj_best_had_brJet_mTcomb = 0.5*adj_best_had_brJet_MT + 0.5 * adj_top_MT;

              if( best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){
                 h1_adj_MT_lept_brJetVec.back()->Fill(adj_best_lept_brJet_MT, evtWeight * scaleMC);
                 h1_adj_MT2_lept_brJetVec.back()->Fill(adj_best_lept_brJet_MT2, evtWeight * scaleMC);
                 h2_adj_MT_vs_MT2_lept_brJetVec.back()->Fill(adj_best_lept_brJet_MT2, adj_best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_adj_MT_lept_vs_top_brJetVec.back()->Fill(adj_top_MT, adj_best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_adj_mTcomb_vs_MT2_lept_brJetVec.back()->Fill(adj_best_lept_brJet_MT2, adj_best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 h2_adj_mTcomb_vs_MT_lept_brJetVec.back()->Fill(adj_best_lept_brJet_MT, adj_best_lept_brJet_mTcomb, evtWeight * scaleMC);

                 h1_adj_MT_had_brJetVec.back()->Fill(adj_best_had_brJet_MT, evtWeight * scaleMC);
                 h1_adj_MT2_had_brJetVec.back()->Fill(adj_best_had_brJet_MT2, evtWeight * scaleMC);
                 h2_adj_MT_vs_MT2_had_brJetVec.back()->Fill(adj_best_had_brJet_MT2, adj_best_had_brJet_MT, evtWeight * scaleMC);
                 h2_adj_MT_had_vs_top_brJetVec.back()->Fill(adj_top_MT, adj_best_had_brJet_MT, evtWeight * scaleMC);
                 h2_adj_mTcomb_vs_MT2_had_brJetVec.back()->Fill(adj_best_had_brJet_MT2, adj_best_had_brJet_mTcomb, evtWeight * scaleMC);
                 h2_adj_mTcomb_vs_MT_had_brJetVec.back()->Fill(adj_best_had_brJet_MT, adj_best_had_brJet_mTcomb, evtWeight * scaleMC);

                 if( adj_best_had_brJet_MT2 != adj_best_lept_brJet_MT2 ) h2_adj_MT2_lept_vs_had_brJetVec.back()->Fill(adj_best_had_brJet_MT2, adj_best_lept_brJet_MT2, evtWeight * scaleMC);
                 else h1_adj_MT2_same_lept_had_brJetVec.back()->Fill(adj_best_lept_brJet_MT2, evtWeight * scaleMC);
                 if( adj_best_had_brJet_MT != adj_best_lept_brJet_MT ){
                    h2_adj_MT_lept_vs_had_brJetVec.back()->Fill(adj_best_had_brJet_MT, adj_best_lept_brJet_MT, evtWeight * scaleMC);
                    h2_adj_mTcomb_lept_vs_had_brJetVec.back()->Fill(adj_best_had_brJet_mTcomb, adj_best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 } else{
                    h1_adj_MT_same_lept_had_brJetVec.back()->Fill(adj_best_lept_brJet_MT, evtWeight * scaleMC);
                    h1_adj_mTcomb_same_lept_had_brJetVec.back()->Fill(adj_best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 }
                 h2_adj_MT2_lept_vs_MT_had_brJetVec.back()->Fill(adj_best_had_brJet_MT, adj_best_lept_brJet_MT2, evtWeight * scaleMC);
                 h2_adj_MT_lept_vs_MT2_had_brJetVec.back()->Fill(adj_best_had_brJet_MT2, adj_best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_adj_MT_lept_vs_mTcomb_had_brJetVec.back()->Fill(adj_best_had_brJet_mTcomb, adj_best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_adj_mTcomb_lept_vs_MT_had_brJetVec.back()->Fill(adj_best_had_brJet_MT, adj_best_lept_brJet_mTcomb, evtWeight * scaleMC);
              }
           }
           if( best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){
//              if( bestTopMass2 != -1 )
//                 h1_cutFlowVec.back()->Fill("adj_cuts", evtWeight * scaleMC * ( adj_best_had_brJet_mTcomb >= 450 && (adj_best_lept_brJet_MT >= 200 || adj_best_had_brJet_mTcomb >=600) ) );
//              else h1_cutFlowVec.back()->Fill("adj_cuts", evtWeight * scaleMC );
           }
// end adj

           std::map<double, TLorentzVector> tJet_for_mtLVecMap;
           std::map<double, TLorentzVector> tDiJets_for_mtLVecMap;
           for(unsigned int it=0; it<topCombIdxVec.size(); it++){
              int idx = topCombIdxVec.at(it);
              double dPhi = metLVec.DeltaPhi(jetsLVec_forTagger.at(idx));
              if( std::abs(dPhi) < 1.5 ){
                 double tJet_MT = type3Ptr->calcMT(jetsLVec_forTagger.at(idx), metLVec);
                 tJet_for_mtLVecMap[tJet_MT] = jetsLVec_forTagger.at(idx);
              }
              for(unsigned int jt=it+1; jt<topCombIdxVec.size(); jt++){
                 int jdx = topCombIdxVec.at(jt);
                 TLorentzVector diJet = jetsLVec_forTagger.at(idx) + jetsLVec_forTagger.at(jdx);
                 double tDiJets_MT = type3Ptr->calcMT(diJet, metLVec);
                 tDiJets_for_mtLVecMap[tDiJets_MT] = diJet;
              }
           }

           if( !badGenInfo ){
//              h1_cutFlowVec.back()->Fill("1b", evtWeight * scaleMC * (cntHadDecay ==1 && cntLeptDecay == 1 & nbJets ==1 && (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) ) );
//              h1_cutFlowVec.back()->Fill("1b_faked", evtWeight * scaleMC * (pickedfakebJet_idx != -1 && cntHadDecay ==1 && cntLeptDecay == 1 & nbJets ==1 && (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) ) );

              if( pickedfakebJet_idx != -1 && cntHadDecay ==1 && cntLeptDecay == 1 && nbJets ==1 && (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) ){
                 double fakeb_csvs = recoJetsBtag_forTagger.at(pickedfakebJet_idx);
                 h1_csvs_fakebVec.back()->Fill(fakeb_csvs, evtWeight * scaleMC);
              }
              for(unsigned int id=0; id<topDecayTypeVec.size(); id++){
                 if( pickedfakebJet_idx != -1 && cntHadDecay ==1 && cntLeptDecay == 1 && nbJets ==1 && (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) ){
                    TLorentzVector fakebLVec = jetsLVec_forTagger.at(pickedfakebJet_idx);
                    TLorentzVector bLVec = topDausLVec[id][0];
                    double dR_fb = fakebLVec.DeltaR(bLVec);
//                    h1_cutFlowVec.back()->Fill("mh_ftb", evtWeight * scaleMC * (topDecayTypeVec[id] == 0 && dR_fb < 0.4 ) );
//                    h1_cutFlowVec.back()->Fill("mh_frb", evtWeight * scaleMC * (topDecayTypeVec[id] == 1 && dR_fb < 0.4 ) );
                    int cntmatch =0;
                    if( topDecayTypeVec[id] == 0 && dR_fb < 0.4 ) cntmatch++;
                    if( topDecayTypeVec[id] == 1 && dR_fb < 0.4 ) cntmatch++;
//                    h1_cutFlowVec.back()->Fill("mh_0fxb", evtWeight * scaleMC * (cntmatch == 0 ));
                 }
                 if( topDecayTypeVec[id] == 0 ){
                    TLorentzVector leadLVec, secLVec, thrLVec;
                    TLorentzVector bLVec = topDausLVec[id][0];
                    TLorentzVector Wdau1LVec = topDausLVec[id][1];
                    TLorentzVector Wdau2LVec = topDausLVec[id][2];
                    if( Wdau1LVec.Pt() < Wdau2LVec.Pt() ){ TLorentzVector tmpLVec = Wdau1LVec; Wdau1LVec = Wdau2LVec; Wdau2LVec = tmpLVec; }
                    if( bLVec.Pt() > Wdau1LVec.Pt() ){ leadLVec = bLVec; secLVec = Wdau1LVec; thrLVec = Wdau2LVec; }
                    else if( bLVec.Pt() > Wdau2LVec.Pt() ){ leadLVec = Wdau1LVec; secLVec = bLVec; thrLVec = Wdau2LVec; }
                    else{ leadLVec = Wdau1LVec; secLVec = Wdau2LVec; thrLVec = bLVec; }
                    TLorentzVector b1LVec = bLVec + Wdau1LVec;
                    TLorentzVector b2LVec = bLVec + Wdau2LVec;
                    double mb1_sq = b1LVec.M() * b1LVec.M(); double mb2_sq = b2LVec.M() * b2LVec.M();
                    double mb1 = b1LVec.M(); double mb2 = b2LVec.M();
                    h2_dalitzVec.back()->Fill(mb2, mb1, evtWeight * scaleMC);
                    
                    double m12 = (leadLVec + secLVec).M(), m13 = (leadLVec + thrLVec).M(), m23 = (secLVec + thrLVec).M();
                    h2_gen_m23overm123vsarctanm13overm12Vec.back()->Fill( TMath::ATan(m13/m12), m23/mTop_, evtWeight * scaleMC);
//                    if( cntHadDecay ==1 && cntLeptDecay ==1 ){
                    if( cntHadDecay ==1 && cntLeptDecay ==1 && topCombIdxVec.size() == 3 ){
                       int matchlead =0, matchsec =0, matchthr = 0;
                       int permatchb =0, permatchd1 = 0, permatchd2 = 0;
                       int permatchj1 = 0, permatchj2 = 0, permatchj3 = 0;
                       int permatchd1_aJet = 0, permatchd2_aJet = 0;
                       for(unsigned int ic=0; ic<topCombIdxVec.size(); ic++){
                          const int idx = topCombIdxVec.at(ic);
                          TLorentzVector perJet = jetsLVec_forTagger.at(idx);
                          double dR_jl = perJet.DeltaR(leadLVec);
                          double dR_js = perJet.DeltaR(secLVec);
                          double dR_jt = perJet.DeltaR(thrLVec);
                          if( dR_jl < 0.4 ) matchlead = 1;
                          if( dR_js < 0.4 ) matchsec = 1;
                          if( dR_jt < 0.4 ) matchthr = 1;
                          double dR_jb = perJet.DeltaR(bLVec);
                          double dR_jd1 = perJet.DeltaR(Wdau1LVec);
                          double dR_jd2 = perJet.DeltaR(Wdau2LVec);
                          if( dR_jb < 0.4 ) permatchb = 1;
                          if( dR_jd1 < 0.4 ) permatchd1 = 1;
                          if( dR_jd2 < 0.4 ) permatchd2 = 1;

                          if( ic == 0 && (dR_jb < 0.4 || dR_jd1 < 0.4 || dR_jd2 < 0.4) ) permatchj1 = 1;
                          if( ic == 1 && (dR_jb < 0.4 || dR_jd1 < 0.4 || dR_jd2 < 0.4) ) permatchj2 = 1;
                          if( ic == 2 && (dR_jb < 0.4 || dR_jd1 < 0.4 || dR_jd2 < 0.4) ) permatchj3 = 1;
                          if( dR_jd1 < 0.4 ) permatchd1_aJet = 1;
                          if( dR_jd2 < 0.4 ) permatchd2_aJet = 1;
                       }
                       int permatchb_rJet =0, permatchb_rbJet = 0;
                       for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                          TLorentzVector &perJet = remainJetsLVec.at(ir);
                          double dR_rJet = perJet.DeltaR(bLVec);
                          if( dR_rJet < 0.4 ) permatchb_rJet = 1;
                          double dR_jd1 = perJet.DeltaR(Wdau1LVec);
                          double dR_jd2 = perJet.DeltaR(Wdau2LVec);
                          if( dR_jd1 < 0.4 ) permatchd1_aJet = 1;
                          if( dR_jd2 < 0.4 ) permatchd2_aJet = 1;
                       }
                       for(unsigned int ir=0; ir<remainbJetsLVec.size(); ir++){
                          TLorentzVector &perJet = remainbJetsLVec.at(ir);
                          double dR_rJet = perJet.DeltaR(bLVec);
                          if( dR_rJet < 0.4 ) permatchb_rbJet = 1;
                          double dR_jd1 = perJet.DeltaR(Wdau1LVec);
                          double dR_jd2 = perJet.DeltaR(Wdau2LVec);
                          if( dR_jd1 < 0.4 ) permatchd1_aJet = 1;
                          if( dR_jd2 < 0.4 ) permatchd2_aJet = 1;
                       }
                       if( best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){
                          h1_pt_gentbVec.back()->Fill(bLVec.Pt(), evtWeight * scaleMC);
                          h1_eta_gentbVec.back()->Fill(bLVec.Eta(), evtWeight * scaleMC);
                       }
/*                       
                       int match0 = 0, matchb = 0, matchd1 = 0, matchd2 = 0, matchbd1 = 0, matchbd2 = 0, matchd12 = 0, matchbd12 = 0;
                       if( ! permatchb && ! permatchd1 && ! permatchd2 ) match0 = 1;
                       if( permatchb && ! permatchd1 && ! permatchd2 ) matchb = 1;
                       if( ! permatchb && permatchd1 && ! permatchd2 ) matchd1 = 1;
                       if( ! permatchb && ! permatchd1 && permatchd2 ) matchd2 = 1;
                       if( permatchb && permatchd1 && ! permatchd2 ) matchbd1 = 1;
                       if( permatchb && ! permatchd1 && permatchd2 ) matchbd2 = 1;
                       if( ! permatchb && permatchd1 && permatchd2 ) matchd12 = 1;
                       if( permatchb && permatchd1 && permatchd2 ) matchbd12 = 1;
                       h1_cutFlowVec.back()->Fill("1h_1l", evtWeight * scaleMC * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mh0", evtWeight * scaleMC * match0 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhb", evtWeight * scaleMC * matchb * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhd1", evtWeight * scaleMC * matchd1 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhd2", evtWeight * scaleMC * matchd2 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhbd1", evtWeight * scaleMC * matchbd1 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhbd2", evtWeight * scaleMC * matchbd2 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhd12", evtWeight * scaleMC * matchd12 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("mhbd12", evtWeight * scaleMC * matchbd12 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));

                       h1_cutFlowVec.back()->Fill("pMT_1h_1l", evtWeight * scaleMC * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mh0", evtWeight * scaleMC * match0 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhb", evtWeight * scaleMC * matchb * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhd1", evtWeight * scaleMC * matchd1 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhd2", evtWeight * scaleMC * matchd2 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhbd1", evtWeight * scaleMC * matchbd1 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhbd2", evtWeight * scaleMC * matchbd2 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhd12", evtWeight * scaleMC * matchd12 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("pMT_mhbd12", evtWeight * scaleMC * matchbd12 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));

                       h1_cutFlow_auxVec.back()->Fill("1h_1l", evtWeight * scaleMC * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_1h_1l", evtWeight * scaleMC * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
/*
                       int match0js = 0, matchj1 = 0, matchj2 = 0, matchj3 = 0, matchj12 = 0, matchj13 = 0, matchj23 = 0, matchj123 = 0;
                       if( !permatchj1 && !permatchj2 && !permatchj3 ) match0js =1;
                       if( permatchj1 && !permatchj2 && !permatchj3 ) matchj1 =1;
                       if( !permatchj1 && permatchj2 && !permatchj3 ) matchj2 =1;
                       if( !permatchj1 && !permatchj2 && permatchj3 ) matchj3 =1;
                       if( permatchj1 && permatchj2 && !permatchj3 ) matchj12 =1;
                       if( permatchj1 && !permatchj2 && permatchj3 ) matchj13 =1;
                       if( !permatchj1 && permatchj2 && permatchj3 ) matchj23 =1;
                       if( permatchj1 && permatchj2 && permatchj3 ) matchj123 =1;
                       h1_cutFlow_auxVec.back()->Fill("mh0", evtWeight * scaleMC * match0js * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj1", evtWeight * scaleMC * matchj1 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj2", evtWeight * scaleMC * matchj2 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj3", evtWeight * scaleMC * matchj3 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj12", evtWeight * scaleMC * matchj12 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj13", evtWeight * scaleMC * matchj13 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj23", evtWeight * scaleMC * matchj23 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("mhj123", evtWeight * scaleMC * matchj123 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));

                       h1_cutFlow_auxVec.back()->Fill("pMT_mh0", evtWeight * scaleMC * match0js * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj1", evtWeight * scaleMC * matchj1 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj2", evtWeight * scaleMC * matchj2 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj3", evtWeight * scaleMC * matchj3 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj12", evtWeight * scaleMC * matchj12 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj13", evtWeight * scaleMC * matchj13 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj23", evtWeight * scaleMC * matchj23 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_mhj123", evtWeight * scaleMC * matchj123 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
                       int genb_matchtj=0, genb_matchrb=0, genb_matchrj=0, genb_match0 =0;
                       int match0_aJet =0, matchb_aJet =0, matchd1_aJet =0, matchd2_aJet =0, matchbd1_aJet =0, matchbd2_aJet =0, matchd12_aJet =0, matchbd12_aJet =0;
                       if( permatchb ) genb_matchtj = 1;
                       if( permatchb_rbJet ) genb_matchrb = 1;
                       if( permatchb_rJet ) genb_matchrj = 1;
                       if( !permatchb && !permatchb_rbJet && !permatchb_rJet ) genb_match0 = 1;

                       bool ismatchb_aJet = false;
                       if( permatchb || permatchb_rbJet || permatchb_rJet ) ismatchb_aJet=true;
                       if( !ismatchb_aJet && !permatchd1_aJet && !permatchd2_aJet ) match0_aJet =1;
                       if( ismatchb_aJet && !permatchd1_aJet && !permatchd2_aJet ) matchb_aJet =1;
                       if( !ismatchb_aJet && permatchd1_aJet && !permatchd2_aJet ) matchd1_aJet =1;
                       if( !ismatchb_aJet && !permatchd1_aJet && permatchd2_aJet ) matchd2_aJet =1;
                       if( ismatchb_aJet && permatchd1_aJet && !permatchd2_aJet ) matchbd1_aJet =1;
                       if( ismatchb_aJet && !permatchd1_aJet && permatchd2_aJet ) matchbd2_aJet =1;
                       if( !ismatchb_aJet && permatchd1_aJet && permatchd2_aJet ) matchd12_aJet =1;
                       if( ismatchb_aJet && permatchd1_aJet && permatchd2_aJet ) matchbd12_aJet =1;
/*                    
                       h1_cutFlow_auxVec.back()->Fill("gb_mhtj", evtWeight * scaleMC * genb_matchtj * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("gb_mhrb", evtWeight * scaleMC * genb_matchrb * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("gb_mhrj", evtWeight * scaleMC * genb_matchrj * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("gb_mh0", evtWeight * scaleMC * genb_match0 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_gb_mhtj", evtWeight * scaleMC * genb_matchtj * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_gb_mhrb", evtWeight * scaleMC * genb_matchrb * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_gb_mhrj", evtWeight * scaleMC * genb_matchrj * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_gb_mh0", evtWeight * scaleMC * genb_match0 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
/*
                       h1_cutFlow_auxVec.back()->Fill("amh0", evtWeight * scaleMC * match0_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhb", evtWeight * scaleMC * matchb_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhd1", evtWeight * scaleMC * matchd1_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhd2", evtWeight * scaleMC * matchd2_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhbd1", evtWeight * scaleMC * matchbd1_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhbd2", evtWeight * scaleMC * matchbd2_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("amhbd12", evtWeight * scaleMC * matchbd12_aJet * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
/*
                       h1_cutFlow_auxVec.back()->Fill("pMT_amh0", evtWeight * scaleMC * match0_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhb", evtWeight * scaleMC * matchb_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhd1", evtWeight * scaleMC * matchd1_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhd2", evtWeight * scaleMC * matchd2_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhbd1", evtWeight * scaleMC * matchbd1_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhbd2", evtWeight * scaleMC * matchbd2_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_amhbd12", evtWeight * scaleMC * matchbd12_aJet * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
/*
                       int match0 = 0, match1 = 0, match2 = 0, match3 = 0, match12 = 0, match13 = 0, match23 = 0, match123 = 0;
                       if( ! matchlead && ! matchsec && ! matchthr ) match0 = 1;
                       if( matchlead && ! matchsec && ! matchthr ) match1 = 1;
                       if( ! matchlead && matchsec && ! matchthr ) match2 = 1;
                       if( ! matchlead && ! matchsec && matchthr ) match3 = 1;
                       if( matchlead && matchsec && ! matchthr ) match12 = 1;
                       if( matchlead && ! matchsec && matchthr ) match13 = 1;
                       if( !matchlead && matchsec && matchthr ) match23 = 1;
                       if( matchlead && matchsec && matchthr ) match123 = 1;
                       h1_cutFlowVec.back()->Fill("match0", evtWeight * scaleMC * match0 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match1", evtWeight * scaleMC * match1 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match2", evtWeight * scaleMC * match2 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match3", evtWeight * scaleMC * match3 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match12", evtWeight * scaleMC * match12 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match13", evtWeight * scaleMC * match13 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match23", evtWeight * scaleMC * match23 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("match123", evtWeight * scaleMC * match123 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));

                       h1_cutFlowVec.back()->Fill("passMT_match0", evtWeight * scaleMC * match0 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match1", evtWeight * scaleMC * match1 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match2", evtWeight * scaleMC * match2 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match3", evtWeight * scaleMC * match3 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match12", evtWeight * scaleMC * match12 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match13", evtWeight * scaleMC * match13 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match23", evtWeight * scaleMC * match23 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlowVec.back()->Fill("passMT_match123", evtWeight * scaleMC * match123 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
                    }
                 }
                 if( topDecayTypeVec[id] == 1 ){
                    TLorentzVector bLVec = topDausLVec[id][0];
                    TLorentzVector Wdau1LVec = topDausLVec[id][1];
                    TLorentzVector Wdau2LVec = topDausLVec[id][2];
                    if( Wdau1LVec.Pt() < Wdau2LVec.Pt() ){ TLorentzVector tmpLVec = Wdau1LVec; Wdau1LVec = Wdau2LVec; Wdau2LVec = tmpLVec; }
                    if( cntHadDecay ==1 && cntLeptDecay ==1 ){
                       int permatchb_tJet = 0, permatchb_rbJet = 0, permatchb_rrJet = 0;
                       for(unsigned int ic=0; ic<topCombIdxVec.size(); ic++){
                          const int idx = topCombIdxVec.at(ic);
                          const TLorentzVector & perJet = jetsLVec_forTagger.at(idx);
                          double dR_jb = perJet.DeltaR(bLVec);
                          if( dR_jb < 0.4 ) permatchb_tJet = 1;
                       }
                       for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                          const TLorentzVector &perJet = remainJetsLVec.at(ir);
                          double dR_rJet = perJet.DeltaR(bLVec);
                          if( dR_rJet < 0.4 ) permatchb_rrJet = 1;
                       }
                       for(unsigned int ir=0; ir<remainbJetsLVec.size(); ir++){
                          const TLorentzVector &perJet = remainbJetsLVec.at(ir);
                          double dR_rJet = perJet.DeltaR(bLVec);
                          if( dR_rJet < 0.4 ) permatchb_rbJet = 1;
                       }

                       int rgb_match0 = 0, rgb_matchtj = 0, rgb_matchrb = 0, rgb_matchrj = 0;
                       if( !permatchb_tJet && !permatchb_rbJet && !permatchb_rrJet ) rgb_match0 = 1;
                       if( permatchb_tJet && !permatchb_rbJet && !permatchb_rrJet ) rgb_matchtj = 1;
                       if( !permatchb_tJet && permatchb_rbJet && !permatchb_rrJet ) rgb_matchrb = 1;
                       if( !permatchb_tJet && !permatchb_rbJet && permatchb_rrJet ) rgb_matchrj = 1;
/*
                       h1_cutFlow_auxVec.back()->Fill("rgb_mhtj", evtWeight * scaleMC * rgb_matchtj * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("rgb_mhrb", evtWeight * scaleMC * rgb_matchrb * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("rgb_mhrj", evtWeight * scaleMC * rgb_matchrj * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("rgb_mh0", evtWeight * scaleMC * rgb_match0 * !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_rgb_mhtj", evtWeight * scaleMC * rgb_matchtj * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_rgb_mhrb", evtWeight * scaleMC * rgb_matchrb * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_rgb_mhrj", evtWeight * scaleMC * rgb_matchrj * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
                       h1_cutFlow_auxVec.back()->Fill("pMT_rgb_mh0", evtWeight * scaleMC * rgb_match0 * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)));
*/
                       if( best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){
                          h1_pt_genrbVec.back()->Fill(bLVec.Pt(), evtWeight * scaleMC);
                          h1_eta_genrbVec.back()->Fill(bLVec.Eta(), evtWeight * scaleMC);
                          if( rgb_match0 ){
                             h1_pt_genrb_match0Vec.back()->Fill(bLVec.Pt(), evtWeight * scaleMC);
                             h1_eta_genrb_match0Vec.back()->Fill(bLVec.Eta(), evtWeight * scaleMC);
                          }
                       }
                    }
                 }
              }
           }

//           h1_cutFlowVec.back()->Fill("aftMTcut", evtWeight * scaleMC * (best_had_brJet_MT >= 200 && best_lept_brJet_MT >= 200) );
//           h1_cutFlowVec.back()->Fill("aftMTcombMTcut", evtWeight * scaleMC * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
//           h1_cutFlowVec.back()->Fill("rTop", evtWeight * scaleMC * (best_rTop_mass != -1 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           if( !badGenInfo && !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450)) ){
/*
              h1_cutFlowVec.back()->Fill("triplet", evtWeight * scaleMC * (passStatusVec.size() == 3) );
              h1_cutFlowVec.back()->Fill("pass1", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass1 );
              h1_cutFlowVec.back()->Fill("pass2", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass2 );
              h1_cutFlowVec.back()->Fill("pass3", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass3 );
              h1_cutFlowVec.back()->Fill("pass12", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass12 );
              h1_cutFlowVec.back()->Fill("pass13", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass13 );
              h1_cutFlowVec.back()->Fill("pass23", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass23 );
*/
           }

           if( !badGenInfo ){
//           if( !badGenInfo && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){

//              h1_cutFlowVec.back()->Fill("goodGenInfo", evtWeight * scaleMC);
/*
              h1_cutFlowVec.back()->Fill("aftMT_triplet", evtWeight * scaleMC * (passStatusVec.size() == 3) );
              h1_cutFlowVec.back()->Fill("aftMT_pass1", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass1 );
              h1_cutFlowVec.back()->Fill("aftMT_pass2", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass2 );
              h1_cutFlowVec.back()->Fill("aftMT_pass3", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass3 );
              h1_cutFlowVec.back()->Fill("aftMT_pass12", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass12 );
              h1_cutFlowVec.back()->Fill("aftMT_pass13", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass13 );
              h1_cutFlowVec.back()->Fill("aftMT_pass23", evtWeight * scaleMC * (passStatusVec.size() == 3) * pass23 );
*/
              if(genTopLVec.size()>=2 ){
                 double dPhi_genTops = genTopLVec[0].DeltaPhi(genTopLVec[1]);
                 double dPhi_genbs = topDausLVec[0].front().DeltaPhi(topDausLVec[1].front());
                 double dPhi_genb_genTop = genTopLVec[0].DeltaPhi(topDausLVec[1].front());
                 h1_dPhi_genTopsVec.back()->Fill(dPhi_genTops, evtWeight * scaleMC);
                 h1_dPhi_genbsVec.back()->Fill(dPhi_genbs, evtWeight * scaleMC);
                 h1_dPhi_genb_genTopVec.back()->Fill(dPhi_genb_genTop, evtWeight * scaleMC);
              }

              if( !tJet_for_mtLVecMap.empty() ){
//                 h1_cutFlowVec.back()->Fill("tJet_MT", evtWeight * scaleMC);
                 h1_MT_tJetVec.back()->Fill(tJet_for_mtLVecMap.begin()->first, evtWeight * scaleMC);
              }
              if( !tDiJets_for_mtLVecMap.empty() ){
//                 h1_cutFlowVec.back()->Fill("tDiJets_MT", evtWeight * scaleMC);
                 h1_MT_tDiJetsVec.back()->Fill(tDiJets_for_mtLVecMap.begin()->first, evtWeight * scaleMC);
              }

              double dR_b_top = topLVec.DeltaR(remainbJetsLVec.front());
              double dPhi_b_top = topLVec.DeltaPhi(remainbJetsLVec.front());
              double mt_b = type3Ptr->calcMT(remainbJetsLVec.front(), metLVec);
              h1_dR_b_topVec.back()->Fill(dR_b_top, evtWeight * scaleMC);
              h1_dPhi_b_topVec.back()->Fill(dPhi_b_top, evtWeight * scaleMC);
              h2_mt_b_vs_dPhi_b_topVec.back()->Fill(dPhi_b_top, mt_b, evtWeight * scaleMC);
   
              if( cntHadDecay ==1 && cntLeptDecay ==1 ){
                 h1_cutFlowVec.back()->Fill("1had_1lept", evtWeight * scaleMC);
//              if( cntHadDecay ==2 ){
//                 h1_cutFlowVec.back()->Fill("2hads", evtWeight * scaleMC);
                 int rJets_cnt = (int)remainbJetsLVec.size() + (int)remainJetsLVec.size();
//                 if(rJets_cnt == 1) h1_cutFlowVec.back()->Fill("rJets_EQ1", evtWeight * scaleMC);
//                 if(rJets_cnt == 2) h1_cutFlowVec.back()->Fill("rJets_EQ2", evtWeight * scaleMC);
//                 if(rJets_cnt >= 3) h1_cutFlowVec.back()->Fill("rJets_LE3", evtWeight * scaleMC);

                 const std::vector<TLorentzVector> & allJetsLVec = tr->getVec<TLorentzVector>("jetsLVec");

                 for(unsigned int id=0; id<topDecayTypeVec.size(); id++){
                    std::vector<int> perDauIdxVec = topDausIdxVec.at(id);
                    std::vector<TLorentzVector> perDauLVec = topDausLVec.at(id);
                    TLorentzVector pergenTopLVec = genTopLVec.at(id);
                    TLorentzVector pickedLeptLVec; int pickedLeptPdgId;
                    TLorentzVector pickedNuLVec; int pickedNuPdgId;
                    TLorentzVector pickedbLVec; int pickedbPdgId;
                    if( topDecayTypeVec[id] == 1 ){
                       for(unsigned int idau=0; idau< perDauIdxVec.size(); idau++){
                          int pdgId = genDecayPdgIdVec.at(perDauIdxVec.at(idau));
                          if( std::abs(pdgId) == 11 || std::abs(pdgId) == 13 || std::abs(pdgId) == 15){
                             pickedLeptLVec = perDauLVec.at(idau);
                             pickedLeptPdgId = pdgId;
                          }
                          if( std::abs(pdgId) == 12 || std::abs(pdgId) == 14 || std::abs(pdgId) == 16){
                             pickedNuLVec = perDauLVec.at(idau);
                             pickedNuPdgId = pdgId;
                          }
                       }
// 0: out-of-acceptance   1: not reconstructed   2: not isolated
                       int lostlept_cat = -1;
                       if( pickedLeptLVec.Pt() < 5.0 || std::abs(pickedLeptLVec.Eta()) > 2.4 ) lostlept_cat = 0;
                       else{
                          int isRecoMuon = 0, isRecoEle = 0, isRecoIsoTrk = 0;
                          for(unsigned int im=0; im<muonsLVec.size(); im++){
                             double dR = pickedLeptLVec.DeltaR(muonsLVec.at(im));
                             if( dR < 0.3 ) isRecoMuon ++;
                          }
                          for(unsigned int ie=0; ie<elesLVec.size(); ie++){
                             double dR = pickedLeptLVec.DeltaR(elesLVec.at(ie));
                             if( dR < 0.3 ) isRecoEle ++;
                          }
                          for(unsigned int is=0; is<loose_isoTrksLVec.size(); is++){
                             double dR = pickedLeptLVec.DeltaR(loose_isoTrksLVec.at(is));
                             if( dR < 0.3 ) isRecoIsoTrk ++;
                          }
                          if( !isRecoMuon && !isRecoEle && !isRecoIsoTrk ) lostlept_cat = 1;
                          else lostlept_cat = 2;
                       }

                       h1_cutFlowVec.back()->Fill("ll_acc", (lostlept_cat ==0) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("ll_rec", (lostlept_cat ==1) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("ll_iso", (lostlept_cat ==2) * evtWeight * scaleMC);

                       double dPhi_lept_met = metLVec.DeltaPhi(pickedLeptLVec);
                       double dPhi_nu_met = metLVec.DeltaPhi(pickedNuLVec);
                       double dPhi_sum_lept_nu_met = metLVec.DeltaPhi( (pickedLeptLVec+pickedNuLVec) );
                       h1_dPhi_lept_metVec.back()->Fill(dPhi_lept_met, evtWeight * scaleMC);
                       h1_dPhi_nu_metVec.back()->Fill(dPhi_nu_met, evtWeight * scaleMC);
                       h1_dPhi_sum_lept_nu_metVec.back()->Fill(dPhi_sum_lept_nu_met, evtWeight * scaleMC);

                       double minDR_rJet_lept = 999.0; TLorentzVector pickedrJetLVec;
                       for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                          double dR = pickedLeptLVec.DeltaR(remainJetsLVec.at(ir));
                          if( minDR_rJet_lept > dR ){ minDR_rJet_lept = dR; pickedrJetLVec = remainJetsLVec.at(ir); }
                       }
                       for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                          double dR = pickedLeptLVec.DeltaR(remainbJetsLVec.at(ib));
                          if( minDR_rJet_lept > dR ){ minDR_rJet_lept = dR; pickedrJetLVec = remainbJetsLVec.at(ib); }
                       }

                       h1_cutFlowVec.back()->Fill("close_rJet_lept", (minDR_rJet_lept < 0.4) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_rJet_acc", (minDR_rJet_lept < 0.4) * (lostlept_cat == 0) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_rJet_rec", (minDR_rJet_lept < 0.4) * (lostlept_cat == 1) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_rJet_iso", (minDR_rJet_lept < 0.4) * (lostlept_cat == 2) * evtWeight * scaleMC);

                       if( minDR_rJet_lept < 0.4 ){
                          double dPhi_rJet_met = pickedrJetLVec.DeltaPhi(metLVec);
                          double dPhi_rJet_b = remainbJetsLVec.front().DeltaPhi(pickedrJetLVec);
                          h1_dPhi_rJet_metVec.back()->Fill(dPhi_rJet_met, evtWeight * scaleMC);
                          h1_dPhi_rJet_bVec.back()->Fill(dPhi_rJet_b, evtWeight * scaleMC);
                          h2_dPhi_rJet_b_vs_metVec.back()->Fill(dPhi_rJet_met, dPhi_rJet_b, evtWeight * scaleMC);
                       }
   
                       double minDR_tJet_lept = 999.0; TLorentzVector pickedtJetLVec;
                       for(unsigned int it=0; it<topCombIdxVec.size(); it++){
                          int idx = topCombIdxVec.at(it);
                          double dR = pickedLeptLVec.DeltaR(jetsLVec_forTagger.at(idx));
                          if( minDR_tJet_lept > dR ){ minDR_tJet_lept = dR; pickedtJetLVec = jetsLVec_forTagger.at(idx); }
                       }

                       h1_cutFlowVec.back()->Fill("close_tJet_lept", (minDR_tJet_lept < 0.4) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_tJet_acc", (minDR_tJet_lept < 0.4) * (lostlept_cat == 0) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_tJet_rec", (minDR_tJet_lept < 0.4) * (lostlept_cat == 1) * evtWeight * scaleMC);
                       h1_cutFlowVec.back()->Fill("close_tJet_iso", (minDR_tJet_lept < 0.4) * (lostlept_cat == 2) * evtWeight * scaleMC);

                       double dPhi_top_met = metLVec.DeltaPhi(topLVec);
                       double minDphi_tJet_met = 999.0; TLorentzVector pickedtJetLVec_minDphi;
                       for(unsigned int it=0; it<topCombIdxVec.size(); it++){
                          int idx = topCombIdxVec.at(it);
                          double dPhi = metLVec.DeltaPhi(jetsLVec_forTagger.at(idx));
                          if( std::abs(minDphi_tJet_met) > std::abs(dPhi) ){ minDphi_tJet_met = dPhi; pickedtJetLVec_minDphi = jetsLVec_forTagger.at(idx); }
                       }
                       h1_dPhi_top_metVec.back()->Fill(dPhi_top_met, evtWeight * scaleMC);
                       h1_minDphi_tJet_metVec.back()->Fill(minDphi_tJet_met, evtWeight * scaleMC);
                    }
                    if( topDecayTypeVec[id] == 0 ){
/*
                       TLorentzVector genWdau1LVec = genDecayLVec.at(perDauIdxVec.at(1));
                       TLorentzVector genWdau2LVec = genDecayLVec.at(perDauIdxVec.at(2));
                       for(unsigned int idau = 0; idau< perDauIdxVec.size(); idau++){
                          int pdgId = genDecayPdgIdVec.at(perDauIdxVec.at(idau));
                          if( std::abs(pdgId) == 5 ){
                             pickedbLVec = perDauLVec.at(idau);
                             pickedbPdgId = pdgId;
                          }
                       }
                       bool ismatched = false;
                       for(unsigned int ib=0; ib<remainbJetsLVec.size(); ib++){
                          double dR = pickedbLVec.DeltaR(remainbJetsLVec.at(ib));
                          if( dR < 0.4 ) ismatched = true;
                       }
                       if( ismatched ){
                          h1_cutFlowVec.back()->Fill("rbJet_mtch_b", evtWeight * scaleMC);
                          if( rJets_cnt ==1 ){
                             h1_pt_genWdau1_nrJets_EQ1Vec.back()->Fill(genWdau1LVec.Pt(), evtWeight * scaleMC);
                             h1_eta_genWdau1_nrJets_EQ1Vec.back()->Fill(genWdau1LVec.Eta(), evtWeight * scaleMC);
                             h1_pt_genWdau2_nrJets_EQ1Vec.back()->Fill(genWdau2LVec.Pt(), evtWeight * scaleMC);
                             h1_eta_genWdau2_nrJets_EQ1Vec.back()->Fill(genWdau2LVec.Eta(), evtWeight * scaleMC);
                          }
                          if( rJets_cnt == 2 ){
                             TLorentzVector rJet = remainJetsLVec.front();
                             double dR_rJet_dau1 = rJet.DeltaR(genWdau1LVec);
                             double dR_rJet_dau2 = rJet.DeltaR(genWdau2LVec);
                             int is_rJet_matched = 0;
                             if( dR_rJet_dau1 < 0.4 ){
                                is_rJet_matched++;
                             }
                             if( dR_rJet_dau2 < 0.4 ){
                                is_rJet_matched++;
                             }
                             if( is_rJet_matched ==0 ) h1_cutFlowVec.back()->Fill("mtch0_rJets_EQ2", evtWeight * scaleMC);
                             if( is_rJet_matched ==1 ) h1_cutFlowVec.back()->Fill("mtch1_rJets_EQ2", evtWeight * scaleMC);
                             if( is_rJet_matched ==0 ){
                                TLorentzVector rbJet = remainbJetsLVec.front() + rJet;
                                double dR_rbJet = remainbJetsLVec.front().DeltaR(rJet);
                                double dPhi_rJet_met = metLVec.DeltaPhi(rJet);
                                h1_mass_rbJet_mtch0_nrJets_EQ2Vec.back()->Fill(rbJet.M(), evtWeight * scaleMC);
                                h1_dR_rbJet_mtch0_nrJets_EQ2Vec.back()->Fill(dR_rbJet, evtWeight * scaleMC);
                                h1_dPhi_rJet_met_mtch0_nrJets_EQ2Vec.back()->Fill(dPhi_rJet_met, evtWeight * scaleMC);
                             }
                             if( is_rJet_matched ==1 ){
                                TLorentzVector rbJet = remainbJetsLVec.front() + rJet;
                                double dR_rbJet = remainbJetsLVec.front().DeltaR(rJet);
                                double dPhi_rJet_met = metLVec.DeltaPhi(rJet);
                                h1_mass_rbJet_mtch1_nrJets_EQ2Vec.back()->Fill(rbJet.M(), evtWeight * scaleMC);
                                h1_dR_rbJet_mtch1_nrJets_EQ2Vec.back()->Fill(dR_rbJet, evtWeight * scaleMC);
                                h1_dPhi_rJet_met_mtch1_nrJets_EQ2Vec.back()->Fill(dPhi_rJet_met, evtWeight * scaleMC);
                             }
                          }
                          if( rJets_cnt >=3 ){
                             int mtch_genWdau_cnt = 0;
                             int isFirstMatched = 0;
                             for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                                if( genWdau1LVec.DeltaR(remainJetsLVec.at(ir)) < 0.4 ) isFirstMatched = 1;
                             }
                             int isSecondMatched = 0;
                             for(unsigned int ir=0; ir<remainJetsLVec.size(); ir++){
                                if( genWdau2LVec.DeltaR(remainJetsLVec.at(ir)) < 0.4 ) isSecondMatched = 1;
                             }
                             if( isFirstMatched ) mtch_genWdau_cnt++;
                             if( isSecondMatched ) mtch_genWdau_cnt++;
                             if( mtch_genWdau_cnt ==0 ) h1_cutFlowVec.back()->Fill("mtch0_genWdaus_LE3", evtWeight * scaleMC);
                             if( mtch_genWdau_cnt ==1 ) h1_cutFlowVec.back()->Fill("mtch1_genWdaus_LE3", evtWeight * scaleMC);
                             if( mtch_genWdau_cnt ==2 ) h1_cutFlowVec.back()->Fill("mtch2_genWdaus_LE3", evtWeight * scaleMC);
                          }
                       }
                       if( !ismatched ){
                          h1_cutFlowVec.back()->Fill("rbJet_NOT_mtch_b", evtWeight * scaleMC);
                          double dR_gen_reco_top = topLVec.DeltaR(pergenTopLVec);
                          h1_dR_gen_reco_topVec.back()->Fill(dR_gen_reco_top, evtWeight * scaleMC);
                       }
*/
                    }
                 }
              }
           }

           TLorentzVector sumLVec = -(topLVec + remainbJetsLVec[0] + metLVec);
           double dPhi_sum_MET = metLVec.DeltaPhi(sumLVec); 
           double MT_sum = type3Ptr->calcMT(sumLVec, metLVec);
           h1_dPhi_sum_METVec.back()->Fill(dPhi_sum_MET, evtWeight * scaleMC);
           h1_MT_sumVec.back()->Fill(MT_sum, evtWeight * scaleMC);

           h1_cutFlow_auxVec.back()->Fill("cntbJetLE1", (cntbJet >=1) * evtWeight * scaleMC);
           h1_cutFlow_auxVec.back()->Fill("pass_MTb", (best_bJet_MT >= 200) * evtWeight * scaleMC);
           h1_cutFlow_auxVec.back()->Fill("pass_MTrb", (best_lept_brJet_MT >= 200) * evtWeight * scaleMC);
           h1_cutFlow_auxVec.back()->Fill("aftMTcombMTcut", evtWeight * scaleMC * (best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("MTcombMT_MT2", evtWeight * scaleMC * (best_had_brJet_MT2 >= 300 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("+MET300", evtWeight * scaleMC * (metLVec.Pt()>300 && best_had_brJet_MT2 >= 300 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("+MET350", evtWeight * scaleMC * (metLVec.Pt()>350 && best_had_brJet_MT2 >= 300 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("+MET400", evtWeight * scaleMC * (metLVec.Pt()>400 && best_had_brJet_MT2 >= 300 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("+MET450", evtWeight * scaleMC * (metLVec.Pt()>450 && best_had_brJet_MT2 >= 300 && best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("pca_MTcombMT", evtWeight * scaleMC * (best_had_brJet_MT >=200 && best_aft_PCA_had_brJet_MT >= 0.3 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ) );
           h1_cutFlow_auxVec.back()->Fill("ori_mTcomb", evtWeight * scaleMC * (mTcomb >= AnaConsts::mTcombcut_) );
           h1_cutFlow_auxVec.back()->Fill("ori_MTcombMT", evtWeight * scaleMC * (mTcomb >= AnaConsts::mTcombcut_ && MT2 >= AnaConsts::MT2cut_) );
           h1_cutFlow_auxVec.back()->Fill("ori+MET300", evtWeight * scaleMC * (metLVec.Pt()>300 && mTcomb >= AnaConsts::mTcombcut_ && MT2 >= AnaConsts::MT2cut_) );
           h1_cutFlow_auxVec.back()->Fill("ori+MET350", evtWeight * scaleMC * (metLVec.Pt()>350 && mTcomb >= AnaConsts::mTcombcut_ && MT2 >= AnaConsts::MT2cut_) );
           h1_cutFlow_auxVec.back()->Fill("ori+MET400", evtWeight * scaleMC * (metLVec.Pt()>400 && mTcomb >= AnaConsts::mTcombcut_ && MT2 >= AnaConsts::MT2cut_) );
           h1_cutFlow_auxVec.back()->Fill("ori+MET450", evtWeight * scaleMC * (metLVec.Pt()>450 && mTcomb >= AnaConsts::mTcombcut_ && MT2 >= AnaConsts::MT2cut_) );

           if( cntbJet >= 1 ){
//              if( best_had_brJet_MT >= 200 && best_lept_brJet_MT >= 200
//               && ( ( best_lept_brJet_mTcomb == best_had_brJet_mTcomb && best_lept_brJet_mTcomb >= 600 ) || 
//               ( best_lept_brJet_mTcomb != best_had_brJet_mTcomb && best_lept_brJet_mTcomb > 500 && best_had_brJet_mTcomb > 800 ) ) 
//              ){
//              if( best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450) ){

              if( !ori_lept_brJetLVecMap.empty() ){
                 TLorentzVector ori_best_lept_brJet = ori_lept_brJetLVecMap.begin()->second;
                 double ori_best_lept_brJet_MT = type3Ptr->calcMT(ori_best_lept_brJet, metLVec);
                 h1_ori_MT_lept_brJetVec.back()->Fill(ori_best_lept_brJet_MT, evtWeight * scaleMC);
              }

              if( !ori_had_brJetLVecMap.empty() ){
                 TLorentzVector ori_best_had_brJet = ori_had_brJetLVecMap.begin()->second;
                 double ori_best_had_brJet_MT2 = type3Ptr->calcMT2(topLVec, ori_best_had_brJet, metLVec);
                 h1_ori_MT2_had_brJetVec.back()->Fill(ori_best_had_brJet_MT2, evtWeight * scaleMC);
              }

              if( !lept_brJetLVecMap.empty() ){
                 h1_MT_lept_brJetVec.back()->Fill(best_lept_brJet_MT, evtWeight * scaleMC);
                 h1_MT2_lept_brJetVec.back()->Fill(best_lept_brJet_MT2, evtWeight * scaleMC);
                 h2_MT_vs_MT2_lept_brJetVec.back()->Fill(best_lept_brJet_MT2, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_MT_lept_vs_top_brJetVec.back()->Fill(top_MT, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_mTcomb_vs_MT2_lept_brJetVec.back()->Fill(best_lept_brJet_MT2, best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 h2_mTcomb_vs_MT_lept_brJetVec.back()->Fill(best_lept_brJet_MT, best_lept_brJet_mTcomb, evtWeight * scaleMC);
              }
              if( !had_brJetLVecMap.empty() ){
                 h1_MT_had_brJetVec.back()->Fill(best_had_brJet_MT, evtWeight * scaleMC);
                 h1_MT2_had_brJetVec.back()->Fill(best_had_brJet_MT2, evtWeight * scaleMC);
                 h2_MT_vs_MT2_had_brJetVec.back()->Fill(best_had_brJet_MT2, best_had_brJet_MT, evtWeight * scaleMC);
                 h2_MT_had_vs_top_brJetVec.back()->Fill(top_MT, best_had_brJet_MT, evtWeight * scaleMC);
                 h2_mTcomb_vs_MT2_had_brJetVec.back()->Fill(best_had_brJet_MT2, best_had_brJet_mTcomb, evtWeight * scaleMC);
                 h2_mTcomb_vs_MT_had_brJetVec.back()->Fill(best_had_brJet_MT, best_had_brJet_mTcomb, evtWeight * scaleMC);
                 h2_aft_PCA_MT_had_vs_top_brJetVec.back()->Fill(best_aft_PCA_top_MT, best_aft_PCA_had_brJet_MT, evtWeight * scaleMC); 
                 h2_ori_MT2_vs_MT2_had_brJetVec.back()->Fill(best_had_brJet_MT2, MT2, evtWeight * scaleMC);
              }
              if( !lept_brJetLVecMap.empty() && !had_brJetLVecMap.empty() ){
                 if( best_had_brJet_MT2 != best_lept_brJet_MT2 ) h2_MT2_lept_vs_had_brJetVec.back()->Fill(best_had_brJet_MT2, best_lept_brJet_MT2, evtWeight * scaleMC);
                 else h1_MT2_same_lept_had_brJetVec.back()->Fill(best_lept_brJet_MT2, evtWeight * scaleMC);
                 if( best_had_brJet_MT != best_lept_brJet_MT ){
                    h2_MT_lept_vs_had_brJetVec.back()->Fill(best_had_brJet_MT, best_lept_brJet_MT, evtWeight * scaleMC);
                    h2_mTcomb_lept_vs_had_brJetVec.back()->Fill(best_had_brJet_mTcomb, best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 } else{
                    h1_MT_same_lept_had_brJetVec.back()->Fill(best_lept_brJet_MT, evtWeight * scaleMC);
                    h1_mTcomb_same_lept_had_brJetVec.back()->Fill(best_lept_brJet_mTcomb, evtWeight * scaleMC);
                 }
                 h2_MT2_lept_vs_MT_had_brJetVec.back()->Fill(best_had_brJet_MT, best_lept_brJet_MT2, evtWeight * scaleMC);
                 h2_MT_lept_vs_MT2_had_brJetVec.back()->Fill(best_had_brJet_MT2, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_MT_lept_vs_mTcomb_had_brJetVec.back()->Fill(best_had_brJet_mTcomb, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_MT_lept_vs_aft_PCA_had_brJetVec.back()->Fill(best_aft_PCA_had_brJet_MT, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_mTcomb_lept_vs_MT_had_brJetVec.back()->Fill(best_had_brJet_MT, best_lept_brJet_mTcomb, evtWeight * scaleMC);

                 h1_MT_bJetVec.back()->Fill(best_bJet_MT, evtWeight * scaleMC);
                 h2_MT_lept_vs_MT_bJetVec.back()->Fill(best_bJet_MT, best_lept_brJet_MT, evtWeight * scaleMC);
                 h2_MT_bJet_vs_mTcomb_had_brJetVec.back()->Fill(best_had_brJet_mTcomb, best_bJet_MT, evtWeight * scaleMC);
/*
                    h1_cutFlowVec.back()->Fill("aftMTcut", evtWeight * scaleMC);
                    if( nbJets == 1 ) h1_cutFlowVec.back()->Fill("_nb1", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("_nb1", 0);
                    if( nbJets == 2 ) h1_cutFlowVec.back()->Fill("_nb2", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("_nb2", 0);
                    if( nbJets >= 3 ) h1_cutFlowVec.back()->Fill("_nb3", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("_nb3", 0);
                    if( best_had_brJet_MT2 != best_lept_brJet_MT2 ) h2_MT2_lept_vs_had_aft_MT_cuts_brJetVec.back()->Fill(best_had_brJet_MT2, best_lept_brJet_MT2, evtWeight * scaleMC);
                    else h1_MT2_same_lept_had_aft_MT_cuts_brJetsVec.back()->Fill(best_lept_brJet_MT2, evtWeight * scaleMC);

                    if( ( best_lept_brJet_MT2 == best_had_brJet_MT2 && best_lept_brJet_MT2 >= 300) || 
                        ( best_lept_brJet_MT2 != best_had_brJet_MT2 && (best_lept_brJet_MT2 >= 300 || best_had_brJet_MT2 >= 300) ) ){
                       h1_cutFlowVec.back()->Fill("aftMTMT2cuts", evtWeight * scaleMC);
                       if( nbJets == 1 ) h1_cutFlowVec.back()->Fill("__nb1", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("__nb1", 0);
                       if( nbJets == 2 ) h1_cutFlowVec.back()->Fill("__nb2", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("__nb2", 0);
                       if( nbJets >= 3 ) h1_cutFlowVec.back()->Fill("__nb3", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("__nb3", 0);
                    }

                    if( ( best_lept_brJet_mTcomb == best_had_brJet_mTcomb && best_lept_brJet_mTcomb >= 700 ) || 
                        ( best_lept_brJet_mTcomb != best_had_brJet_mTcomb && (best_lept_brJet_mTcomb > 700 || best_had_brJet_mTcomb > 700) ) ){
                       h1_cutFlowVec.back()->Fill("aftMTcombcut", evtWeight * scaleMC);
                       if( nbJets == 1 ) h1_cutFlowVec.back()->Fill("___nb1", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("___nb1", 0);
                       if( nbJets == 2 ) h1_cutFlowVec.back()->Fill("___nb2", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("___nb2", 0);
                       if( nbJets >= 3 ) h1_cutFlowVec.back()->Fill("___nb3", evtWeight * scaleMC); else h1_cutFlowVec.back()->Fill("___nb3", 0);
                    }
*/
              }else{
                 std::cout<<"WARNING ... NOT empty?! lept_brJetLVecMap.size : "<<lept_brJetLVecMap.size()<<"  had_brJetLVecMap.size : "<<had_brJetLVecMap.size()<<std::endl;
              }
//              } // MT cut
              h1_nComb_lept_brJetVec.back()->Fill(lept_brJetIdxMap.size(), evtWeight * scaleMC);
              h1_nComb_had_brJetVec.back()->Fill(had_brJetIdxMap.size(), evtWeight * scaleMC);
              h2_nComb_lept_vs_had_brJetVec.back()->Fill(had_brJetIdxMap.size(), lept_brJetIdxMap.size(), evtWeight * scaleMC);
           }
        }

        bool passMore = true;
/*
        if( nTops ==1 // in a specific nTop search bin
           && doMT2mTcombCuts && (MT2 <= AnaConsts::MT2cut_ || mTcomb <= AnaConsts::mTcombcut_) // if doMT2mTcombCuts is true and either MT2 or mTcomb not pass 
        ) passMore = false;
*/
/*
        if( nTops ==1 // in a specific nTop search bin
           && doMT2mTcombCuts && !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450))
        ) passMore = false;
*/
        if( (nTops ==1 || nTops ==2 ) // in a specific nTop search bin
           && doMT2mTcombCuts && !(best_had_brJet_MT >=200 && best_had_brJet_mTcomb >= 400 && (best_lept_brJet_MT >= 200 || best_had_brJet_mTcomb >= 450))
        ) passMore = false;

        if( passMore ){ cnt_passMore_WeightedScaledMC += evtWeight * scaleMC; cnt_passMore_WeightedErrorScaledMC += evtWeight * evtWeight * scaleMC * scaleMC; }

// lepton always means e or mu; AND only had. tau is separated out as of now
        int cntlepton =0, cnttauHad =0;
        std::vector<int> idxWjetsVec, idxTopjetsVec;
        std::vector<int> valCSVSWjetsVec, valCSVSTopjetsVec;
        std::vector<int> bTagInfoWjetsVec, bTagInfoTopjetsVec;

        for(unsigned int ij=0; ij<jetsLVec_forTagger.size(); ij++){
           TLorentzVector perJetLVec = jetsLVec_forTagger.at(ij);
           double perCSVS = recoJetsBtag_forTagger.at(ij);

           std::vector<TLorentzVector> dummyLVec; std::vector<double> dummyBtagVec;
           dummyLVec.push_back(perJetLVec); dummyBtagVec.push_back(perCSVS);

           if( !AnaFunctions::countJets(dummyLVec, AnaConsts::pt30Eta24Arr) ) continue;

           int isBtagged = AnaFunctions::countCSVS(dummyLVec, dummyBtagVec, AnaConsts::cutCSVS, AnaConsts::bTagArr);

           if( perJetLVec.M() >= AnaConsts::lowWjetMass_ && perJetLVec.M() <= AnaConsts::highWjetMass_ ){
              double minDR_otherbJets = 999.0;
              for(unsigned int jj=0; jj<jetsLVec_forTagger.size(); jj++){
                 if( jj == ij ) continue;
                 dummyLVec.clear(); dummyBtagVec.clear();
                 dummyLVec.push_back(jetsLVec_forTagger.at(jj)); dummyBtagVec.push_back(recoJetsBtag_forTagger.at(jj));
                 if( AnaFunctions::countCSVS(dummyLVec, dummyBtagVec, AnaConsts::cutCSVS, AnaConsts::bTagArr) ){
                    double deltaR_otherbJets = perJetLVec.DeltaR(jetsLVec_forTagger.at(jj));
                    if( minDR_otherbJets > deltaR_otherbJets ) minDR_otherbJets = deltaR_otherbJets;
                 }
              }
              idxWjetsVec.push_back(ij);
              valCSVSWjetsVec.push_back(perCSVS);
              bTagInfoWjetsVec.push_back(isBtagged);
              h1_Wjets_bTagCatsVec.back()->Fill(isBtagged, evtWeight*scaleMC);
              double minDR_Wjet_genb = 999.0, minDR_Wjet_genW = 999.0;
              int minDR_genb_idx = -1, minDR_genW_idx = -1;
              std::vector<int> allDRless0p4_genb_idxVec, allDRless0p4_genW_idxVec;
              for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
                 const int pdgId = genDecayPdgIdVec.at(ig);
                 if( std::abs(pdgId) != 5 && std::abs(pdgId) != 24 ) continue;
                 double perDR = perJetLVec.DeltaR(genDecayLVec.at(ig));
                 if( std::abs(pdgId) == 5 ){
                    if( minDR_Wjet_genb > perDR ){ minDR_Wjet_genb = perDR; minDR_genb_idx = ig; } 
                    if( perDR < 0.5 ) allDRless0p4_genb_idxVec.push_back(ig);
                 }
                 if( std::abs(pdgId) == 24 ){
                    if( minDR_Wjet_genW > perDR ){ minDR_Wjet_genW = perDR; minDR_genW_idx = ig; } 
                    if( perDR < 0.5 ) allDRless0p4_genW_idxVec.push_back(ig);
                 }
              }
              if( isBtagged ){
                 h2_bTagged_Wjets_minDR_genb_vs_minDR_genWVec.back()->Fill(minDR_Wjet_genW, minDR_Wjet_genb, evtWeight * scaleMC);
                 h2_bTagged_Wjets_mass_vs_minDR_genbVec.back()->Fill(minDR_Wjet_genb, perJetLVec.M(), evtWeight * scaleMC);
                 h1_bTagged_Wjets_CSVVec.back()->Fill(perCSVS, evtWeight * scaleMC);
                 h1_bTagged_Wjets_massVec.back()->Fill(perJetLVec.M(), evtWeight * scaleMC);
                 if( minDR_Wjet_genb < 0.5 && minDR_Wjet_genW < 0.5 ){
                    bool printDecay = false;
                    h1_bTagged_Wjets_minDR_otherbJetsVec.back()->Fill(minDR_otherbJets, evtWeight * scaleMC);

                    int pickedMomIdx = -1; std::vector<int> pickedDauIdxVec;
                    for(unsigned int jg=0; jg<genDecayPdgIdVec.size(); jg++){
                       if( genDecayIdxVec.at(jg) == genDecayMomIdxVec[minDR_genb_idx] ) pickedMomIdx = jg;
                       if( genDecayMomIdxVec.at(jg) == genDecayIdxVec[minDR_genW_idx] ) pickedDauIdxVec.push_back(jg);
                    }

                    int idxDau1 = -1, idxDau2 = -1, idxDaub = minDR_genb_idx; 
                    double deltaR1b_Wjet_genDaus = -1, deltaR2b_Wjet_genDaus = -1, deltaR12_Wjet_genDaus = -1;
                    if( pickedDauIdxVec.size()!= 2 ){ /*std::cout<<"WARNING ... "<<pickedDauIdxVec.size()<<"  of decay product of a top quark with idx : "<<genDecayIdxVec.at(pickedMomIdx)<<"!"<<std::endl; printDecay = true; */ }
                    else{
                       idxDau1 = pickedDauIdxVec[0]; idxDau2 = pickedDauIdxVec[1];
                       TLorentzVector dau1LVec = genDecayLVec.at(idxDau1), dau2LVec = genDecayLVec.at(idxDau2), daubLVec = genDecayLVec.at(idxDaub);
                       if( dau1LVec.Pt() < dau2LVec.Pt() ){
                          int idxtmp = idxDau1; idxDau1 = idxDau2; idxDau2 = idxtmp;
                          TLorentzVector tmpdau1LVec = dau1LVec; dau1LVec = dau2LVec; dau2LVec = tmpdau1LVec;
                       }
                       deltaR1b_Wjet_genDaus = dau1LVec.DeltaR(daubLVec);
                       deltaR2b_Wjet_genDaus = dau2LVec.DeltaR(daubLVec);
                       deltaR12_Wjet_genDaus = dau1LVec.DeltaR(dau2LVec);
                    }

                    if( genDecayMomIdxVec[minDR_genb_idx] != genDecayMomIdxVec[minDR_genW_idx] ){
                       bool ifAnySameMother = false;
                       for(unsigned int ib=0; ib<allDRless0p4_genb_idxVec.size(); ib++){
                          for(unsigned int iw=0; iw<allDRless0p4_genW_idxVec.size(); iw++){
                             if( genDecayMomIdxVec[allDRless0p4_genb_idxVec[ib]] == genDecayMomIdxVec[allDRless0p4_genW_idxVec[iw]] ) ifAnySameMother = true;
                          }
                       }
//                       std::cout<<"WARNING ... not same mother?! for genb mom : "<<genDecayMomIdxVec[minDR_genb_idx]<<"  genW mom : "<<genDecayMomIdxVec[minDR_genW_idx]<<"  allDRless0p4_genb_idxVec.size : "<<allDRless0p4_genb_idxVec.size()<<"  allDRless0p4_genW_idxVec.size : "<<allDRless0p4_genW_idxVec.size()<<"  ifAnySameMother : "<<ifAnySameMother<<std::endl; 
//                       printDecay = true;

                       h1_bTagged_rndmComb_Wjets_ptVec.back()->Fill(perJetLVec.Pt(), evtWeight * scaleMC);

                    }else{
                       h1_bTagged_sameTop_Wjets_ptVec.back()->Fill(perJetLVec.Pt(), evtWeight * scaleMC);

                       if( pickedMomIdx == -1 ){ std::cout<<"ERROR ... mom idx is not stored for the "<<minDR_genb_idx<<"th gen particle!"<<std::endl; printDecay = true; }
                       else if( std::abs(genDecayPdgIdVec.at(pickedMomIdx)) == 6 ){
                          if(pickedDauIdxVec.size()== 2 ){
                              h1_bTagged_Wjets_deltaR1b_genDausVec.back()->Fill(deltaR1b_Wjet_genDaus, evtWeight * scaleMC); 
                              h1_bTagged_Wjets_deltaR2b_genDausVec.back()->Fill(deltaR2b_Wjet_genDaus, evtWeight * scaleMC); 
                              h1_bTagged_Wjets_deltaR12_genDausVec.back()->Fill(deltaR12_Wjet_genDaus, evtWeight * scaleMC); 
                          }
                       }
                    }
                    if( printDecay ) std::cout<<"INFO ... "<<tr->getVec<std::string>("genDecayStrVec").front().c_str()<<std::endl<<std::endl;
                 }
              }else{
                 h2_nobTagged_Wjets_minDR_genb_vs_minDR_genWVec.back()->Fill(minDR_Wjet_genW, minDR_Wjet_genb, evtWeight * scaleMC);
                 h1_nobTagged_Wjets_massVec.back()->Fill(perJetLVec.M(), evtWeight * scaleMC);
              }
           } 
           if( perJetLVec.M() >= AnaConsts::lowTopjetMass_ && perJetLVec.M() <= AnaConsts::highTopjetMass_ ){
              idxTopjetsVec.push_back(ij);
              valCSVSTopjetsVec.push_back(recoJetsBtag_forTagger.at(ij));
              bTagInfoTopjetsVec.push_back(isBtagged);
              h1_Topjets_bTagCatsVec.back()->Fill(isBtagged, evtWeight*scaleMC);
              double minDR_Topjet_genb = 999.0, minDR_Topjet_genTop = 999.0;
              for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
                 const int pdgId = genDecayPdgIdVec.at(ig);
                 if( std::abs(pdgId) != 5 && std::abs(pdgId) != 6 ) continue;
                 double perDR = perJetLVec.DeltaR(genDecayLVec.at(ig));
                 if( std::abs(pdgId) == 5 ){ if( minDR_Topjet_genb > perDR ) minDR_Topjet_genb = perDR; }
                 if( std::abs(pdgId) == 6 ){ if( minDR_Topjet_genTop > perDR ) minDR_Topjet_genTop = perDR; }
              }
              if( isBtagged ){
                 h2_bTagged_Topjets_minDR_genb_vs_minDR_genTopVec.back()->Fill(minDR_Topjet_genTop, minDR_Topjet_genb, evtWeight * scaleMC);
                 h1_bTagged_Topjets_massVec.back()->Fill(perJetLVec.M(), evtWeight * scaleMC);
              }else{
                 h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTopVec.back()->Fill(minDR_Topjet_genTop, minDR_Topjet_genb, evtWeight * scaleMC);
                 h1_nobTagged_Topjets_massVec.back()->Fill(perJetLVec.M(), evtWeight * scaleMC);
              }
           } 
        }

        cntlepton = W_emuVec.size() + W_tau_emuVec.size();
        for(unsigned int ig=0; ig<genDecayPdgIdVec.size(); ig++){
           int pdgId = genDecayPdgIdVec.at(ig);
           if( abs(pdgId) == 15 ){
              bool hasprongs = false;
              for(unsigned int it=0; it<W_tau_prongsVec.size(); it++){
                 if( find_mother(ig, W_tau_prongsVec.at(it), genDecayIdxVec, genDecayMomIdxVec) ) hasprongs = true;
              }
              if( hasprongs) cnttauHad ++;
           }
        }
        TLorentzVector pickedMaxPtLeptORprongLVec;
        for(unsigned int ip=0; ip<W_emuVec.size(); ip++){
           int idx = W_emuVec.at(ip);
           if( pickedMaxPtLeptORprongLVec.Pt() < genDecayLVec.at(idx).Pt() ) pickedMaxPtLeptORprongLVec = genDecayLVec.at(idx);
        }

        for(unsigned int ip=0; ip<W_tau_emuVec.size(); ip++){
           int idx = W_tau_emuVec.at(ip);
           if( pickedMaxPtLeptORprongLVec.Pt() < genDecayLVec.at(idx).Pt() ) pickedMaxPtLeptORprongLVec = genDecayLVec.at(idx);
        }

        for(unsigned int ip=0; ip<W_tau_prongsVec.size(); ip++){
           int idx = W_tau_prongsVec.at(ip);
           if( pickedMaxPtLeptORprongLVec.Pt() < genDecayLVec.at(idx).Pt() ) pickedMaxPtLeptORprongLVec = genDecayLVec.at(idx);
        }

        if( !keyStringT.Contains("Signal") ) if( cntlepton + cnttauHad > 2 ) std::cout<<"WARNING ... cntlepton : "<<cntlepton<<"  cnttauHad : "<<cnttauHad<<std::endl;

        if( passMore ){

           for(int iSR=0; iSR<nSR; iSR++){
              if(    (nbJets_SR_lo[iSR] == -1 || nbJets >= nbJets_SR_lo[iSR]) && (nbJets_SR_hi[iSR] == -1 || nbJets <= nbJets_SR_hi[iSR])
                  && (nTops_SR_lo[iSR] == -1 || nTops >= nTops_SR_lo[iSR]) && (nTops_SR_hi[iSR] == -1 || nTops <= nTops_SR_hi[iSR])
                ){
   
                 int nJetsRsys = 0;
                 for(unsigned int ij=0; ij<jetsLVec_forTagger.size(); ij++){
                    std::vector<TLorentzVector> dummyLVec; dummyLVec.push_back(jetsLVec_forTagger.at(ij));
                    std::vector<double> dummyCSVSVec; dummyCSVSVec.push_back(recoJetsBtag_forTagger.at(ij));
                    if( AnaFunctions::countCSVS(dummyLVec, dummyCSVSVec, AnaConsts::cutCSVS, AnaConsts::bTagArr) ) continue;
                    if( !AnaFunctions::countJets(dummyLVec, AnaConsts::pt30Eta24Arr) ) continue;
                    
                    if( std::find(type3Ptr->allIdxCacheVec.begin(), type3Ptr->allIdxCacheVec.end(), ij) != type3Ptr->allIdxCacheVec.end() ) continue;

                    nJetsRsys++;
                 }

                 h1_nJetsVec[iSR].back()->Fill(tr->getVar<int>("cntNJetsPt30Eta24"), evtWeight*scaleMC);
                 h1_nJetsRsysVec[iSR].back()->Fill(nJetsRsys, evtWeight*scaleMC);
                 h1_metVec[iSR].back()->Fill(tr->getVar<double>("met"), evtWeight*scaleMC);
                 if( nTops <=1 ){
//                    h1_MT2Vec[iSR].back()->Fill(tr->getVar<double>("MT2"), evtWeight*scaleMC);
                    h1_MT2Vec[iSR].back()->Fill(best_had_brJet_MT2, evtWeight*scaleMC);
                    h2_MT2_vs_metVec[iSR].back()->Fill(tr->getVar<double>("met"), best_had_brJet_MT2, evtWeight*scaleMC);
                    h2_coarse_bin_MT2_vs_metVec[iSR].back()->Fill(tr->getVar<double>("met"), best_had_brJet_MT2, evtWeight*scaleMC);
                 }else{
                    const int combIdx0 = type3Ptr->pickedTopCandSortedVec[0];
                    const std::vector<int> topCombIdxVec0 = type3Ptr->finalCombfatJets[combIdx0];
                    const TLorentzVector topLVec0 = type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec0);
                   
                    const int combIdx1 = type3Ptr->pickedTopCandSortedVec[1];
                    const std::vector<int> topCombIdxVec1 = type3Ptr->finalCombfatJets[combIdx1];
                    const TLorentzVector topLVec1 = type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec1);
                   
                    double altMT2 = type3Ptr->calcMT2(topLVec0, topLVec1, metLVec);
                    h1_MT2Vec[iSR].back()->Fill(altMT2, evtWeight*scaleMC);
                    h2_MT2_vs_metVec[iSR].back()->Fill(tr->getVar<double>("met"), altMT2, evtWeight*scaleMC);
                    h2_coarse_bin_MT2_vs_metVec[iSR].back()->Fill(tr->getVar<double>("met"), altMT2, evtWeight*scaleMC);
                 }
                 h1_mTcombVec[iSR].back()->Fill(tr->getVar<double>("mTcomb"), evtWeight*scaleMC);
                 h1_HTVec[iSR].back()->Fill(tr->getVar<double>("HT"), evtWeight*scaleMC);
                 h2_met_vs_nJetsVec[iSR].back()->Fill(tr->getVar<int>("cntNJetsPt30Eta24"), tr->getVar<double>("met"), evtWeight*scaleMC);

                 for(int iSR_met = 0; iSR_met < nSR_met; iSR_met ++ ){
                    if(    (met_SR_lo[iSR_met] == -1 || tr->getVar<double>("met") >= met_SR_lo[iSR_met] )
                        && (met_SR_hi[iSR_met] == -1 || tr->getVar<double>("met") < met_SR_hi[iSR_met] )
                      ){
                       cnt_passSRmet_WeightedScaledMCVec.back()[iSR][iSR_met] += evtWeight*scaleMC;
                       cnt_passSRmet_WeightedErrorScaledMCVec.back()[iSR][iSR_met] += pow(evtWeight*scaleMC, 2.0);
                       if( !keyStringT.Contains("Signal") ){
                          cnt_passSRmet_sumSM_WeightedScaledMCVec[iSR][iSR_met] += evtWeight*scaleMC;
                          cnt_passSRmet_sumSM_WeightedErrorScaledMCVec[iSR][iSR_met] += pow(evtWeight*scaleMC, 2.0);
                       }
                    }
                 }
              }
           }

           if( cntlepton == 0 && cnttauHad == 0 ){ cnt_allHad_WeightedScaledMC += evtWeight * scaleMC; cnt_allHad_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }
           if( cntlepton == 0 && cnttauHad == 1 ){ cnt_Wtauhad_WeightedScaledMC += evtWeight * scaleMC; cnt_Wtauhad_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }
           if( cntlepton == 0 && cnttauHad >= 2 ){ cnt_diWtauhads_WeightedScaledMC += evtWeight * scaleMC; cnt_diWtauhads_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }
           if( cntlepton == 1 && cnttauHad == 0 ){ cnt_onelepton_WeightedScaledMC += evtWeight * scaleMC; cnt_onelepton_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }
           if( cntlepton == 1 && cnttauHad >= 1 ){ cnt_dileptons_inc_Wtauhad_WeightedScaledMC += evtWeight * scaleMC; cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }
           if( cntlepton >= 2 && cnttauHad == 0 ){ cnt_dileptons_WeightedScaledMC += evtWeight * scaleMC; cnt_dileptons_WeightedErrorScaledMC += pow(evtWeight * scaleMC, 2.0); }

           std::vector<TLorentzVector> topCandLVec, bestWLVec, bestbLVec, bestLeptLVec;
           std::vector<int> bestbJetIdxVec;
           std::vector<double> mTtopVec, mTWVec, mTbVec;
           std::vector<double> minDphiLeptMETVec;
           std::vector<double> bestW_dausDRVec;
           std::vector<double> deltaR12Vec, deltaR1bVec, deltaR2bVec;
           std::vector<double> relPt12Vec, relPt1bVec, relPt2bVec;

           for(unsigned int ip=0; ip<type3Ptr->pickedTopCandSortedVec.size(); ip++){
              const int combIdx = type3Ptr->pickedTopCandSortedVec[ip];
              const std::vector<int> topCombIdxVec = type3Ptr->finalCombfatJets[combIdx];
              const TLorentzVector topLVec = type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec);

              double bestWmass = -1; int bestbJetIdx = -1; TLorentzVector perbestWLVec, perbestbLVec, perLeptLVec;
              int pickedLeptIdx = -1; double minDphiLeptMET = -1;
              double perbestW_dausDR = -1;
              double perdeltaR12 = -1, perdeltaR1b = -1, perdeltaR2b = -1;
              double perrelPt12 = -1, perrelPt1b = -1, perrelPt2b = -1;
              for(unsigned int ic=0; ic<topCombIdxVec.size(); ic++){
                 const int idxJet = topCombIdxVec.at(ic);
                 const TLorentzVector perJetLVec = jetsLVec_forTagger.at(idxJet);
                 std::vector<int> cachedIdxVec = topCombIdxVec; cachedIdxVec.erase(cachedIdxVec.begin()+ic);
                 const TLorentzVector perDiJetsLVec = type3Ptr->buildLVec(jetsLVec_forTagger, cachedIdxVec);

                 const double dPhiJetMET = metLVec.DeltaPhi(perJetLVec);
                 
                 std::vector<TLorentzVector> dummyLVec; std::vector<double> dummyBtagVec;
                 dummyLVec.push_back(perJetLVec); dummyBtagVec.push_back(recoJetsBtag_forTagger.at(idxJet));
                 int isBtagged = AnaFunctions::countCSVS(dummyLVec, dummyBtagVec, AnaConsts::cutCSVS, AnaConsts::bTagArr);

                 if( topCombIdxVec.size() == 1 ){
// One ak4 jet is a top quark cand. Set a special value of W mass (0)
// Note that bestbJetIdx could be -1 and perbestWLVec, perbestbLVec could be default (0) TLreontzVector!
                    bestWmass = mW_; if( isBtagged ){ bestbJetIdx = idxJet; perbestbLVec = perJetLVec; }
                    perbestWLVec.SetVectM(perJetLVec.Vect(), mW_);
                    perbestW_dausDR = -1;
//                    perdeltaR12 = 0.4; perdeltaR1b = 0.4; perdeltaR2b = 0.4; 
                    perdeltaR12 = -1.0; perdeltaR1b = -1.0; perdeltaR2b = -1.0; 
                    perrelPt12 = -1.0; perrelPt1b = -1.0; perrelPt2b = -1.0; 
                    break;
                 }
                 if( isBtagged ){
// If in triplet or doublet, there is a btagged jet then it's done: W would be the combination without b jet!
// Note that no >=2 b jets cases in top candidates!
                    bestbJetIdx = idxJet;
                    bestWmass = perDiJetsLVec.M();
                    perbestWLVec = perDiJetsLVec;
                    perbestbLVec = perJetLVec;
                    if( cachedIdxVec.size() == 1 ){ 
                       perbestW_dausDR = -1;
//                       perdeltaR12 = 0.4; 
                       perdeltaR12 = -1.0; 
                       perdeltaR1b = perdeltaR2b = perbestbLVec.DeltaR(perbestWLVec);
                       perrelPt12 = -1.0; perrelPt1b = -1.0; perrelPt2b = -1.0;
                    }else{
                       perbestW_dausDR = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(jetsLVec_forTagger.at(cachedIdxVec[1]));
                       perdeltaR12 = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(jetsLVec_forTagger.at(cachedIdxVec[1]));
                       perdeltaR1b = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(perbestbLVec);
                       perdeltaR2b = jetsLVec_forTagger.at(cachedIdxVec[1]).DeltaR(perbestbLVec);
                       perrelPt12 = jetsLVec_forTagger.at(cachedIdxVec[1]).Pt()/jetsLVec_forTagger.at(cachedIdxVec[0]).Pt();
                       perrelPt1b = jetsLVec_forTagger.at(cachedIdxVec[0]).Pt()/perbestbLVec.Pt();
                       perrelPt2b = jetsLVec_forTagger.at(cachedIdxVec[1]).Pt()/perbestbLVec.Pt();
                    }
                    break;
                 }
// If other cases, find the best W mass closer to mW_
                 if( bestWmass == -1 || std::abs(bestWmass - mW_) > std::abs(perDiJetsLVec.M() - mW_) ){
                    bestbJetIdx = idxJet; 
                    bestWmass = perDiJetsLVec.M(); 
                    perbestWLVec = perDiJetsLVec;
                    perbestbLVec = perJetLVec;
                    if( cachedIdxVec.size() == 1 ){
                       perbestW_dausDR = -1;
//                       perdeltaR12 = 0.4;
                       perdeltaR12 = -1.0;
                       perdeltaR1b = perdeltaR2b = perbestbLVec.DeltaR(perbestWLVec);
                       perrelPt12 = -1.0; perrelPt1b = -1.0; perrelPt2b = -1.0;
                    }else{
                       perbestW_dausDR = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(jetsLVec_forTagger.at(cachedIdxVec[1]));
                       perdeltaR12 = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(jetsLVec_forTagger.at(cachedIdxVec[1]));
                       perdeltaR1b = jetsLVec_forTagger.at(cachedIdxVec[0]).DeltaR(perbestbLVec);
                       perdeltaR2b = jetsLVec_forTagger.at(cachedIdxVec[1]).DeltaR(perbestbLVec);
                       perrelPt12 = jetsLVec_forTagger.at(cachedIdxVec[1]).Pt()/jetsLVec_forTagger.at(cachedIdxVec[0]).Pt();
                       perrelPt1b = jetsLVec_forTagger.at(cachedIdxVec[0]).Pt()/perbestbLVec.Pt();
                       perrelPt2b = jetsLVec_forTagger.at(cachedIdxVec[1]).Pt()/perbestbLVec.Pt();
                    }
                 }
                 if( !isBtagged ){
                    if( pickedLeptIdx == -1 || (minDphiLeptMET > std::abs(dPhiJetMET) ) ){ pickedLeptIdx = idxJet; minDphiLeptMET = std::abs(dPhiJetMET); perLeptLVec = perJetLVec; }
                 }
              }
              double mTtop = calcMT(topLVec, metLVec);
              double mTW = calcMT(perbestWLVec, metLVec);
              double mTb = calcMT(perbestbLVec, metLVec);
              topCandLVec.push_back(topLVec); bestWLVec.push_back(perbestWLVec); bestbLVec.push_back(perbestbLVec);
              bestbJetIdxVec.push_back(bestbJetIdx);
              mTtopVec.push_back(mTtop); mTWVec.push_back(mTW); mTbVec.push_back(mTb);
              if( pickedLeptIdx == -1 ){ perLeptLVec = topLVec; minDphiLeptMET = std::abs(topLVec.DeltaPhi(metLVec)); }
//              if( pickedLeptIdx == bestbJetIdx ){ perLeptLVec.SetPxPyPzE(0, 0, 0, 0); }
              bestLeptLVec.push_back(perLeptLVec); minDphiLeptMETVec.push_back(minDphiLeptMET);
              bestW_dausDRVec.push_back(perbestW_dausDR);
              deltaR12Vec.push_back(perdeltaR12); deltaR1bVec.push_back(perdeltaR1b); deltaR2bVec.push_back(perdeltaR2b);
              relPt12Vec.push_back(perrelPt12); relPt1bVec.push_back(perrelPt1b); relPt2bVec.push_back(perrelPt2b);
           }

           double MT2_nTopsEQ2 = 0, MT2_TopAndbLept_nTopsEQ2 = 0;
           if( nTops == 2 ){
              MT2_nTopsEQ2 = type3Ptr->calcMT2(topCandLVec[0], topCandLVec[1], metLVec);
              TLorentzVector bPlusLeptLVec = bestbLVec[1] + bestLeptLVec[1]; 
              MT2_TopAndbLept_nTopsEQ2 = type3Ptr->calcMT2(topCandLVec[0], bPlusLeptLVec, metLVec); 
              if( minDphiLeptMETVec[1] > 1.0 ) MT2_TopAndbLept_nTopsEQ2 = MT2_nTopsEQ2;
           }

           for(unsigned int ic=0; ic<topCandLVec.size(); ic++){
              if( ic+1 > nTopCandToPlot ) break;

              const int combIdx = type3Ptr->ori_pickedTopCandSortedVec[ic];
              const std::vector<int> topCombIdxVec = type3Ptr->finalCombfatJets[combIdx];
              h1_topCand_nJetsVec[ic].back()->Fill(topCombIdxVec.size(), evtWeight*scaleMC);

              h1_topCand_MVec[ic].back()->Fill(type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec).M(), evtWeight*scaleMC);
              h1_topCand_PtVec[ic].back()->Fill(type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec).Pt(), evtWeight*scaleMC);
              h1_topCand_EtaVec[ic].back()->Fill(type3Ptr->buildLVec(jetsLVec_forTagger, topCombIdxVec).Eta(), evtWeight*scaleMC);

//              h1_topCand_MVec[ic].back()->Fill(topCandLVec[ic].M(), evtWeight*scaleMC);
//              h1_topCand_PtVec[ic].back()->Fill(topCandLVec[ic].Pt(), evtWeight*scaleMC);
//              h1_topCand_EtaVec[ic].back()->Fill(topCandLVec[ic].Eta(), evtWeight*scaleMC);

              h1_W_MVec[ic].back()->Fill(bestWLVec[ic].M(), evtWeight*scaleMC);
              h1_W_PtVec[ic].back()->Fill(bestWLVec[ic].Pt(), evtWeight*scaleMC);
              h1_W_EtaVec[ic].back()->Fill(bestWLVec[ic].Eta(), evtWeight*scaleMC);
              h1_W_dausDRVec[ic].back()->Fill(bestW_dausDRVec[ic], evtWeight*scaleMC);
              h2_W_dausDR_versus_W_MVec[ic].back()->Fill(bestWLVec[ic].M(), bestW_dausDRVec[ic], evtWeight*scaleMC);

              h1_b_MVec[ic].back()->Fill(bestbLVec[ic].M(), evtWeight*scaleMC);
              h1_b_PtVec[ic].back()->Fill(bestbLVec[ic].Pt(), evtWeight*scaleMC);
              h1_b_EtaVec[ic].back()->Fill(bestbLVec[ic].Eta(), evtWeight*scaleMC);

              h1_ratio_mW_over_mTopVec[ic].back()->Fill(bestWLVec[ic].M()/topCandLVec[ic].M(), evtWeight*scaleMC);

              h2_mW_versus_mTopVec[ic].back()->Fill(topCandLVec[ic].M(), bestWLVec[ic].M(), evtWeight*scaleMC);
              h2_ratio_mW_over_mTop_versus_mTopVec[ic].back()->Fill(topCandLVec[ic].M(), bestWLVec[ic].M()/topCandLVec[ic].M(), evtWeight*scaleMC);

              h1_mTtopVec[ic].back()->Fill(mTtopVec[ic], evtWeight*scaleMC);
              h1_mTWVec[ic].back()->Fill(mTWVec[ic], evtWeight*scaleMC);
              h1_mTbVec[ic].back()->Fill(mTbVec[ic], evtWeight*scaleMC);

              h1_deltaR12Vec[ic].back()->Fill(deltaR12Vec[ic], evtWeight*scaleMC);
              h1_deltaR1bVec[ic].back()->Fill(deltaR1bVec[ic], evtWeight*scaleMC);
              h1_deltaR2bVec[ic].back()->Fill(deltaR2bVec[ic], evtWeight*scaleMC);
              h2_deltaR12_vs_deltaR1bVec[ic].back()->Fill(deltaR1bVec[ic], deltaR12Vec[ic], evtWeight*scaleMC);
              h2_deltaR12_vs_deltaR2bVec[ic].back()->Fill(deltaR2bVec[ic], deltaR12Vec[ic], evtWeight*scaleMC);
              h2_deltaR1b_vs_deltaR2bVec[ic].back()->Fill(deltaR2bVec[ic], deltaR1bVec[ic], evtWeight*scaleMC);

              h1_relPt12Vec[ic].back()->Fill(relPt12Vec[ic], evtWeight*scaleMC);
              h1_relPt1bVec[ic].back()->Fill(relPt1bVec[ic], evtWeight*scaleMC);
              h1_relPt2bVec[ic].back()->Fill(relPt2bVec[ic], evtWeight*scaleMC);
              h2_relPt12_vs_relPt1bVec[ic].back()->Fill(relPt1bVec[ic], relPt12Vec[ic], evtWeight*scaleMC);
              h2_relPt12_vs_relPt2bVec[ic].back()->Fill(relPt2bVec[ic], relPt12Vec[ic], evtWeight*scaleMC);
              h2_relPt1b_vs_relPt2bVec[ic].back()->Fill(relPt2bVec[ic], relPt1bVec[ic], evtWeight*scaleMC);
           }
           if( nTops == 2 ){
              h1_MT2_nTopsEQ2Vec.back()->Fill(MT2_nTopsEQ2, evtWeight*scaleMC);
              h1_MT2_TopAndbLept_nTopsEQ2Vec.back()->Fill(MT2_TopAndbLept_nTopsEQ2, evtWeight*scaleMC);
              h1_minDphiLeptMET_nTopsEQ2Vec.back()->Fill(minDphiLeptMETVec[1], evtWeight*scaleMC);

              h2_mTW1_vs_mTtop1_nTopsEQ2Vec.back()->Fill(mTtopVec[0], mTWVec[0], evtWeight*scaleMC);
              h2_mTb1_vs_mTtop1_nTopsEQ2Vec.back()->Fill(mTtopVec[0], mTbVec[0], evtWeight*scaleMC);
              h2_mTb1_vs_mTW1_nTopsEQ2Vec.back()->Fill(mTWVec[0], mTbVec[0], evtWeight*scaleMC);

              h2_mTW2_vs_mTtop2_nTopsEQ2Vec.back()->Fill(mTtopVec[1], mTWVec[1], evtWeight*scaleMC);
              h2_mTb2_vs_mTtop2_nTopsEQ2Vec.back()->Fill(mTtopVec[1], mTbVec[1], evtWeight*scaleMC);
              h2_mTb2_vs_mTW2_nTopsEQ2Vec.back()->Fill(mTWVec[1], mTbVec[1], evtWeight*scaleMC);

              h2_mTtop2_vs_mTtop1_nTopsEQ2Vec.back()->Fill(mTtopVec[0], mTtopVec[1], evtWeight*scaleMC);
              h2_mTW2_vs_mTtop1_nTopsEQ2Vec.back()->Fill(mTtopVec[0], mTWVec[1], evtWeight*scaleMC);
              h2_mTb2_vs_mTtop1_nTopsEQ2Vec.back()->Fill(mTtopVec[0], mTbVec[1], evtWeight*scaleMC);

              h2_mTtop2_vs_mTW1_nTopsEQ2Vec.back()->Fill(mTWVec[0], mTtopVec[1], evtWeight*scaleMC);
              h2_mTW2_vs_mTW1_nTopsEQ2Vec.back()->Fill(mTWVec[0], mTWVec[1], evtWeight*scaleMC);
              h2_mTb2_vs_mTW1_nTopsEQ2Vec.back()->Fill(mTWVec[0], mTbVec[1], evtWeight*scaleMC);

              h2_mTtop2_vs_mTb1_nTopsEQ2Vec.back()->Fill(mTbVec[0], mTtopVec[1], evtWeight*scaleMC);
              h2_mTW2_vs_mTb1_nTopsEQ2Vec.back()->Fill(mTbVec[0], mTWVec[1], evtWeight*scaleMC);
              h2_mTb2_vs_mTb1_nTopsEQ2Vec.back()->Fill(mTbVec[0], mTbVec[1], evtWeight*scaleMC);
           }
  
           bool isGenLepton = false; if( (cntlepton != 0 || cnttauHad != 0) ) isGenLepton = true; 
           if( nTops ==2 ){
              double minDR_triplets = -1, minDphi_triplets = -1;
              double minDR_Rsys = -1, minDphi_Rsys = -1;
              double minMTj = 999, minDphi_met = 999;
              int pickedIdx_triplets = -1, pickedIdx_Rsys = -1, pickedIdx_minMTj = -1, pickedIdx_minDphi_met = -1;
              for(unsigned int ij=0; ij<jetsLVec_forTagger.size(); ij++){
                 double deltaR = jetsLVec_forTagger.at(ij).DeltaR(pickedMaxPtLeptORprongLVec);
                 double deltaPhi = std::abs(jetsLVec_forTagger.at(ij).DeltaPhi(pickedMaxPtLeptORprongLVec));
                 double MTj = calcMT(jetsLVec_forTagger.at(ij), metLVec);
                 double deltaPhi_met = std::abs(jetsLVec_forTagger.at(ij).DeltaPhi(metLVec));
//                 if( std::find(type3Ptr->allIdxCacheVec.begin(), type3Ptr->allIdxCacheVec.end(), ij) != type3Ptr->allIdxCacheVec.end() ){
                 if( std::find(type3Ptr->bestTopJetComb.begin(), type3Ptr->bestTopJetComb.end(), ij) != type3Ptr->bestTopJetComb.end() ){
                    if( minDR_triplets == -1 || minDR_triplets > deltaR ){ minDR_triplets = deltaR; pickedIdx_triplets = ij; }
                    if( minDphi_triplets == -1 || minDphi_triplets > deltaPhi ) minDphi_triplets = deltaPhi;
                 }else{
                    bool inATopCand = false;
                    if( std::find(type3Ptr->allIdxCacheVec.begin(), type3Ptr->allIdxCacheVec.end(), ij) != type3Ptr->allIdxCacheVec.end() ){
                       if( minDR_Rsys == -1 || minDR_Rsys > deltaR ){ minDR_Rsys = deltaR; pickedIdx_Rsys = ij; }
                       if( minDphi_Rsys == -1 || minDphi_Rsys > deltaPhi ) minDphi_Rsys = deltaPhi;
                       inATopCand = true;
                    }
                    std::vector<TLorentzVector> dummyLVec; std::vector<double> dummyBtagVec;
                    double minDR_Rsys_triplets = -1;
                    for(unsigned int it=0; it<type3Ptr->bestTopJetComb.size(); it++){
                       double perDR_Rsys_triplets = jetsLVec_forTagger.at(ij).DeltaR( jetsLVec_forTagger.at(type3Ptr->bestTopJetComb[it]) );
//                    for(unsigned int it=0; it<type3Ptr->allIdxCacheVec.size(); it++){
//                       double perDR_Rsys_triplets = jetsLVec_forTagger.at(ij).DeltaR( jetsLVec_forTagger.at(type3Ptr->allIdxCacheVec[it]) );
                       if( minDR_Rsys_triplets == -1 || minDR_Rsys_triplets > perDR_Rsys_triplets ) minDR_Rsys_triplets = perDR_Rsys_triplets;
                    }
                    dummyLVec.push_back(jetsLVec_forTagger.at(ij)); dummyBtagVec.push_back(recoJetsBtag_forTagger.at(ij));
                    int cntCSVS = AnaFunctions::countCSVS(dummyLVec, dummyBtagVec, AnaConsts::cutCSVS, AnaConsts::bTagArr);
//                    if( !cntCSVS && minDR_Rsys_triplets > 1.5 ){
//                    if( minDR_Rsys_triplets > 1.5 ){
                    if( !inATopCand ){
                       if( minMTj > MTj ){ minMTj = MTj; pickedIdx_minMTj = ij; }
                       if( minDphi_met > deltaPhi_met ){ minDphi_met = deltaPhi_met; pickedIdx_minDphi_met = ij; }
                    }
                 }
              }

              h1_minMTj_lepJetVec.back()->Fill(minMTj, evtWeight*scaleMC);
              h1_minDphi_met_lepJetVec.back()->Fill(minDphi_met, evtWeight*scaleMC);

//              if( isGenLepton ){
              if( cntlepton + cnttauHad == 1 ){
                 if( pickedIdx_triplets != -1 ){
                    double relPt_triplets = pickedMaxPtLeptORprongLVec.Pt()/jetsLVec_forTagger.at(pickedIdx_triplets).Pt();
                    h1_relPt_genLeptORprongOVERjetPt_tripletsVec.back()->Fill(relPt_triplets, evtWeight*scaleMC);
                    h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPtVec.back()->Fill(pickedMaxPtLeptORprongLVec.Pt(), relPt_triplets, evtWeight*scaleMC);
                 }
                 if( pickedIdx_Rsys != -1 ){
                    double relPt_Rsys = pickedMaxPtLeptORprongLVec.Pt()/jetsLVec_forTagger.at(pickedIdx_Rsys).Pt();
                    h1_relPt_genLeptORprongOVERjetPt_RsysVec.back()->Fill(relPt_Rsys, evtWeight*scaleMC);
                    h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPtVec.back()->Fill(pickedMaxPtLeptORprongLVec.Pt(), relPt_Rsys, evtWeight*scaleMC);
                 }

                 double deltaPhi_met = std::abs(pickedMaxPtLeptORprongLVec.DeltaPhi(metLVec));
      
                 h1_genLeptORprong_PtVec.back()->Fill(pickedMaxPtLeptORprongLVec.Pt(), evtWeight*scaleMC);
                 h1_genLeptORprong_EtaVec.back()->Fill(pickedMaxPtLeptORprongLVec.Eta(), evtWeight*scaleMC);
                 h1_genLeptORprong_PhiVec.back()->Fill(pickedMaxPtLeptORprongLVec.Phi(), evtWeight*scaleMC);
         
                 h1_genLeptORprong_minDR_tripletsVec.back()->Fill(minDR_triplets, evtWeight*scaleMC);
                 h1_genLeptORprong_minDphi_tripletsVec.back()->Fill(minDphi_triplets, evtWeight*scaleMC);
         
                 h1_genLeptORprong_minDR_RsysVec.back()->Fill(minDR_Rsys, evtWeight*scaleMC);
                 h1_genLeptORprong_minDphi_RsysVec.back()->Fill(minDphi_Rsys, evtWeight*scaleMC);
         
                 h1_genLeptORprong_dPhi_metVec.back()->Fill(deltaPhi_met, evtWeight*scaleMC);
      
                 h2_genLeptORprong_minDR_Rsys_vs_tripletsVec.back()->Fill(minDR_triplets, minDR_Rsys, evtWeight*scaleMC);
                 h2_genLeptORprong_minDphi_Rsys_vs_tripletsVec.back()->Fill(minDphi_triplets, minDphi_Rsys, evtWeight*scaleMC);
    
                 h2_minMTj_pickedIdx_vs_minDR_pickedIdxVec.back()->Fill(pickedIdx_Rsys, pickedIdx_minMTj, evtWeight*scaleMC);
   
                 h2_minDphi_met_pickedIdx_vs_minDR_pickedIdxVec.back()->Fill(pickedIdx_Rsys, pickedIdx_minDphi_met, evtWeight*scaleMC);
              }
           }
// Search region in bins of
// nbJets: 1, 2, >=3
// nTops : 1, 2, >=3
           int nbJetsCopy = nbJets, nTopsCopy = nTops;
           if( nbJetsCopy >=3 ) nbJetsCopy = 3; if( nTopsCopy >=3 ) nTopsCopy = 3;

           h2_evtCnt_nbJets_vs_nTopsVec.back()->Fill(nTopsCopy, nbJetsCopy, evtWeight*scaleMC);
           if( !keyStringT.Contains("Signal") ){
              h2_evtCnt_sumSM_nbJets_vs_nTops->Fill(nTopsCopy, nbJetsCopy, evtWeight*scaleMC);
           }
        }
     }
     pcaVec.back()->MakePrincipals();
     pcaVec.back()->MakeHistograms();
//     pcaVec.back()->MakeCode("pca_"+keyWordVec.back());
  }

  for(int iSR = 0 ; iSR < nSR; iSR++){
     for(int iSR_met =0; iSR_met < nSR_met; iSR_met++){
        cnt_passSRmet_WeightedErrorScaledMCVec.back()[iSR][iSR_met] = sqrt(cnt_passSRmet_WeightedErrorScaledMCVec.back()[iSR][iSR_met]);
     }
  }

  cnt_passLeptVeto_WeightedErrorScaledMC = sqrt(cnt_passLeptVeto_WeightedErrorScaledMC);
  cnt_passnJets_WeightedErrorScaledMC = sqrt(cnt_passnJets_WeightedErrorScaledMC);
  cnt_passdPhis_WeightedErrorScaledMC = sqrt(cnt_passdPhis_WeightedErrorScaledMC);
  cnt_passBJets_WeightedErrorScaledMC = sqrt(cnt_passBJets_WeightedErrorScaledMC);
  cnt_passMET_WeightedErrorScaledMC = sqrt(cnt_passMET_WeightedErrorScaledMC);
  cnt_passTagger_WeightedErrorScaledMC = sqrt(cnt_passTagger_WeightedErrorScaledMC);
  cnt_passBaseline_WeightedErrorScaledMC = sqrt(cnt_passBaseline_WeightedErrorScaledMC);
  cnt_passMore_WeightedErrorScaledMC = sqrt(cnt_passMore_WeightedErrorScaledMC);

  cnt_onelepton_WeightedErrorScaledMC = sqrt(cnt_onelepton_WeightedErrorScaledMC);
  cnt_Wtauhad_WeightedErrorScaledMC = sqrt(cnt_Wtauhad_WeightedErrorScaledMC);
  cnt_dileptons_WeightedErrorScaledMC = sqrt(cnt_dileptons_WeightedErrorScaledMC);
  cnt_diWtauhads_WeightedErrorScaledMC = sqrt(cnt_diWtauhads_WeightedErrorScaledMC);
  cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC = sqrt(cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC);
  cnt_allHad_WeightedErrorScaledMC = sqrt(cnt_allHad_WeightedErrorScaledMC);
                               
  std::cout<<"\n\n"<<sampleKeyString.c_str()<<"_cnt_passLeptVeto : "<<cnt_passLeptVeto_WeightedScaledMC<<" +- "<<cnt_passLeptVeto_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passnJets : "<<cnt_passnJets_WeightedScaledMC<<" +- "<<cnt_passnJets_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passdPhis : "<<cnt_passdPhis_WeightedScaledMC<<" +- "<<cnt_passdPhis_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passBJets : "<<cnt_passBJets_WeightedScaledMC<<" +- "<<cnt_passBJets_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passMET : "<<cnt_passMET_WeightedScaledMC<<" +- "<<cnt_passMET_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passTagger : "<<cnt_passTagger_WeightedScaledMC<<" +- "<<cnt_passTagger_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passBaseline : "<<cnt_passBaseline_WeightedScaledMC<<" +- "<<cnt_passBaseline_WeightedErrorScaledMC<<std::endl;
  std::cout<<sampleKeyString.c_str()<<"_cnt_passMore : "<<cnt_passMore_WeightedScaledMC<<" +- "<<cnt_passMore_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_onelepton : "<<cnt_onelepton_WeightedScaledMC<<" +- "<<cnt_onelepton_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_Wtauhad : "<<cnt_Wtauhad_WeightedScaledMC<<" +- "<<cnt_Wtauhad_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_dileptons : "<<cnt_dileptons_WeightedScaledMC<<" +- "<<cnt_dileptons_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_diWtauhads : "<<cnt_diWtauhads_WeightedScaledMC<<" +- "<<cnt_diWtauhads_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_dileptons_inc_Wtauhad : "<<cnt_dileptons_inc_Wtauhad_WeightedScaledMC<<" +- "<<cnt_dileptons_inc_Wtauhad_WeightedErrorScaledMC<<std::endl;
  std::cout<<"  --> "<<sampleKeyString.c_str()<<"_cnt_allHad : "<<cnt_allHad_WeightedScaledMC<<" +- "<<cnt_allHad_WeightedErrorScaledMC<<std::endl;
/*  
  std::cout<<"==> In search bins ... "<<std::endl;
  for(int iSR = 0; iSR < nSR; iSR ++ ){
     std::cout<<"  --> nbJets"<<keyStr_nbJets_SR[iSR].c_str()<<"_nTops"<<keyStr_nTops_SR[iSR].c_str()<<":"<<std::endl;
     for(int iSR_met = 0; iSR_met < nSR_met; iSR_met ++){
        std::cout<<"      "<<sampleKeyString.c_str()<<"_MET_"<<keyStr_met_SR[iSR_met]<<" : "<<cnt_passSRmet_WeightedScaledMCVec.back()[iSR][iSR_met]<<" +- "<<cnt_passSRmet_WeightedErrorScaledMCVec.back()[iSR][iSR_met]<<std::endl;
     }
  }
*/
}

void explore(){

   initMCscales();

   AnaFunctions::prepareTopTagger();

   NTupleReader *tr = 0;

   std::vector<TTree*> treeVec;
   std::vector<std::string> subSampleKeysVec;
 
   char filenames[500], names[500];
   std::vector<std::string> filesDataVec;
   std::vector<std::string> filesTTJetsVec,
                            filesWJetsToLNu_HT_600toInfVec, filesWJetsToLNu_HT_400to600Vec, filesWJetsToLNu_HT_200to400Vec, filesWJetsToLNu_HT_100to200Vec,
                            filesZJetsToNuNu_HT_600toInfVec, filesZJetsToNuNu_HT_400to600Vec, filesZJetsToNuNu_HT_200to400Vec, filesZJetsToNuNu_HT_100to200Vec,
                            filesDYJetsToLL_HT_600toInfVec, filesDYJetsToLL_HT_400to600Vec, filesDYJetsToLL_HT_200to400Vec, filesDYJetsToLL_HT_100to200Vec,
                            filesQCD_HT_1000toInfVec, filesQCD_HT_500to1000Vec, filesQCD_HT_250to500Vec,
                            filesT_tWVec, filesTbar_tWVec, filesTTZVec;
   std::vector<std::string> filesSignal_T1tttt_mGluino1200_mLSP800Vec,
                            filesSignal_T1tttt_mGluino1500_mLSP100Vec,
                            filesSignal_T5tttt_mGluino1300_mStop300_mChi280Vec, filesSignal_T5tttt_mGluino1300_mStop300_mCh285Vec,
                            filesSignal_T1bbbb_mGluino1000_mLSP900Vec,
                            filesSignal_T1bbbb_mGluino1500_mLSP100Vec,
                            filesSignal_T2tt_mStop425_mLSP325Vec,
                            filesSignal_T2tt_mStop500_mLSP325Vec,
                            filesSignal_T2tt_mStop650_mLSP325Vec,
                            filesSignal_T2tt_mStop850_mLSP100Vec,
                            filesSignal_T2bb_mSbottom600_mLSP580Vec,
                            filesSignal_T2bb_mSbottom900_mLSP100Vec;
 
   ifstream finData("rootlist_Data_METPD.txt"); while( finData.getline(filenames, 500) ){ filesDataVec.push_back(filenames); }
 
   ifstream finTTJets("rootlist_TTJets.txt"); while( finTTJets.getline(filenames, 500) ){ filesTTJetsVec.push_back(filenames); }
 
   ifstream finWJetsToLNu_HT_600toInf("rootlist_WJetsToLNu_HT_600toInf.txt"); while( finWJetsToLNu_HT_600toInf.getline(filenames, 500) ){ filesWJetsToLNu_HT_600toInfVec.push_back(filenames); }
   ifstream finWJetsToLNu_HT_400to600("rootlist_WJetsToLNu_HT_400to600.txt"); while( finWJetsToLNu_HT_400to600.getline(filenames, 500) ){ filesWJetsToLNu_HT_400to600Vec.push_back(filenames); }
   ifstream finWJetsToLNu_HT_200to400("rootlist_WJetsToLNu_HT_200to400.txt"); while( finWJetsToLNu_HT_200to400.getline(filenames, 500) ){ filesWJetsToLNu_HT_200to400Vec.push_back(filenames); }
   ifstream finWJetsToLNu_HT_100to200("rootlist_WJetsToLNu_HT_100to200.txt"); while( finWJetsToLNu_HT_100to200.getline(filenames, 500) ){ filesWJetsToLNu_HT_100to200Vec.push_back(filenames); }
 
   ifstream finZJetsToNuNu_HT_600toInf("rootlist_ZJetsToNuNu_HT_600toInf.txt"); while( finZJetsToNuNu_HT_600toInf.getline(filenames, 500) ){ filesZJetsToNuNu_HT_600toInfVec.push_back(filenames); }
   ifstream finZJetsToNuNu_HT_400to600("rootlist_ZJetsToNuNu_HT_400to600.txt"); while( finZJetsToNuNu_HT_400to600.getline(filenames, 500) ){ filesZJetsToNuNu_HT_400to600Vec.push_back(filenames); }
   ifstream finZJetsToNuNu_HT_200to400("rootlist_ZJetsToNuNu_HT_200to400.txt"); while( finZJetsToNuNu_HT_200to400.getline(filenames, 500) ){ filesZJetsToNuNu_HT_200to400Vec.push_back(filenames); }
   ifstream finZJetsToNuNu_HT_100to200("rootlist_ZJetsToNuNu_HT_100to200.txt"); while( finZJetsToNuNu_HT_100to200.getline(filenames, 500) ){ filesZJetsToNuNu_HT_100to200Vec.push_back(filenames); }
 
   ifstream finDYJetsToLL_HT_600toInf("rootlist_DYJetsToLL_HT_600toInf.txt"); while( finDYJetsToLL_HT_600toInf.getline(filenames, 500) ){ filesDYJetsToLL_HT_600toInfVec.push_back(filenames); }
   ifstream finDYJetsToLL_HT_400to600("rootlist_DYJetsToLL_HT_400to600.txt"); while( finDYJetsToLL_HT_400to600.getline(filenames, 500) ){ filesDYJetsToLL_HT_400to600Vec.push_back(filenames); }
   ifstream finDYJetsToLL_HT_200to400("rootlist_DYJetsToLL_HT_200to400.txt"); while( finDYJetsToLL_HT_200to400.getline(filenames, 500) ){ filesDYJetsToLL_HT_200to400Vec.push_back(filenames); }
   ifstream finDYJetsToLL_HT_100to200("rootlist_DYJetsToLL_HT_100to200.txt"); while( finDYJetsToLL_HT_100to200.getline(filenames, 500) ){ filesDYJetsToLL_HT_100to200Vec.push_back(filenames); }
 
   ifstream finQCD_HT_1000toInf("rootlist_QCD_HT_1000toInf.txt"); while( finQCD_HT_1000toInf.getline(filenames, 500) ){ filesQCD_HT_1000toInfVec.push_back(filenames); }
   ifstream finQCD_HT_500to1000("rootlist_QCD_HT_500to1000.txt"); while( finQCD_HT_500to1000.getline(filenames, 500) ){ filesQCD_HT_500to1000Vec.push_back(filenames); }
   ifstream finQCD_HT_250to500("rootlist_QCD_HT_250to500.txt"); while( finQCD_HT_250to500.getline(filenames, 500) ){ filesQCD_HT_250to500Vec.push_back(filenames); }
 
   ifstream finT_tW("rootlist_T_tW.txt"); while( finT_tW.getline(filenames, 500) ){ filesT_tWVec.push_back(filenames); }
   ifstream finTbar_tW("rootlist_Tbar_tW.txt"); while( finTbar_tW.getline(filenames, 500) ){ filesTbar_tWVec.push_back(filenames); }
   ifstream finTTZ("rootlist_TTZ.txt"); while( finTTZ.getline(filenames, 500) ){ filesTTZVec.push_back(filenames); }
 
   ifstream finSignal_T1tttt_mGluino1200_mLSP800("rootlist_T1tttt_mGluino1200_mLSP800.txt"); while(finSignal_T1tttt_mGluino1200_mLSP800.getline(filenames, 500) ){ filesSignal_T1tttt_mGluino1200_mLSP800Vec.push_back(filenames); }
   ifstream finSignal_T1tttt_mGluino1500_mLSP100("rootlist_T1tttt_mGluino1500_mLSP100.txt"); while(finSignal_T1tttt_mGluino1500_mLSP100.getline(filenames, 500) ){ filesSignal_T1tttt_mGluino1500_mLSP100Vec.push_back(filenames); }
   ifstream finSignal_T5tttt_mGluino1300_mStop300_mChi280("rootlist_T5tttt_mGluino1300_mStop300_mChi280.txt"); while(finSignal_T5tttt_mGluino1300_mStop300_mChi280.getline(filenames, 500) ){filesSignal_T5tttt_mGluino1300_mStop300_mChi280Vec.push_back(filenames); }
   ifstream finSignal_T5tttt_mGluino1300_mStop300_mCh285("rootlist_T5tttt_mGluino1300_mStop300_mCh285.txt"); while(finSignal_T5tttt_mGluino1300_mStop300_mCh285.getline(filenames, 500) ){ filesSignal_T5tttt_mGluino1300_mStop300_mCh285Vec.push_back(filenames); }
   ifstream finSignal_T1bbbb_mGluino1000_mLSP900("rootlist_T1bbbb_mGluino1000_mLSP900.txt"); while(finSignal_T1bbbb_mGluino1000_mLSP900.getline(filenames, 500) ){ filesSignal_T1bbbb_mGluino1000_mLSP900Vec.push_back(filenames); }
   ifstream finSignal_T1bbbb_mGluino1500_mLSP100("rootlist_T1bbbb_mGluino1500_mLSP100.txt"); while(finSignal_T1bbbb_mGluino1500_mLSP100.getline(filenames, 500) ){ filesSignal_T1bbbb_mGluino1500_mLSP100Vec.push_back(filenames); }
   ifstream finSignal_T2tt_mStop425_mLSP325("rootlist_T2tt_mStop425_mLSP325.txt"); while(finSignal_T2tt_mStop425_mLSP325.getline(filenames, 500) ){ filesSignal_T2tt_mStop425_mLSP325Vec.push_back(filenames); }
   ifstream finSignal_T2tt_mStop500_mLSP325("rootlist_T2tt_mStop500_mLSP325.txt"); while(finSignal_T2tt_mStop500_mLSP325.getline(filenames, 500) ){ filesSignal_T2tt_mStop500_mLSP325Vec.push_back(filenames); }
   ifstream finSignal_T2tt_mStop650_mLSP325("rootlist_T2tt_mStop650_mLSP325.txt"); while(finSignal_T2tt_mStop650_mLSP325.getline(filenames, 500) ){ filesSignal_T2tt_mStop650_mLSP325Vec.push_back(filenames); }
   ifstream finSignal_T2tt_mStop850_mLSP100("rootlist_T2tt_mStop850_mLSP100.txt"); while(finSignal_T2tt_mStop850_mLSP100.getline(filenames, 500) ){ filesSignal_T2tt_mStop850_mLSP100Vec.push_back(filenames); }
   ifstream finSignal_T2bb_mSbottom600_mLSP580("rootlist_T2bb_mSbottom600_mLSP580.txt"); while(finSignal_T2bb_mSbottom600_mLSP580.getline(filenames, 500) ){ filesSignal_T2bb_mSbottom600_mLSP580Vec.push_back(filenames); }
   ifstream finSignal_T2bb_mSbottom900_mLSP100("rootlist_T2bb_mSbottom900_mLSP100.txt"); while(finSignal_T2bb_mSbottom900_mLSP100.getline(filenames, 500) ){ filesSignal_T2bb_mSbottom900_mLSP100Vec.push_back(filenames); }

   std::cout<<"\n"<<std::endl; timer.Print(); timer.Start();

   h2_evtCnt_sumSM_nbJets_vs_nTops = new TH2D("h2_evtCnt_sumSM_nbJets_vs_nTops", "SumSM: event counts nbJets versus nTops; nTops; nbJets", 4, 0, 4, 3, 1, 4);

   cnt_passSRmet_sumSM_WeightedScaledMCVec.resize(nSR); cnt_passSRmet_sumSM_WeightedErrorScaledMCVec.resize(nSR);
   for(int iSR=0; iSR < nSR; iSR++){
      cnt_passSRmet_sumSM_WeightedScaledMCVec[iSR].resize(nSR_met);  cnt_passSRmet_sumSM_WeightedErrorScaledMCVec[iSR].resize(nSR_met);
   }

#include "bkgSamples.h"

   std::cout<<"==> In search bins ... "<<std::endl;
   for(int iSR = 0; iSR < nSR; iSR ++ ){
      std::cout<<"  --> nbJets"<<keyStr_nbJets_SR[iSR].c_str()<<"_nTops"<<keyStr_nTops_SR[iSR].c_str()<<":"<<std::endl;
      for(int iSR_met = 0; iSR_met < nSR_met; iSR_met ++){
         cnt_passSRmet_sumSM_WeightedErrorScaledMCVec[iSR][iSR_met] = sqrt(cnt_passSRmet_sumSM_WeightedErrorScaledMCVec[iSR][iSR_met]);
         std::cout<<"    MET_"<<keyStr_met_SR[iSR_met]<<" : "<<cnt_passSRmet_sumSM_WeightedScaledMCVec[iSR][iSR_met]<<" +- "<<cnt_passSRmet_sumSM_WeightedErrorScaledMCVec[iSR][iSR_met]<<std::endl;
      }
   }

#include "sigSamples.h"

// Plotting
   setTDRStyle();

  //  For the Global title:
   tdrStyle->SetOptTitle(1);
   tdrStyle->SetTitleFont(42,"");
   tdrStyle->SetTitleColor(1);
   tdrStyle->SetTitleTextColor(1);
   tdrStyle->SetTitleFillColor(0);
   tdrStyle->SetTitleFontSize(0.1);
   tdrStyle->SetTitleAlign(13);
   tdrStyle->SetTitleX(0.6);
   tdrStyle->SetTitleH(0.05);
   tdrStyle->SetTitleBorderSize(0);
   tdrStyle->SetTitleAlign(13);
   tdrStyle->SetTitleX(0.19);
   tdrStyle->SetTitleH(0.038);
                                                                         
     //  For the frame
   tdrStyle->SetFrameBorderMode(0);
   tdrStyle->SetFrameBorderSize(1);
   tdrStyle->SetFrameFillColor(kBlack);
   tdrStyle->SetFrameFillStyle(0);
   tdrStyle->SetFrameLineColor(kBlack);
   tdrStyle->SetFrameLineStyle(0);
   tdrStyle->SetFrameLineWidth(1);
 
      //  For the Pad
   tdrStyle->SetPadBorderMode(0);
   tdrStyle->SetPadColor(kWhite);
   tdrStyle->SetPadGridX(false);
   tdrStyle->SetPadGridY(false);
   tdrStyle->SetGridColor(0);
   tdrStyle->SetGridStyle(3);
   tdrStyle->SetGridWidth(1);
 
      //  For the axis
   tdrStyle->SetAxisColor(1,"XYZ");
   tdrStyle->SetTickLength(0.03,"XYZ");
   tdrStyle->SetNdivisions(505,"XYZ");
   tdrStyle->SetPadTickX(1);
   tdrStyle->SetPadTickY(1);
   tdrStyle->SetStripDecimals(kFALSE);
 
   tdrStyle->SetLabelSize(0.050, "XYZ");
 
   tdrStyle->SetPadTopMargin(0.1); tdrStyle->SetPadBottomMargin(0.15);
   tdrStyle->SetPadLeftMargin(0.15); tdrStyle->SetPadRightMargin(0.15);
 
   tdrStyle->SetOptStat(1111);
 
   tdrStyle->SetHistLineWidth(1);
 
   tdrStyle->SetPaintTextFormat("4.2f");

   tdrStyle->SetTitleXOffset(5.50); tdrStyle->SetTitleYOffset(6.50);

   TCanvas *cs = new TCanvas("cs", "cs", 1200, 900);
   int divW=4, divH=3;
   cs->Divide(divW, divH);

   cs->Print("allINone_"+treeStrTtype+".pdf[");
   for(unsigned int ic=0; ic<h2_evtCnt_nbJets_vs_nTopsVec.size(); ic++){ h2_evtCnt_nbJets_vs_nTopsVec[ic]->SetMarkerSize(h2_evtCnt_nbJets_vs_nTopsVec[ic]->GetMarkerSize()*2.0); }
   draw2DallINone(cs, divW*divH, h2_evtCnt_nbJets_vs_nTopsVec, "colz text"); cs->Print("allINone_"+treeStrTtype+".pdf");
   for(unsigned int ic=0; ic<h1_cutFlowVec.size(); ic++){ h1_cutFlowVec[ic]->LabelsDeflate(); h1_cutFlowVec[ic]->LabelsOption("v"); h1_cutFlowVec[ic]->SetMarkerSize(h1_cutFlowVec[ic]->GetMarkerSize()*1.5); }
   for(unsigned int ic=0; ic<h1_cutFlow_auxVec.size(); ic++){ h1_cutFlow_auxVec[ic]->LabelsDeflate(); h1_cutFlow_auxVec[ic]->LabelsOption("v"); h1_cutFlow_auxVec[ic]->SetMarkerSize(h1_cutFlow_auxVec[ic]->GetMarkerSize()*1.5); }
   draw1DallINone(cs, divW*divH, h1_cutFlowVec, 1, "text"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_cutFlow_auxVec, 1, "text"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_dalitzVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_gen_m23overm123vsarctanm13overm12Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_m23overm123vsarctanm13overm12Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_pt_gentbVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_eta_gentbVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_pt_genrbVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_eta_genrbVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_csvs_fakebVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_pt_genrb_match0Vec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_eta_genrb_match0Vec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_genTopsVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_genb_genTopVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_genbsVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_top_metVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw1DallINone(cs, divW*divH, h1_minDphi_tJet_metVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw1DallINone(cs, divW*divH, h1_dR_gen_reco_topVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dR_b_topVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_b_topVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mt_b_vs_dPhi_b_topVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_lept_metVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_nu_metVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_sum_lept_nu_metVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_rJet_metVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw1DallINone(cs, divW*divH, h1_dPhi_rJet_bVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw2DallINone(cs, divW*divH, h2_dPhi_rJet_b_vs_metVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_MT_tJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw1DallINone(cs, divW*divH, h1_MT_tDiJetsVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf"); 

   draw1DallINone(cs, divW*divH, h1_pt_genWdau1_nrJets_EQ1Vec, 4, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_eta_genWdau1_nrJets_EQ1Vec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_pt_genWdau2_nrJets_EQ1Vec, 4, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_eta_genWdau2_nrJets_EQ1Vec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_mass_rbJet_mtch0_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_mass_rbJet_mtch1_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dR_rbJet_mtch0_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dR_rbJet_mtch1_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_rJet_met_mtch0_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_dPhi_rJet_met_mtch1_nrJets_EQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_dPhi_sum_METVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT_sumVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_nComb_top_brJetsVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_nComb_lept_brJetVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_mass_rTopVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_ori_MT_lept_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT_lept_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT_bJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_MT_bJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT2_lept_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_nComb_had_brJetVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_ori_MT2_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT2_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_nComb_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_top_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_vs_MT_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_vs_MT2_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_vs_MT2_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_had_vs_top_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_aft_PCA_MT_had_vs_top_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_ori_MT2_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT2_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT2_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_mTcomb_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_aft_PCA_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_bJet_vs_mTcomb_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_lept_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT2_lept_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT_lept_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTcomb_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_mTcomb_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_MT2_lept_vs_had_aft_MT_cuts_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT2_same_lept_had_aft_MT_cuts_brJetsVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");


   draw1DallINone(cs, divW*divH, h1_adj_MT_lept_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT2_lept_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT2_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_adj_MT_lept_vs_top_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_vs_MT_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_vs_MT2_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_vs_MT2_lept_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_had_vs_top_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT2_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT2_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_lept_vs_mTcomb_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_lept_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT2_lept_vs_MT_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT_lept_vs_MT2_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_mTcomb_lept_vs_had_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_mTcomb_same_lept_had_brJetVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_adj_MT2_lept_vs_had_aft_MT_cuts_brJetVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_adj_MT2_same_lept_had_aft_MT_cuts_brJetsVec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");


   draw1DallINone(cs, divW*divH, h1_gen1b_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen1b_deltaRVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen1MET_deltaPhiVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen2b_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen2b_deltaRVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen2MET_deltaPhiVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_gen12_deltaRVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_Wjets_bTagCatsVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_Topjets_bTagCatsVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_CSVVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_bTagged_Wjets_minDR_genb_vs_minDR_genWVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_bTagged_Wjets_mass_vs_minDR_genbVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_nobTagged_Wjets_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_bTagged_rndmComb_Wjets_ptVec, 1, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_bTagged_sameTop_Wjets_ptVec, 1, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_minDR_otherbJetsVec, 1, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf"); 
   draw2DallINone(cs, divW*divH, h2_nobTagged_Wjets_minDR_genb_vs_minDR_genWVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_bTagged_Topjets_minDR_genb_vs_minDR_genTopVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_nobTagged_Topjets_minDR_genb_vs_minDR_genTopVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_bTagged_Topjets_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_nobTagged_Topjets_massVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_deltaR1b_genDausVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_deltaR2b_genDausVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_bTagged_Wjets_deltaR12_genDausVec, 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw1DallINone(cs, divW*divH, h1_MT2_nTopsEQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_minDphiLeptMET_nTopsEQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_MT2_TopAndbLept_nTopsEQ2Vec, 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_mTW1_vs_mTtop1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb1_vs_mTtop1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb1_vs_mTW1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_mTW2_vs_mTtop2_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb2_vs_mTtop2_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb2_vs_mTW2_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_mTtop2_vs_mTtop1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTW2_vs_mTtop1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb2_vs_mTtop1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_mTtop2_vs_mTW1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTW2_vs_mTW1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb2_vs_mTW1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   draw2DallINone(cs, divW*divH, h2_mTtop2_vs_mTb1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTW2_vs_mTb1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_mTb2_vs_mTb1_nTopsEQ2Vec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

   for(unsigned int ic=0; ic< h1_topCand_MVec.size(); ic++){
      draw1DallINone(cs, divW*divH, h1_topCand_MVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_topCand_PtVec[ic], 2, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_topCand_EtaVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_topCand_nJetsVec[ic], 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_W_MVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_W_PtVec[ic], 2, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_W_EtaVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_W_dausDRVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_W_dausDR_versus_W_MVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_b_MVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_b_PtVec[ic], 2, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_b_EtaVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_ratio_mW_over_mTopVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_mW_versus_mTopVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_ratio_mW_over_mTop_versus_mTopVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_mTtopVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_mTWVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_mTbVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_deltaR12Vec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_deltaR1bVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_deltaR2bVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw2DallINone(cs, divW*divH, h2_deltaR12_vs_deltaR1bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_deltaR12_vs_deltaR2bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_deltaR1b_vs_deltaR2bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw1DallINone(cs, divW*divH, h1_relPt12Vec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_relPt1bVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_relPt2bVec[ic], 2, ""); cs->Print("allINone_"+treeStrTtype+".pdf");

      draw2DallINone(cs, divW*divH, h2_relPt12_vs_relPt1bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_relPt12_vs_relPt2bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_relPt1b_vs_relPt2bVec[ic], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   }
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_PtVec, 4, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_EtaVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_PhiVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_minDR_tripletsVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_minDphi_tripletsVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_relPt_genLeptORprongOVERjetPt_tripletsVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_relPt_genLeptORprongOVERjetPt_triplets_vs_genLeptORprongPtVec); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_minDR_RsysVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_minDphi_RsysVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_relPt_genLeptORprongOVERjetPt_RsysVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_relPt_genLeptORprongOVERjetPt_Rsys_vs_genLeptORprongPtVec); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_genLeptORprong_dPhi_metVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_genLeptORprong_minDR_Rsys_vs_tripletsVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_genLeptORprong_minDphi_Rsys_vs_tripletsVec, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_minMTj_lepJetVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_minMTj_pickedIdx_vs_minDR_pickedIdxVec); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw1DallINone(cs, divW*divH, h1_minDphi_met_lepJetVec, 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   draw2DallINone(cs, divW*divH, h2_minDphi_met_pickedIdx_vs_minDR_pickedIdxVec); cs->Print("allINone_"+treeStrTtype+".pdf");
   for(int iSR=0; iSR<nSR; iSR++){
      draw1DallINone(cs, divW*divH, h1_nJetsVec[iSR], 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_nJetsRsysVec[iSR], 1, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_metVec[iSR], 4, "SetLogy"); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_MT2Vec[iSR], 4, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_MT2_vs_metVec[iSR], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
//      for(unsigned int is=0; is<h2_coarse_bin_MT2_vs_metVec[iSR].size(); is++){ h2_coarse_bin_MT2_vs_metVec[iSR][is]->SetMarkerSize(h2_coarse_bin_MT2_vs_metVec[iSR][is]->GetMarkerSize()*2.0); }
      draw2DallINone(cs, divW*divH, h2_coarse_bin_MT2_vs_metVec[iSR], "colz text"); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_mTcombVec[iSR], 3, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw1DallINone(cs, divW*divH, h1_HTVec[iSR], 3, ""); cs->Print("allINone_"+treeStrTtype+".pdf");
      draw2DallINone(cs, divW*divH, h2_met_vs_nJetsVec[iSR], ""); cs->Print("allINone_"+treeStrTtype+".pdf");
   }
   cs->Print("allINone_"+treeStrTtype+".pdf]");

   TCanvas *ct = new TCanvas("ct", "ct", 800, 800);
   ct->cd();
   h2_evtCnt_sumSM_nbJets_vs_nTops->SetMarkerSize(h2_evtCnt_sumSM_nbJets_vs_nTops->GetMarkerSize()*2.0);
   h2_evtCnt_sumSM_nbJets_vs_nTops->Draw("colz text");
   ct->Print("evtCnt_sumSM_nbJets_vs_nTops.pdf");

   TFile *exploreFile = new TFile("explore.root", "RECREATE");
   TDirectory * srDir = (TDirectory*) exploreFile->mkdir("SR");
   srDir->cd();
   for(unsigned int is = 0; is< h2_evtCnt_nbJets_vs_nTopsVec.size(); is++){
      h2_evtCnt_nbJets_vs_nTopsVec[is]->Write();
      for(int iSR = 0; iSR < nSR; iSR++){
         h1_nJetsVec[iSR][is]->Write();
         h1_metVec[iSR][is]->Write();
         h1_MT2Vec[iSR][is]->Write();
         h1_mTcombVec[iSR][is]->Write();
         h1_HTVec[iSR][is]->Write();
      }
      h2_MT_had_vs_top_brJetVec[is]->Write();
   }
   exploreFile->Write(); exploreFile->Close();

}

bool find_mother(int momIdx, int dauIdx, const std::vector<int> &genDecayIdxVec, const std::vector<int> &genDecayMomIdxVec){
   if( momIdx == -1 || dauIdx == -1 ) return false;

   if( dauIdx == momIdx ) return true;

   int thisIdx = dauIdx;
   while( thisIdx >=0 ){
      int momGenIdx = genDecayMomIdxVec[thisIdx];
      thisIdx = find_idx(momGenIdx, genDecayIdxVec);
      if( thisIdx == momIdx ) return true;
   }

   return false;
}

int find_idx(int genIdx, const std::vector<int> &genDecayIdxVec){
   for(int ig=0; ig<(int)genDecayIdxVec.size(); ig++){
      if( genDecayIdxVec[ig] == genIdx ) return ig;
   }
   return -1;
}

double calcMT(const TLorentzVector &objLVec, const TLorentzVector &metLVec){

   const double objMass = objLVec.M(), objPt = objLVec.Pt(), objPx = objLVec.Px(), objPy = objLVec.Py();
   const double met = metLVec.Pt(), metphi = metLVec.Phi();

   double mt = sqrt( objMass*objMass + 2*( met*sqrt(objMass*objMass + objPt*objPt) -( met*cos(metphi)*objPx + met*sin(metphi)*objPy ) ) );

   return mt;
}

void draw2DallINone(TCanvas *cs, const int lastPadIdx, const std::vector<TH2D*> &h2_inputVec, const TString optDrawStr){

  int cntPadIdx = 0;
  unsigned int ntype = keyStringCachedVec.size();

  for(unsigned int it=0; it<ntype; it++){
     if( it == 0 ) cntPadIdx = 0;
     TString keyStringT(keyStringCachedVec[it]);
     if( keyStringT.Contains("DYJets") || keyStringCachedVec[it] == "T_t" || keyStringT.Contains("Data") ) continue;
     cntPadIdx ++;
     if( cntPadIdx >= lastPadIdx ){ /*std::cout<<"Overwritten the last pad with index : "<<lastPadIdx<<"! Please make sure enough pads are created!"<<std::endl;*/ cntPadIdx = lastPadIdx; }

     cs->cd(cntPadIdx); TPad * pad = (TPad*) cs->GetPad(cntPadIdx); pad->Clear(); pad->SetLogy(0);
     if( optDrawStr.Length() == 0 ){
        h2_inputVec[it]->Draw("colz");
     }else{
        h2_inputVec[it]->Draw(optDrawStr);
     }
  }

  for(int ic=cntPadIdx+1; ic<=lastPadIdx; ic++){ cs->cd(ic); TPad * pad = (TPad*) cs->GetPad(ic); pad->Clear(); }
}

void draw1DallINone(TCanvas *cs, const int lastPadIdx, const std::vector<TH1D*> &h1_inputVec, const int nRebin, const TString optDrawStr){

  int statusLogy = 0;
  if( optDrawStr.Contains("SetLogy") ) statusLogy =1;
  int doNormalization = 0;
  if( optDrawStr.Contains("normalization") ) doNormalization =1;
  int drawText = 0;
  if( optDrawStr.Contains("text") ) drawText =1;

  int cntPadIdx = 0;
  std::vector<TH1D*> h1_stackedVec; TH1D * h1_dataCached =0, * h1_signal1Cached =0, *h1_signal2Cached =0;

  std::vector<TString> stackedKeyStringTVec;

  unsigned int ntype = keyStringCachedVec.size();

  for(unsigned int it=0; it<ntype; it++){

     TString keyStringT(keyStringCachedVec[it]);

     if( it == 0 ){ cntPadIdx = 0; h1_stackedVec.clear(); }

//     if( keyStringT.Contains("DYJets") || keyStringCachedVec[it] == "T_t" ) continue;
     if( keyStringT.Contains("DYJets") || keyStringCachedVec[it] == "T_t" || keyStringT.Contains("Data") ) continue;

     cntPadIdx ++;
     if( cntPadIdx >= lastPadIdx ){ /*std::cout<<"Overwritten the last pad with index : "<<lastPadIdx<<"! Please make sure enough pads are created!"<<std::endl;*/ cntPadIdx = lastPadIdx; }

     cs->cd(cntPadIdx); TPad * pad = (TPad*) cs->GetPad(cntPadIdx); pad->Clear(); pad->SetLogy(statusLogy);

     TH1D *tmpHist = (TH1D*) h1_inputVec[it]->Clone();

     tmpHist->Rebin(nRebin); tmpHist->Scale(scaleMCCachedVec[it]);

     tmpHist->SetLineColor(colorCachedVec[it]); tmpHist->SetMarkerColor(colorCachedVec[it]); tmpHist->SetMarkerStyle(6); if( !drawText) tmpHist->SetMarkerSize(0.25);

     drawOverFlowBin(tmpHist);

     if( keyStringT.Contains("Data") ){ tmpHist->SetLineColor(kBlack); tmpHist->SetMarkerColor(kBlack); }

     if( !statusLogy ){
        tmpHist->SetAxisRange(0, tmpHist->GetMaximum()*1.5, "Y");
     }else{
        tmpHist->SetAxisRange(tmpHist->GetMaximum()/1e03, tmpHist->GetMaximum()*32*5, "Y");
     }

     if( keyStringT.Contains("Data") || keyStringT.Contains("Signal") ){
        if( !drawText ) tmpHist->Draw("e");
        else tmpHist->Draw("e text");
     }else{
        if( !drawText) tmpHist->Draw("hist");
        else tmpHist->Draw("hist text");
     }

     if( !keyStringT.Contains("Data") && !keyStringT.Contains("Signal") ){
        stackedKeyStringTVec.push_back(keyStringT);
        if( h1_stackedVec.empty() ){
           h1_stackedVec.push_back( (TH1D*) tmpHist->Clone() );
           h1_stackedVec.back()->SetFillColor( colorCachedVec[it] );
           h1_stackedVec.back()->SetLineColor( colorCachedVec[it] );
           h1_stackedVec.back()->SetMarkerColor( colorCachedVec[it] );
        }else{
           TH1D *tmpSumHist = (TH1D*) tmpHist->Clone();
           tmpSumHist->Add(h1_stackedVec.back(), 1);
           h1_stackedVec.push_back( (TH1D*) tmpSumHist->Clone() );
           h1_stackedVec.back()->SetFillColor( colorCachedVec[it] );
           h1_stackedVec.back()->SetLineColor( colorCachedVec[it] );
           h1_stackedVec.back()->SetMarkerColor( colorCachedVec[it] );
        }
     }
     if( keyStringT.Contains("Data") ){ h1_dataCached = (TH1D*) tmpHist->Clone(); }
     if( keyStringT.Contains("Signal") ){
        if( !h1_signal1Cached ) h1_signal1Cached = (TH1D*) tmpHist->Clone();
        else h1_signal2Cached = (TH1D*) tmpHist->Clone();
     }

  }

  for(int ic=cntPadIdx+1; ic<=lastPadIdx; ic++){ cs->cd(ic); TPad * pad = (TPad*) cs->GetPad(ic); pad->Clear(); pad->SetLogy(statusLogy); }

  Float_t legendX1 = .60;
  Float_t legendX2 = .85;
  Float_t legendY1 = .55;
  Float_t legendY2 = .85;
  TLegend* catLeg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);
  catLeg1->SetTextSize(0.030);

  cs->cd(lastPadIdx);

  if( h1_dataCached ){
     double dataIntegral = h1_dataCached->Integral();
     double sumTotIntegral = h1_stackedVec.back()->Integral();
     double normScale = dataIntegral/sumTotIntegral;

     h1_dataCached->Draw("e");
     catLeg1->AddEntry(h1_dataCached, "Data");
     for(unsigned int is=0; is<h1_stackedVec.size(); is++){
        unsigned int ivsIdx = h1_stackedVec.size()-is-1;
        TH1D * tmpStacked = (TH1D*) h1_stackedVec[ivsIdx]->Clone();
        if( doNormalization ) tmpStacked->Scale(normScale);
        tmpStacked->Draw("hist same");
        catLeg1->AddEntry(tmpStacked, stackedKeyStringTVec[ivsIdx]);
     }
     h1_signal1Cached->SetLineColor(kRed); h1_signal1Cached->SetMarkerColor(kRed); h1_signal1Cached->SetLineStyle(3);
     h1_signal1Cached->Draw("same");
     h1_dataCached->Draw("same e");

     catLeg1->AddEntry(h1_signal1Cached, "Signal");
     catLeg1->SetFillColor(kWhite);
     catLeg1->SetBorderSize(0);
     catLeg1->Draw();

     TPad *lastPad = (TPad*) cs->GetPad(lastPadIdx);
     lastPad->RedrawAxis();
  }

}

int main(){
   explore();
}

