#include <ROOT/RDataFrame.hxx>
#include <cmath>
#include <cstdlib>
#include <random>

#include "ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TVector3.h"

using namespace std;

const double m_Z = 91.2;
const double m_W = 80.419;
const double Gamma_Z = 2.494;
const double cW = 0.8819;
const double sW = 0.4713;
const double gHZZ = m_Z / (cW * sW);
const double cL = -0.27775;
const double cR = 0.22225;  // ?

double ComputeInvariantMass(
    float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2)
{
  TLorentzVector p1;
  TLorentzVector p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, mass1);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, mass2);
  return (p1 + p2).M();
}

double K_factor(float zmass1, float zmass2, float higgsmass)
{
  return (higgsmass * higgsmass - zmass1 * zmass1 - zmass2 * zmass2) / (2 * zmass1 * zmass2);
}

double P_factor(float mass)
{
  return (2 * gHZZ * mass * mass)
      / ((mass * mass - m_Z * m_Z) * (mass * mass - m_Z * m_Z) + Gamma_Z * Gamma_Z * m_Z * m_Z);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

double Asq_LL(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return 4 * K * K * P1 * P2 * (cL * cL + cR * cR) * (cL * cL + cR * cR) * (sin(theta1) * sin(theta1))
      * (sin(theta2) * sin(theta2));
}

double Asq_pp(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return P1 * P2
      * (pow(cL, 4) * pow((1 + cos(theta1)), 2) * pow((1 + cos(theta2)), 2)
          + pow(cR, 4) * pow((1 - cos(theta1)), 2) * pow((1 - cos(theta2)), 2)
          + pow(cL, 2) * pow(cR, 2) * pow((1 + cos(theta1)), 2) * pow((1 - cos(theta2)), 2)
          + cL * cL * cR * cR * pow((1 - cos(theta1)), 2) * pow((1 + cos(theta2)), 2));
}

double Asq_mm(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return P1 * P2
      * (pow(cL, 4) * pow((1 - cos(theta1)), 2) * pow((1 - cos(theta2)), 2)
          + pow(cR, 4) * pow((1 + cos(theta1)), 2) * pow((1 + cos(theta2)), 2)
          + pow(cL, 2) * pow(cR, 2) * pow((1 + cos(theta1)), 2) * pow((1 - cos(theta2)), 2)
          + cL * cL * cR * cR * pow((1 - cos(theta1)), 2) * pow((1 + cos(theta2)), 2));
}

double Aint_LLpp(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return -4 * K * P1 * P2
      * (pow(cL, 4) * (1 + cos(theta1)) * (1 + cos(theta2)) + pow(cR, 4) * (1 - cos(theta1)) * (1 - cos(theta2))
          - cL * cL * cR * cR * (1 + cos(theta1)) * (1 - cos(theta2))
          - cL * cL * cR * cR * (1 - cos(theta1)) * (1 + cos(theta2)))
      * sin(theta1) * sin(theta2) * cos(phi1 - phi2);
}

double Aint_LLmm(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return -4 * K * P1 * P2
      * (pow(cL, 4) * (1 - cos(theta1)) * (1 - cos(theta2)) + pow(cR, 4) * (1 + cos(theta1)) * (1 + cos(theta2))
          - cL * cL * cR * cR * (1 + cos(theta1)) * (1 - cos(theta2))
          - cL * cL * cR * cR * (1 - cos(theta1)) * (1 + cos(theta2)))
      * sin(theta1) * sin(theta2) * cos(phi1 - phi2);
}

double Aint_ppmm(float theta1, float phi1, float theta2, float phi2, double P1, double P2, double K)
{
  return 2 * P1 * P2 * pow((cL * cL + cR * cR), 2) * pow((sin(theta1) * sin(theta2)), 2) * cos(2 * (phi1 - phi2));
}

#pragma GCC diagnostic pop

int Reweight()
{
  gSystem->Load("ExRootAnalysis/libExRootAnalysis.so");
  //gSystem->Load("libPhysics");
  TString filepath = "root/";
  TString filename = "h2zz.root";
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df_sig("LHEF", filepath + filename);

  auto df_sig_boost =
      df_sig
          .Define("angle",
              [](const ROOT::RVec<double> &pt, const ROOT::RVec<double> &eta, const ROOT::RVec<double> &phi,
                  const ROOT::RVec<int> &pid, const ROOT::RVec<int> &mother) {
                // Select the id of leptons.
                double z1_l1_id = -1, z1_l2_id = -1, z2_l1_id = -1, z2_l2_id = -1;
                std::vector<int> id;
                for(int i = 2; i < (int)pid.size(); i++) {
                  if(fabs(pid[i]) != 11 && fabs(pid[i]) != 13) continue;
                  id.push_back(i);
                }
                if(id.size() != 4) { cout << "\t!!!Error: Incorrect particles number!!!" << endl; }
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(0, 3);
                int id_dominate = id[dis(gen)];
                z1_l1_id = id_dominate;
                //z1_l1_id = 2;
                for(int i = 0; i < (int)id.size(); i++) {
                  if(id[i] == z1_l1_id) continue;
                  if(mother[id[i]] == mother[z1_l1_id]) { z1_l2_id = id[i]; }
                  if(mother[id[i]] != mother[z1_l1_id]) {
                    if(z2_l1_id == -1)
                      z2_l1_id = id[i];
                    else if(z2_l2_id == -1)
                      z2_l2_id = id[i];
                    else
                      cout << "\tError: invalid mother id!!" << endl;
                  }
                }
                if(z1_l1_id == -1 || z1_l2_id == -1 || z2_l1_id == -1 || z2_l2_id == -1)
                  cout << "\tError: Invalid particle id!!" << endl;
                //cout << "e1:" << z1_l1_id << "; e2:" << z1_l2_id << ";mu1: " << z2_l1_id << "; mu2:" << z2_l2_id
                //     << endl;

                // Z1
                TLorentzVector e1_p4;  // positive
                TLorentzVector e2_p4;  // negative
                if(pid[z1_l1_id] < 0) {
                  e1_p4.SetPtEtaPhiM(pt[z1_l1_id], eta[z1_l1_id], phi[z1_l1_id], 0.0);
                  e2_p4.SetPtEtaPhiM(pt[z1_l2_id], eta[z1_l2_id], phi[z1_l2_id], 0.0);
                } else {
                  e2_p4.SetPtEtaPhiM(pt[z1_l1_id], eta[z1_l1_id], phi[z1_l1_id], 0.0);
                  e1_p4.SetPtEtaPhiM(pt[z1_l2_id], eta[z1_l2_id], phi[z1_l2_id], 0.0);
                }
                TLorentzVector z1_p4 = e1_p4 + e2_p4;
                float THETA = z1_p4.Theta();
                float zmass1 = z1_p4.M();

                // Z2
                TLorentzVector mu1_p4;
                TLorentzVector mu2_p4;
                if(pid[z2_l1_id] < 0) {
                  mu1_p4.SetPtEtaPhiM(pt[z2_l1_id], eta[z2_l1_id], phi[z2_l1_id], 0.0);
                  mu2_p4.SetPtEtaPhiM(pt[z2_l2_id], eta[z2_l2_id], phi[z2_l2_id], 0.0);
                } else {
                  mu2_p4.SetPtEtaPhiM(pt[z2_l1_id], eta[z2_l1_id], phi[z2_l1_id], 0.0);
                  mu1_p4.SetPtEtaPhiM(pt[z2_l2_id], eta[z2_l2_id], phi[z2_l2_id], 0.0);
                }
                TLorentzVector z2_p4 = mu1_p4 + mu2_p4;
                float zmass2 = z2_p4.M();

                // Boosting
                // co
                TLorentzVector higgs_p4 = z1_p4 + z2_p4;
                float higgsmass = higgs_p4.M();

                TVector3 beta_higgs = -higgs_p4.BoostVector();
                z1_p4.Boost(beta_higgs);
                z2_p4.Boost(beta_higgs);
                TLorentzVector e1_p4_cms = e1_p4;
                TLorentzVector mu1_p4_cms = mu1_p4;
                e1_p4_cms.Boost(beta_higgs);
                mu1_p4_cms.Boost(beta_higgs);

                // Boosting
                TVector3 beta_z1 = -z1_p4.BoostVector();
                TVector3 beta_z2 = -z2_p4.BoostVector();
                e1_p4_cms.Boost(beta_z1);
                mu1_p4_cms.Boost(beta_z2);

                TVector3 v3_z1 = z1_p4.Vect();
                TVector3 e1_v3_cms = e1_p4_cms.Vect();
                TVector3 v3_z2 = -z2_p4.Vect();
                TVector3 mu1_v3_cms = mu1_p4_cms.Vect();
                TVector3 v3_zaxis_lab(0, 0, 1);
                TVector3 v3_zaxis_lab2(0, 0, 1);  // Will not affect the results since phi is uniform?

                // Calculate the theta and phi.
                // Z1
                TVector3 v3_n_1 = v3_zaxis_lab.Cross(v3_z1);
                double angle_zaxis_to_z1 = v3_zaxis_lab.Angle(v3_z1);
                v3_n_1 *= (1.0 / sin(angle_zaxis_to_z1));
                TVector3 unit_v3_z1 = v3_z1.Unit();
                TVector3 ProjectTovBoson_e1_z1 = unit_v3_z1;
                ProjectTovBoson_e1_z1 *= (unit_v3_z1.Dot(e1_v3_cms));
                TVector3 PerpendicularComponent_e1_z1 = (e1_v3_cms - ProjectTovBoson_e1_z1);
                float phi1 = v3_n_1.Angle(PerpendicularComponent_e1_z1);
                TVector3 v3_r1 = (1.0 / sin(angle_zaxis_to_z1)) * (v3_zaxis_lab - (cos(angle_zaxis_to_z1)) * v3_z1);
                if((v3_r1.Dot(PerpendicularComponent_e1_z1) < 0)) phi1 *= -1.0;
                float theta1 = v3_z1.Angle(e1_v3_cms);
                // Z2
                TVector3 v3_n_2 = v3_zaxis_lab2.Cross(v3_z2);
                double angle_zaxis_to_z2 = v3_zaxis_lab2.Angle(v3_z2);
                v3_n_2 *= (1.0 / sin(angle_zaxis_to_z2));
                TVector3 unit_v3_z2 = v3_z2.Unit();
                TVector3 ProjectTovBoson_mu1_z2 = unit_v3_z2;
                ProjectTovBoson_mu1_z2 *= (unit_v3_z2.Dot(mu1_v3_cms));
                TVector3 PerpendicularComponent_mu1_z2 = (mu1_v3_cms - ProjectTovBoson_mu1_z2);
                float phi2 = v3_n_2.Angle(PerpendicularComponent_mu1_z2);
                TVector3 v3_r2 = (1.0 / sin(angle_zaxis_to_z2)) * (v3_zaxis_lab2 - (cos(angle_zaxis_to_z2)) * v3_z2);
                if((v3_r2.Dot(PerpendicularComponent_mu1_z2) < 0)) phi2 *= -1.0;
                float theta2 = v3_z2.Angle(mu1_v3_cms);

                return ROOT::RVec<float>{ theta1, phi1, theta2, phi2, THETA, zmass1, zmass2, higgsmass };
              },
              { "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.PID", "Particle.Mother1" })

          .Define("costheta_scatter",
              [](ROOT::RVec<float> &angle) {
                float theta = angle[4];
                return cos(theta);
              },
              { "angle" })

          .Define("costheta1",
              [](ROOT::RVec<float> &angle) {
                float theta = angle[0];
                return cos(theta);
              },
              { "angle" })

          .Define("zmass1", [](ROOT::RVec<float> &angle) { return angle[5]; }, { "angle" })

          .Define("zmass2", [](ROOT::RVec<float> &angle) { return angle[6]; }, { "angle" })

          .Define("delta_phi",
              [](ROOT::RVec<float> &angle) {
                float phi1 = angle[1];
                float phi2 = angle[3];
                return phi1 - phi2;
              },
              { "angle" })

          .Define("higgsmass", [](ROOT::RVec<float> &angle) { return angle[7]; }, { "angle" });

  auto df_sig_weight = df_sig_boost
                           .Define("w_T",
                               [](ROOT::RVec<float> &angle) {
                                 float theta1 = angle[0];
                                 float phi1 = angle[1];
                                 float theta2 = angle[2];
                                 float phi2 = angle[3];
                                 float zmass1 = angle[5];
                                 float zmass2 = angle[6];
                                 float higgsmass = angle[7];
                                 double P1 = P_factor(zmass1);
                                 double P2 = P_factor(zmass2);
                                 double K = K_factor(zmass1, zmass2, higgsmass);
                                 double asq_LL = Asq_LL(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_pp = Asq_pp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_mm = Asq_mm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLpp = Aint_LLpp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLmm = Aint_LLmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_ppmm = Aint_ppmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_tot = asq_LL + asq_mm + asq_pp + aint_LLmm + aint_LLpp + aint_ppmm;
                                 double asq_TT = asq_pp + asq_mm + aint_ppmm;
                                 return asq_TT / asq_tot;
                               },
                               { "angle" })

                           .Define("w_L",
                               [](ROOT::RVec<float> &angle) {
                                 float theta1 = angle[0];
                                 float phi1 = angle[1];
                                 float theta2 = angle[2];
                                 float phi2 = angle[3];
                                 float zmass1 = angle[5];
                                 float zmass2 = angle[6];
                                 float higgsmass = angle[7];
                                 double P1 = P_factor(zmass1);
                                 double P2 = P_factor(zmass2);
                                 double K = K_factor(zmass1, zmass2, higgsmass);
                                 double asq_LL = Asq_LL(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_pp = Asq_pp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_mm = Asq_mm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLpp = Aint_LLpp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLmm = Aint_LLmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_ppmm = Aint_ppmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_tot = asq_LL + asq_mm + asq_pp + aint_LLmm + aint_LLpp + aint_ppmm;
                                 return asq_LL / asq_tot;
                               },
                               { "angle" })

                           .Define("inter",
                               [](ROOT::RVec<float> &angle) {
                                 float theta1 = angle[0];
                                 float phi1 = angle[1];
                                 float theta2 = angle[2];
                                 float phi2 = angle[3];
                                 float zmass1 = angle[5];
                                 float zmass2 = angle[6];
                                 float higgsmass = angle[7];
                                 double P1 = P_factor(zmass1);
                                 double P2 = P_factor(zmass2);
                                 double K = K_factor(zmass1, zmass2, higgsmass);
                                 double asq_LL = Asq_LL(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_pp = Asq_pp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_mm = Asq_mm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLpp = Aint_LLpp(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_LLmm = Aint_LLmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double aint_ppmm = Aint_ppmm(theta1, phi1, theta2, phi2, P1, P2, K);
                                 double asq_tot = asq_LL + asq_mm + asq_pp + aint_LLmm + aint_LLpp + aint_ppmm;
                                 return (aint_LLpp + aint_LLmm) / asq_tot;
                               },
                               { "angle" });

  df_sig_weight.Snapshot("tree", "./test_" + filename);

  return 0;
}
