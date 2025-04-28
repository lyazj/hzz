import awkward as ak
import numpy as np
import vector
import random
from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# Constants
m_Z = 91.2
m_W = 80.419
Gamma_Z = 2.494
cW = 0.8819
sW = 0.4713
gHZZ = m_Z / (cW * sW)
cL = -0.27775
cR = 0.22225

# Functions
def compute_invariant_mass(p1, p2):
    return (p1 + p2).mass

def K_factor(zmass1, zmass2, higgsmass):
    return (higgsmass**2 - zmass1**2 - zmass2**2) / (2 * zmass1 * zmass2)

def P_factor(mass):
    return (2 * gHZZ * mass**2) / (((mass**2 - m_Z**2)**2) + Gamma_Z**2 * m_Z**2)

def Asq_LL(theta1, phi1, theta2, phi2, P1, P2, K):
    return 4 * K**2 * P1 * P2 * (cL**2 + cR**2)**2 * (np.sin(theta1))**2 * (np.sin(theta2))**2

def Asq_pp(theta1, phi1, theta2, phi2, P1, P2, K):
    return P1 * P2 * (
        (cL**4) * (1+np.cos(theta1))**2 * (1+np.cos(theta2))**2 +
        (cR**4) * (1-np.cos(theta1))**2 * (1-np.cos(theta2))**2 +
        (cL**2 * cR**2) * ((1+np.cos(theta1))**2 * (1-np.cos(theta2))**2 + (1-np.cos(theta1))**2 * (1+np.cos(theta2))**2)
    )

def Asq_mm(theta1, phi1, theta2, phi2, P1, P2, K):
    return Asq_pp(theta2, phi2, theta1, phi1, P1, P2, K)

def Aint_LLpp(theta1, phi1, theta2, phi2, P1, P2, K):
    return -4 * K * P1 * P2 * (
        (cL**4) * (1+np.cos(theta1)) * (1+np.cos(theta2)) +
        (cR**4) * (1-np.cos(theta1)) * (1-np.cos(theta2)) -
        (cL**2 * cR**2) * ((1+np.cos(theta1)) * (1-np.cos(theta2)) + (1-np.cos(theta1)) * (1+np.cos(theta2)))
    ) * np.sin(theta1) * np.sin(theta2) * np.cos(phi1-phi2)

def Aint_LLmm(theta1, phi1, theta2, phi2, P1, P2, K):
    return Aint_LLpp(theta2, phi2, theta1, phi1, P1, P2, K)

def Aint_ppmm(theta1, phi1, theta2, phi2, P1, P2, K):
    return 2 * P1 * P2 * (cL**2 + cR**2)**2 * (np.sin(theta1) * np.sin(theta2))**2 * np.cos(2*(phi1-phi2))

class ReweightProcessor(processor.ProcessorABC):
    def process(self, events):
        # Replace custom 'Particle' branch with NanoAOD-compatible LHEPart branch
        particles = events.LHEPart

        # Select leptons (electrons/muons)
        electrons = particles[abs(particles.pdgId) == 11]
        muons = particles[abs(particles.pdgId) == 13]
        del particles
        mask = (ak.num(electrons) == 2) & (ak.num(muons) == 2)
        electrons = electrons[mask]
        muons = muons[mask]
        del mask

        # Build TLorentzVectors
        p4s_e = vector.zip({"pt": electrons.pt, "eta": electrons.eta, "phi": electrons.phi, "mass": 0*electrons.pt})
        p4s_mu = vector.zip({"pt": muons.pt, "eta": muons.eta, "phi": muons.phi, "mass": 0*muons.pt})
        e1, e2 = p4s_e[:,0], p4s_e[:,1]
        mu1, mu2 = p4s_mu[:,0], p4s_mu[:,1]

        z1 = e1 + e2
        z2 = mu1 + mu2
        higgs = z1 + z2

        # Boost to Higgs rest frame
        boostvec = -higgs.to_beta3()
        boostvec_4D = vector.zip({"x": boostvec.x, "y": boostvec.y, "z": boostvec.z, "t": 1.0})
        e1_boost = e1.boost_p4(boostvec_4D)
        mu1_boost = mu1.boost_p4(boostvec_4D)
        z1_boost = z1.boost_p4(boostvec_4D)
        z2_boost = z2.boost_p4(boostvec_4D)

        # Boost to Higgs rest frame
        boostvec_higgs = -higgs.to_beta3()
        boostvec_higgs_4D = vector.zip({"x": boostvec_higgs.x, "y": boostvec_higgs.y, "z": boostvec_higgs.z, "t": 1.0})
        e1_boost = e1.boost_p4(boostvec_higgs_4D)
        mu1_boost = mu1.boost_p4(boostvec_higgs_4D)
        z1_boost = z1.boost_p4(boostvec_higgs_4D)
        z2_boost = z2.boost_p4(boostvec_higgs_4D)
        
        # Boost each lepton into Z rest frames
        boostvec_z1 = -z1_boost.to_beta3()
        boostvec_z1_4D = vector.zip({"x": boostvec_z1.x, "y": boostvec_z1.y, "z": boostvec_z1.z, "t": 1.0})
        e1_final = e1_boost.boost_p4(boostvec_z1_4D)
        boostvec_z2 = -z2_boost.to_beta3()
        boostvec_z2_4D = vector.zip({"x": boostvec_z2.x, "y": boostvec_z2.y, "z": boostvec_z2.z, "t": 1.0})
        mu1_final = mu1_boost.boost_p4(boostvec_z2_4D)

        # Compute angles
        theta1 = z1_boost.deltaangle(e1_final)
        theta2 = (-z2_boost).deltaangle(mu1_final)

        # Define phi angles (azimuths)
        phi1 = e1_final.phi
        phi2 = mu1_final.phi

        # Masses
        zmass1 = z1.mass
        zmass2 = z2.mass
        higgsmass = higgs.mass

        # Compute weights
        P1 = P_factor(zmass1)
        P2 = P_factor(zmass2)
        K = K_factor(zmass1, zmass2, higgsmass)

        asq_LL = Asq_LL(theta1, phi1, theta2, phi2, P1, P2, K)
        asq_pp = Asq_pp(theta1, phi1, theta2, phi2, P1, P2, K)
        asq_mm = Asq_mm(theta1, phi1, theta2, phi2, P1, P2, K)
        aint_LLpp = Aint_LLpp(theta1, phi1, theta2, phi2, P1, P2, K)
        aint_LLmm = Aint_LLmm(theta1, phi1, theta2, phi2, P1, P2, K)
        aint_ppmm = Aint_ppmm(theta1, phi1, theta2, phi2, P1, P2, K)

        asq_tot = asq_LL + asq_pp + asq_mm + aint_LLpp + aint_LLmm + aint_ppmm
        asq_TT = asq_pp + asq_mm + aint_ppmm

        w_T = asq_TT / asq_tot
        w_L = asq_LL / asq_tot
        inter = (aint_LLpp + aint_LLmm) / asq_tot

        return {
            "w_T": w_T,
            "w_L": w_L,
            "inter": inter,
            "zmass1": zmass1,
            "zmass2": zmass2,
            "higgsmass": higgsmass,
            "costheta1": np.cos(theta1),
            "costheta_scatter": np.cos(higgs.theta)
        }

    def postprocess(self, accumulator):
        return accumulator
