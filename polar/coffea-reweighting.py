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
        # Assume 'Particle' branches (custom LHE-like)
        particles = events.Particle

        # Select leptons (electrons/muons)
        is_lepton = (abs(particles.PID) == 11) | (abs(particles.PID) == 13)
        leptons = particles[is_lepton]

        # Random choice of dominant lepton
        indices = ak.local_index(leptons.pt)
        rand_idx = ak.Array(random.choices(range(4), k=len(leptons)))
        dominate = leptons[rand_idx]

        # Pairing
        mothers = leptons.Mother1
        def pairing(leptons, dominate_idx):
            dominate_mother = leptons[dominate_idx].Mother1
            ids = ak.local_index(leptons)
            z1 = ids[leptons.Mother1 == dominate_mother]
            z2 = ids[leptons.Mother1 != dominate_mother]
            return (z1[0], z1[1], z2[0], z2[1])

        z1_l1_idx, z1_l2_idx, z2_l1_idx, z2_l2_idx = ak.unzip([pairing(leptons, i) for i in rand_idx])

        # Build TLorentzVectors
        lep = vector.awkward.zip({"pt": leptons.pt, "eta": leptons.eta, "phi": leptons.phi, "mass": 0*leptons.pt})
        e1, e2, mu1, mu2 = lep[z1_l1_idx], lep[z1_l2_idx], lep[z2_l1_idx], lep[z2_l2_idx]

        z1 = e1 + e2
        z2 = mu1 + mu2
        higgs = z1 + z2

        # Boost to Higgs rest frame
        boostvec = -higgs.to_beta3()
        e1_boost = e1.boost_p4(boostvec)
        mu1_boost = mu1.boost_p4(boostvec)

        z1_boost = z1.boost_p4(boostvec)
        z2_boost = z2.boost_p4(boostvec)

        # Boost each lepton into Z rest frames
        e1_final = e1_boost.boost_p4(-z1_boost.to_beta3())
        mu1_final = mu1_boost.boost_p4(-z2_boost.to_beta3())

        # Compute angles
        theta1 = z1_boost.angle(e1_final)
        theta2 = (-z2_boost).angle(mu1_final)

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

events = NanoEventsFactory.from_root('forran_additional.root', treepath='LHEF').events()
output = ReweightProcessor().process(events)
