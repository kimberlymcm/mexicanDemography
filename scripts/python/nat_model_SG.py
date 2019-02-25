"""
   Defined demographic models to use with dadi
"""


from dadi import Numerics, PhiManip, Integration, Spectrum
import dadi


def Kimberly_method_branching((Tb, Nb, T_s1, T_s2, T_s3, N_Tar, N_Hui, N_Trq, N_Mya),
                              (ns_TarMya, ns_TarTrq, ns_HuiMya, ns_HuiTrq, ns_TarHui, ns_MyaTrq),
                              pts):
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, Tb, nu=Nb)  #bottleneck shared by all 4 populations
    #First the TAR breaks off. 
    phi_T_HMQ = dadi.PhiManip.phi_1D_to_2D(xx, phi) #Splits the grid into 2 populations, this will be the oldest split: Tar on one side, HUI, MYA, TRQ on the other.
    phi_T_HMQ = dadi.Integration.two_pops(phi_T_HMQ, xx, T_s1, nu1=N_Tar, nu2=Nb, m12=0, m21=0)
    #This Phi is then [Tar, HuiTrqMya]
    #This is TAR vs. HUI (the next one to branch)
    phi_T_H = dadi.Integration.two_pops(phi_T_HMQ, xx, T_s2 + T_s3, nu1=N_Tar, nu2=N_Hui, m12=0, m21=0) #Since assuming TarHui comes first
    fs_Tar_Hui = dadi.Spectrum.from_phi(phi_T_H, ns_TarHui, (xx, xx)) #1
    #This is TAR vs. MYA
    phi_T_M = dadi.Integration.two_pops(phi_T_HMQ, xx, T_s2, nu1=N_Tar, nu2=Nb, m12=0, m21=0) #Since assuming TarHui comes first
    phi_T_M = dadi.Integration.two_pops(phi_T_M, xx, T_s3, nu1=N_Tar, nu2=N_Mya, m12=0, m21=0) #Since assuming TarHui comes first
    fs_Tar_Mya = dadi.Spectrum.from_phi(phi_T_M, ns_TarMya, (xx, xx)) #2
    #TAR vs. TRQ
    phi_T_Q = dadi.Integration.two_pops(phi_T_HMQ, xx, T_s2, nu1=N_Tar, nu2=Nb, m12=0, m21=0) #Since assuming TarHui comes first
    phi_T_Q = dadi.Integration.two_pops(phi_T_Q, xx, T_s2, nu1=N_Tar, nu2=N_Trq, m12=0, m21=0) #Since assuming TarHui comes first
    fs_Tar_Trq = dadi.Spectrum.from_phi(phi_T_Q, ns_TarTrq, (xx, xx)) #3
    #Now next split (HUI, (MYA,TRQ))
    phi_ONE = dadi.Integration.one_pop(phi, xx, T_s1, nu=Nb) #Getting the pop freq distribution after the first split for TarHui
    phi_H_MQ = dadi.PhiManip.phi_1D_to_2D(xx, phi_ONE) 
    #HUI vs. MYA
    phi_H_M = dadi.Integration.two_pops(phi_H_MQ, xx, T_s2, nu1=N_Hui, nu2=Nb) 
    phi_H_M = dadi.Integration.two_pops(phi_H_M, xx, T_s3, nu1=N_Hui, nu2=N_Mya) 
    fs_Hui_Mya = dadi.Spectrum.from_phi(phi_H_M, ns_HuiMya, (xx, xx)) #4
    #HUI vs TRQ
    phi_H_Q = dadi.Integration.two_pops(phi_H_MQ, xx, T_s3, nu1=N_Hui, nu2=N_Trq) 
    fs_Hui_Trq = dadi.Spectrum.from_phi(phi_H_Q, ns_HuiTrq, (xx, xx)) #5

    #Now the MYA vs TRQ splot
    phi_TWO = dadi.Integration.one_pop(phi, xx, T_s1 + T_s2, nu=Nb) #Getting the pop freq distribution after the first split for TarHui
    phi_M_Q = dadi.PhiManip.phi_1D_to_2D(xx, phi_TWO)
    phi_M_Q = dadi.Integration.two_pops(phi_M_Q, xx, T_s3, nu1=N_Mya, nu2=N_Trq) 
    fs_Mya_Trq = dadi.Spectrum.from_phi(phi_M_Q, ns_MyaTrq, (xx, xx)) #5
    #return [fs_Tar_Mya.swapaxes(0,1), fs_Hui_Mya, fs_Tar_Trq, fs_Hui_Trq.swapaxes(0,1), fs_MQ, fs_TH] #Flipped if T_TarHui comes second
    return [fs_Tar_Mya, fs_Tar_Trq, fs_Hui_Mya, fs_Hui_Trq, fs_Tar_Hui, fs_Mya_Trq] #fs_TN_Q=[TarHui,Mya], fs_H_QM=[TarHui,Trq], fs_TH=[Tar,Hui], fs_MQ=[Mya,Trq]


def Kimberly_method_fewerSizes((Tb, Nb, T_s1, T_TarHui, T_TrqMya, N_Tar, N_Hui, N_Trq, N_Mya),
                               (ns_TarMya, ns_TarTrq, ns_HuiMya, ns_HuiTrq, ns_TarHui, ns_MyaTrq),
                               pts):
    """
        Model:
         Tb: Time between bottleneck shared by all 4 populations and the deepest split (TarHui vs. TrqMya)
         Nb: Bottleneck size shared by all 4 populations.
             We assume this is the population size of each set of pops until each population is on its own
         T_s1: Split time between first deepest split (TarHui vs. TrqMya) and 2nd split (Tar vs. Hui)
         T_TarHui: Split time between the 2nd split and 3rd split (Trq  and Mya split)
         T_TrqMya: Split time between 3rd split and present
         N_Tar: Size of TAR after all the splits
         N_Hui: Size of HUI after all the splits
         N_Trq: Size of MTRQ after all the splits
         N_Mya: Size of MYA after all the splits  
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, Tb, nu=Nb) # Bottleneck shared by all 4 populations

    # Below splits the grid into 2 populations, this will be the oldest split: Tar & Hui on one side, Trq & Mya on the other
    phi_TH_QM = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # T_s1 is the time between the first, deepest split and second split. So integrate the two new populations for that time.
    phi_TH_QM = dadi.Integration.two_pops(phi_TH_QM, xx, T_s1, nu1=Nb, nu2=Nb, m12=0, m21=0)
    # This Phi is then [TarHui, TrqMya], this brings us to T_TarHui, which we assume comes before T_MyaTrq.

    # For T_TarHui time, we have 3 populations, Tar, Hui, and MyaTrq. Then for T_TrqMya time we have 4 populations: Tar, Hui, Mya, Trq
    phi_T_QM = dadi.Integration.two_pops(phi_TH_QM, xx, T_TarHui, nu1=N_Tar, nu2=Nb, m12=0, m21=0) # Tar vs. MyaTrq for T_TarHui time
    phi_T_M = dadi.Integration.two_pops(phi_T_QM, xx, T_TrqMya, nu1=N_Tar, nu2=N_Mya, m12=0, m21=0) # Tar vs. Mya for T_TrqMya time (to present)
    fs_Tar_Mya = dadi.Spectrum.from_phi(phi_T_M, ns_TarMya, (xx, xx))
    phi_T_M = dadi.Integration.two_pops(phi_T_QM, xx, T_TrqMya, nu1=N_Tar, nu2=N_Trq, m12=0, m21=0) # Tar vs. Trq for T_TrqMya time (to present)
    fs_Tar_Trq = dadi.Spectrum.from_phi(phi_T_M, ns_TarTrq, (xx, xx))

    phi_H_QM = dadi.Integration.two_pops(phi_TH_QM, xx, T_TarHui, nu1=N_Hui, nu2=Nb, m12=0, m21=0) # Hui vs MyaTrq for T_TarHui time
    phi_H_M = dadi.Integration.two_pops(phi_H_QM, xx, T_TrqMya, nu1=N_Hui, nu2=N_Mya, m12=0, m21=0) # Hui vs. Mya for T_TrqMya time (to present)
    fs_Hui_Mya = dadi.Spectrum.from_phi(phi_T_M, ns_HuiMya, (xx, xx))
    phi_H_Q = dadi.Integration.two_pops(phi_H_QM, xx, T_TrqMya, nu1=N_Hui, nu2=N_Trq, m12=0, m21=0) # Hui vs. Trq for T_TrqMya time (to present)
    fs_Hui_Trq = dadi.Spectrum.from_phi(phi_H_Q, ns_HuiTrq, (xx, xx))
    
    # Getting the pop freq distribution after the first split for TarHui
    # Integrates forward T_s1 in time, at which point Tar and Hui split
    phi_ONE = dadi.Integration.one_pop(phi, xx, T_s1, nu=Nb)
    
    # Getting the pop freq distribution after the first split for MyaTrq.
    # Integrates forward T_s1 + T_TarHui in time, at which point Mya and Trq split
    phi_TWO = dadi.Integration.one_pop(phi, xx, T_s1 + T_TarHui, nu=Nb)

    # Now I want to do the last two splits - lets say that the TarHui split comes first.
    phi_TH = dadi.PhiManip.phi_1D_to_2D(xx, phi_ONE) # Splits TarHui at Tb + T_s1
    phi_TH = dadi.Integration.two_pops(phi_TH, xx, T_TarHui + T_TrqMya, nu1=N_Tar, nu2=N_Hui) # Integrates freq spec to present data
    fs_TH = dadi.Spectrum.from_phi(phi_TH, ns_TarHui, (xx, xx)) #fs [Tar, Hui]

    # Now lets say the MyaTrq split came second
    phi_MQ = dadi.PhiManip.phi_1D_to_2D(xx, phi_TWO) 
    phi_MQ = dadi.Integration.two_pops(phi_MQ, xx, T_TrqMya, nu1=N_Mya, nu2=N_Trq) # Integrates freq spec to present data
    fs_MQ = dadi.Spectrum.from_phi(phi_MQ, ns_MyaTrq, (xx, xx)) #fs [Mya, Trq]
    #return [fs_Tar_Mya.swapaxes(0,1), fs_Hui_Mya, fs_Tar_Trq, fs_Hui_Trq.swapaxes(0,1), fs_MQ, fs_TH] #Flipped if T_TarHui comes second
    return [fs_Tar_Mya, fs_Tar_Trq, fs_Hui_Mya, fs_Hui_Trq, fs_TH, fs_MQ] #fs_TN_Q=[TarHui,Mya], fs_H_QM=[TarHui,Trq], fs_TH=[Tar,Hui], fs_MQ=[Mya,Trq]


def Kimberly_convert_params_fewerSizes((Tb, Nb, T_s1, T_TarHui, T_TrqMya, N_Tar, N_Hui, N_Trq, N_Mya),
                                       theta=3813.75, L=3101932, mu=2.36e-8, gen_time=25.):
    """
    Converts parameters from genetic to physical units. The default theta gives Na=10k.
    """
    # L is the 'effective' length accounting for scanned length and chimp
    # success rate.
    theta_per_base = 1.*theta/L
    Na = theta_per_base/(4*mu)
    Tb = Tb*2*Na*gen_time
    T_s1 = T_s1 * 2 * Na * gen_time
    T_TarHui = T_TarHui*2*Na*gen_time
    T_TrqMya = T_TrqMya*2*Na*gen_time
    # Convert the times given as fractional:
    Nb *= Na
    N_Tar *= Na
    N_Hui *= Na
    N_Trq *= Na
    N_Mya *= Na
    return (Na, Tb+T_s1+T_TarHui+T_TrqMya, T_s1+T_TarHui + T_TrqMya, T_TarHui+T_TrqMya, T_TrqMya, Nb, N_Tar, N_Hui, N_Trq, N_Mya)
