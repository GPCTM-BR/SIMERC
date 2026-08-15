# cython: embedsignature=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: language = c++

import cython
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import AbstractState
from libc.math cimport sqrt, fabs, M_PI as PI

# Mapeamento estático das constantes de entrada para evitar buscas no Python
cdef int HmassP_INPUTS = CP.HmassP_INPUTS
cdef int PSmass_INPUTS = CP.PSmass_INPUTS
cdef int PQ_INPUTS = CP.PQ_INPUTS
cdef int PT_INPUTS = CP.PT_INPUTS
cdef int iphase_liquid = CP.iphase_liquid
cdef int iphase_gas = CP.iphase_gas
cdef int iphase_twophase = CP.iphase_twophase

# =============================================================================
# CLASS: MATERIAL STREAM
# =============================================================================
cdef class MaterialStream:
    cdef public object fluid
    cdef public double h, s, T, p, Q, rho, mass_flow

    cdef object _cp_update
    cdef object _cp_hmass
    cdef object _cp_smass
    cdef object _cp_T
    cdef object _cp_p
    cdef object _cp_Q
    cdef object _cp_rhomass

    def __init__(self, object fluid, bint bicubic=False):
        # Substitui os None por 0.0, pro C não reclamar
        self.mass_flow = 1.0
        self.h = 0.0
        self.s = 0.0
        self.T = 0.0
        self.p = 0.0
        self.Q = 0.0
        self.rho = 0.0
        
        # igual na versão pyhton pura, pode ser str ou AbstractState
        if isinstance(fluid, str):
            try:
                if bicubic:
                    self.fluid = AbstractState("BICUBIC&HEOS", fluid)
                else:
                    self.fluid = AbstractState("HEOS", fluid)
            except Exception:
                raise ValueError(f"Error: {fluid} is not a valid fluid name in CoolProp.")
        elif isinstance(fluid, AbstractState):
            self.fluid = fluid
        else:
            raise ValueError("Fluid must be a string or an AbstractState object.")

        self._cp_update = self.fluid.update
        self._cp_hmass = self.fluid.hmass
        self._cp_smass = self.fluid.smass
        self._cp_T = self.fluid.T
        self._cp_p = self.fluid.p
        self._cp_Q = self.fluid.Q
        self._cp_rhomass = self.fluid.rhomass

    cpdef void calculate(self, int flash_type, double property_1, double property_2) except *:
        """
        Calculate the thermodynamic state using two properties and the CoolProp constant.
        """
        self._cp_update(flash_type, property_1, property_2)
        self.h = self._cp_hmass()
        self.s = self._cp_smass()
        self.T = self._cp_T()
        self.p = self._cp_p()
        self.Q = self._cp_Q()
        self.rho = self._cp_rhomass()


# =============================================================================
# CLASS: PUMP
# =============================================================================
cdef class Pump:
    cdef public MaterialStream inlet_stream
    cdef public MaterialStream outlet_stream
    cdef public object fluid
    cdef public double h_1, s_1, mass_flow, isentropic_efficiency
    cdef public double h_2_is, h_2, s_2, work

    # DECLARAÇÃO DOS ATALHOS
    cdef object _cp_update
    cdef object _cp_hmass
    cdef object _cp_smass
    cdef object _cp_T
    cdef object _cp_Q
    cdef object _cp_rhomass

    def __init__(self, MaterialStream inlet_stream, MaterialStream outlet_stream, object fluid):
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.h_1 = inlet_stream.h
        self.s_1 = inlet_stream.s
        self.mass_flow = inlet_stream.mass_flow 
        self.isentropic_efficiency = 0.75

        # 2. CAPTURA DOS ATALHOS NA INICIALIZAÇÃO
        self._cp_update = self.fluid.update
        self._cp_hmass = self.fluid.hmass
        self._cp_smass = self.fluid.smass
        self._cp_T = self.fluid.T
        self._cp_Q = self.fluid.Q
        self._cp_rhomass = self.fluid.rhomass

    def set_isentropic_efficiency(self, double isentropic_efficiency, object percent=False):
        if percent is True or percent == '%':
            self.isentropic_efficiency = isentropic_efficiency / 100.0
        else:
            self.isentropic_efficiency = isentropic_efficiency

    cpdef void calculate(self, double outlet_pressure) except *:
        self.h_1 = self.inlet_stream.h
        self.s_1 = self.inlet_stream.s
        self.mass_flow = self.inlet_stream.mass_flow

        
        self._cp_update(PSmass_INPUTS, outlet_pressure, self.s_1)
        self.h_2_is = self._cp_hmass()
        self.h_2 = self.h_1 + (self.h_2_is - self.h_1) / self.isentropic_efficiency
        
        self._cp_update(HmassP_INPUTS, self.h_2, outlet_pressure)
        self.s_2 = self._cp_smass()
        self.work = self.mass_flow * (self.h_2 - self.h_1)

        # Atualiza a corrente de saída
        self.outlet_stream.mass_flow = self.mass_flow
        self.outlet_stream.p = outlet_pressure
        self.outlet_stream.T = self._cp_T()
        self.outlet_stream.h = self.h_2
        self.outlet_stream.s = self.s_2
        self.outlet_stream.Q = self._cp_Q()
        self.outlet_stream.rho = self._cp_rhomass()


# =============================================================================
# CLASS: VALVE
# =============================================================================
cdef class Valve:
    cdef public MaterialStream inlet_stream
    cdef public MaterialStream outlet_stream
    cdef public object fluid
    cdef public double h_1, P_1, P_2

    # DECLARAÇÃO DOS ATALHOS
    cdef object _cp_update
    cdef object _cp_T
    cdef object _cp_hmass
    cdef object _cp_smass
    cdef object _cp_Q
    cdef object _cp_rhomass

    def __init__(self, MaterialStream inlet_stream, MaterialStream outlet_stream, object fluid):
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.h_1 = inlet_stream.h
        self.P_1 = inlet_stream.p

        # 2. CAPTURA DOS ATALHOS LENTA (PROCURA APENAS UMA VEZ)
        self._cp_update = self.fluid.update
        self._cp_T = self.fluid.T
        self._cp_hmass = self.fluid.hmass
        self._cp_smass = self.fluid.smass
        self._cp_Q = self.fluid.Q
        self._cp_rhomass = self.fluid.rhomass

    cpdef void calculate(self, double outlet_pressure) except *:
        
        self.h_1 = self.inlet_stream.h
        self.P_1 = self.inlet_stream.p
        self.P_2 = outlet_pressure
        
        
        self._cp_update(HmassP_INPUTS, self.h_1, self.P_2)

        self.outlet_stream.mass_flow = self.inlet_stream.mass_flow
        self.outlet_stream.p = outlet_pressure
        self.outlet_stream.T = self._cp_T()
        self.outlet_stream.h = self._cp_hmass()
        self.outlet_stream.s = self._cp_smass()
        self.outlet_stream.Q = self._cp_Q()
        self.outlet_stream.rho = self._cp_rhomass()


# =============================================================================
# CLASS: HEATER / COOLER
# =============================================================================
cdef class HeaterCooler:
    cdef public MaterialStream inlet_stream
    cdef public MaterialStream outlet_stream
    cdef public object fluid
    cdef public double efficiency, mass_flow, h_1, P_1, heat_duty

    cdef object _cp_update
    cdef object _cp_hmass
    cdef object _cp_smass
    cdef object _cp_rhomass
    cdef object _cp_T
    cdef object _cp_Q
    cdef object _cp_p

    def __init__(self, MaterialStream inlet_stream, MaterialStream outlet_stream, object fluid, double efficiency=1.0):
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.efficiency = efficiency
        self.mass_flow = inlet_stream.mass_flow
        self.h_1 = inlet_stream.h
        self.P_1 = inlet_stream.p

        self._cp_update = self.fluid.update
        self._cp_hmass = self.fluid.hmass
        self._cp_smass = self.fluid.smass
        self._cp_rhomass = self.fluid.rhomass
        self._cp_T = self.fluid.T
        self._cp_Q = self.fluid.Q
        self._cp_p = self.fluid.p

    def set_efficiency(self, double efficiency, object percent=False):
        if percent is True or percent == '%':
            self.efficiency = efficiency / 100.0
        else:
            self.efficiency = efficiency

    cpdef void calculate(self, int flash_type, double value_1, double value_2) except *:
        self.mass_flow = self.inlet_stream.mass_flow
        self.h_1 = self.inlet_stream.h
        self.P_1 = self.inlet_stream.p

        self._cp_update(flash_type, value_1, value_2)
        self.heat_duty = self.mass_flow * (self._cp_hmass() - self.h_1) / self.efficiency
        
        self.outlet_stream.mass_flow = self.mass_flow
        self.outlet_stream.p = self._cp_p()
        self.outlet_stream.T = self._cp_T()
        self.outlet_stream.h = self._cp_hmass()
        self.outlet_stream.s = self._cp_smass()
        self.outlet_stream.Q = self._cp_Q()
        self.outlet_stream.rho = self._cp_rhomass()

# =============================================================================
# CLASS: EJECTOR
# =============================================================================
cdef class Ejector:
    """
    A CPM ejector model. It uses the modeling developed by Cardemil and Cole (2012).
    The user must provide the inlet streams (primary and secondary) and the fluid.
    The fluid must be a CoolProp AbstractState object.
    Use the **set_dimensions** method to specify the ejector geometry.
    Use the **set_efficiencies** method to specify the ejector efficiencies and other experimental parameters.

    What the ejector calculates:
    * Primary and secondary mass_flows (**m_p** and **m_s**)
    * Entrainment ratio (**entrainment_ratio**)
    * Critical backpressure (**Pd_crit**)
    * Area Ratio (**Ar**, which is A_3/A_t)
    * Pressure Ratio (**Pr**, which is P_s0/P_p0)
    
    notice that A_3 = A_2 = A_4 = A_const
    """
    cdef object _cp_update
    cdef object _cp_hmass
    cdef object _cp_smass
    cdef object _cp_rhomass
    cdef object _cp_speed_sound
    cdef object _cp_T
    cdef object _cp_Q
    cdef object _cp_cpmass
    cdef object _cp_isobaric_expansion_coefficient
    cdef object _cp_phase
    cdef object _cp_specify_phase
    cdef object _cp_unspecify_phase


    # =========================================================================
    # DECLARAÇÃO ESTÁTICA DE TIPOS
    # =========================================================================
    cdef public double Pt_tol, Pp1_tol, Pconst_tol, rho4_tol, P5_tol
    cdef public int max_iter_Pt, max_iter_Pp1, max_iter_Pconst, max_iter_rho4, max_iter_P5

    # ======CONVERGÊNCIA DE CADA LAÇO SECANTE ======
    cdef public bint converged_Pt, converged_Pp1, converged_Pconst, converged_rho4, converged_P5
    cdef public int iter_Pt, iter_Pp1, iter_Pconst, iter_rho4, iter_P5
    cdef public double residual_Pt, residual_Pp1, residual_Pconst, residual_rho4, residual_P5
    
    cdef public double P_p0, s_p0, h_p0, rho_p0
    cdef public double P_s0, s_s0, h_s0, rho_s0
    
    cdef public double eta_t, eta_m, eta_d, phi_m, psi
    cdef public double d_t, d_p1, d_const, A_t, A_p1, A_const, Ar
    
    cdef public double P_t, P_p1, P_const, rho_4, P_5
    
    cdef public double h_t_is, h_t, V_t, a_t, rho_t, s_t, m_p, deviation_t
    cdef public double s_p1, h_p1, rho_p1, V_p1, V_p1_line, deviation_p1
    cdef public double h_2h, V_2h, rho_2h, A_2h, A_p2, A_s2
    cdef public double h_p2_is, h_s2, rho_s2, h_p2, rho_p2, V_p2, V_s2, m_s, deviation_ms
    cdef public double P_3, P_2, V_3, m, h_3, rho_3, s_3
    cdef public double P_4, h_4, V_4, rho_4_line, s_4, deviation_rho_4
    cdef public double h_5, h_5_is
    
    cdef public double entrainment_ratio, Pd_crit, P_lift_ratio
    cdef public bint converged           
    cdef public double max_residual      
    
    cdef public object fluid
    cdef public MaterialStream primary_inlet_stream
    cdef public MaterialStream secondary_inlet_stream
    cdef public str dimension

    def __init__(self, MaterialStream primary_inlet_stream, MaterialStream secondary_inlet_stream, object fluid):
        
        #====== TOLERANCE OF EACH ITERATIVE LOOP AND MAX ITER ======
        self.Pt_tol = 1.5e-8
        self.Pp1_tol = 1.5e-8
        self.Pconst_tol = 1.5e-8 
        self.rho4_tol = 1.5e-8
        self.P5_tol = 1.5e-8

        self.max_iter_Pt = 50
        self.max_iter_Pp1 = 50
        self.max_iter_Pconst = 60
        self.max_iter_rho4 = 50
        self.max_iter_P5 = 50

        #==============================================

        self.fluid = fluid
        self.primary_inlet_stream = primary_inlet_stream
        self.secondary_inlet_stream = secondary_inlet_stream

        self.P_p0 = primary_inlet_stream.p
        self.s_p0 = primary_inlet_stream.s
        self.h_p0 = primary_inlet_stream.h
        self.rho_p0 = primary_inlet_stream.rho

        self.P_s0 = secondary_inlet_stream.p
        self.s_s0 = secondary_inlet_stream.s
        self.h_s0 = secondary_inlet_stream.h
        self.rho_s0 = secondary_inlet_stream.rho

        #======Initial experimental parameters=========
        self.eta_t = 0.95
        self.eta_m = 0.95
        self.eta_d = 0.95
        self.phi_m = 0.87
        self.psi = 0.88

        #Initial dimensions
        self.d_t = 2.64 / 1000.0
        self.d_p1 = 4.5 / 1000.0
        self.d_const = 6.7 / 1000.0
        self.A_t = PI * (self.d_t / 2.0)**2
        self.A_p1 = PI * (self.d_p1 / 2.0)**2
        self.A_const = PI * (self.d_const / 2.0)**2
        self.Ar = self.A_const / self.A_t

        self._cp_update = self.fluid.update
        self._cp_hmass = self.fluid.hmass
        self._cp_smass = self.fluid.smass
        self._cp_rhomass = self.fluid.rhomass
        self._cp_speed_sound = self.fluid.speed_sound
        self._cp_T = self.fluid.T
        self._cp_Q = self.fluid.Q
        self._cp_cpmass = self.fluid.cpmass
        self._cp_isobaric_expansion_coefficient = self.fluid.isobaric_expansion_coefficient
        self._cp_phase = self.fluid.phase
        self._cp_specify_phase = self.fluid.specify_phase
        self._cp_unspecify_phase = self.fluid.unspecify_phase

    cpdef void set_dimensions(self, str dimension, double value_t, double value_p1, double value_const):
        self.dimension = dimension
        if value_t <= 0 or value_p1 <= 0 or value_const <= 0:
            raise ValueError("Dimensions must be positive values.")
            
        if self.dimension.strip().lower()[0] == "d":
            self.d_t = value_t
            self.d_p1 = value_p1
            self.d_const = value_const
            self.A_t = PI * (self.d_t / 2.0)**2
            self.A_p1 = PI * (self.d_p1 / 2.0)**2
            self.A_const = PI * (self.d_const / 2.0)**2
            self.Ar = self.A_const / self.A_t
        elif self.dimension.strip().lower()[0] == "a":
            self.A_t = value_t
            self.A_p1 = value_p1
            self.A_const = value_const
            self.Ar = self.A_const / self.A_t
            self.d_t = sqrt(self.A_t * 4.0 / PI)
            self.d_p1 = sqrt(self.A_p1 * 4.0 / PI)
            self.d_const = sqrt(self.A_const * 4.0 / PI)
        else:
            raise ValueError("Invalid dimension type. Use 'd' for diameter or 'a' for area.")

    cpdef void set_efficiencies(self, double eta_t, double eta_m, double eta_d, double phi_m, double psi):
        if not (0 < eta_t <= 1 and 0 < eta_m <= 1 and 0 < eta_d <= 1 and 0 < phi_m <= 1 and 0 < psi <= 1):
            raise ValueError("Efficiencies must be between 0 and 1.")
        self.eta_t = eta_t
        self.eta_m = eta_m
        self.eta_d = eta_d
        self.phi_m = phi_m
        self.psi = psi

    # =========================================================================
    # SOLUCIONADORES ROOT-FINDING NATIVOS EM C (Substitutos de opt.fsolve)
    # =========================================================================
    #P.S: eu screvi uma secante para cada um porque ai eu não preciso passar a função, o que iria requerer uma passagem pelo python
    #(do jeito que tá dá pra ver no annotate que tá tudo branquinho, só o except * passa pelo python)
    # NOTA: critério de convergência é o RESÍDUO real fabs(f1).
    # fabs(f1-f0) só mede se o secante parou de se mexer entre iterações — o que
    # pode acontecer por estagnação (f0≈f1 sem nenhum dos dois estar perto de
    # zero), "convergindo" numericamente sem ter achado a raiz de verdade.
    # Ao final, converged_X/iter_X/residual_X ficam registrados como atributos
    # públicos.
    cdef double _secant_P_t(self, double guess) except *:
        cdef double P0 = guess
        cdef double P1 = guess * 0.99
        cdef double f0 = self.nozzle_p_throat_c(P0)
        cdef double f1 = self.nozzle_p_throat_c(P1)
        cdef int i = 0
        cdef double P_next
        while fabs(f1) > self.Pt_tol and i < self.max_iter_Pt:
            if fabs(f1 - f0) < 1e-14:  # denominador ~0: parou de se mexer
                break
            P_next = P1 - f1 * (P1 - P0) / (f1 - f0)
            P0 = P1; f0 = f1
            P1 = P_next
            f1 = self.nozzle_p_throat_c(P1)
            i += 1
        self.iter_Pt = i
        self.residual_Pt = fabs(f1)
        self.converged_Pt = self.residual_Pt <= self.Pt_tol
        return P1

    cdef double _secant_P_p1(self, double guess) except *:
        cdef double P0 = guess
        cdef double P1 = guess * 1.05
        cdef double f0 = self.nozzle_p_exit_c(P0)
        cdef double f1 = self.nozzle_p_exit_c(P1)
        cdef int i = 0
        cdef double P_next
        while fabs(f1) > self.Pp1_tol and i < self.max_iter_Pp1:
            if fabs(f1 - f0) < 1e-14:
                break
            P_next = P1 - f1 * (P1 - P0) / (f1 - f0)
            P0 = P1; f0 = f1
            P1 = P_next
            f1 = self.nozzle_p_exit_c(P1)
            i += 1
        self.iter_Pp1 = i
        self.residual_Pp1 = fabs(f1)
        self.converged_Pp1 = self.residual_Pp1 <= self.Pp1_tol
        return P1

    cdef double _secant_P_const(self, double guess) except *:
        cdef double P0 = guess
        cdef double P1 = guess * 0.95
        cdef double f0 = self.aerodynamic_throat_c(P0)
        cdef double f1 = self.aerodynamic_throat_c(P1)
        cdef int i = 0
        cdef double P_next
        while fabs(f1) > self.Pconst_tol and i < self.max_iter_Pconst:
            if fabs(f1 - f0) < 1e-14:
                break
            P_next = P1 - f1 * (P1 - P0) / (f1 - f0)
            P0 = P1; f0 = f1
            P1 = P_next
            f1 = self.aerodynamic_throat_c(P1)
            i += 1
        self.iter_Pconst = i
        self.residual_Pconst = fabs(f1)
        self.converged_Pconst = self.residual_Pconst <= self.Pconst_tol
        return P1

    cdef double _secant_rho_4(self, double guess) except *:
        cdef double x0 = guess
        cdef double x1 = guess * 1.05
        cdef double f0 = self.shock_c(x0)
        cdef double f1 = self.shock_c(x1)
        cdef int i = 0
        cdef double x_next
        while fabs(f1) > self.rho4_tol and i < self.max_iter_rho4:
            if fabs(f1 - f0) < 1e-14:
                break
            x_next = x1 - f1 * (x1 - x0) / (f1 - f0)
            x0 = x1; f0 = f1
            x1 = x_next
            f1 = self.shock_c(x1)
            i += 1
        self.iter_rho4 = i
        self.residual_rho4 = fabs(f1)
        self.converged_rho4 = self.residual_rho4 <= self.rho4_tol
        return x1

    cdef double _secant_P_5(self, double guess) except *:
        cdef double P0 = guess
        cdef double P1 = guess * 1.02
        cdef double f0 = self.diffuser_pressure_c(P0)
        cdef double f1 = self.diffuser_pressure_c(P1)
        cdef int i = 0
        cdef double P_next
        while fabs(f1) > self.P5_tol and i < self.max_iter_P5:
            if fabs(f1 - f0) < 1e-14:
                break
            P_next = P1 - f1 * (P1 - P0) / (f1 - f0)
            P0 = P1; f0 = f1
            P1 = P_next
            f1 = self.diffuser_pressure_c(P1)
            i += 1
        self.iter_P5 = i
        self.residual_P5 = fabs(f1)
        self.converged_P5 = self.residual_P5 <= self.P5_tol
        return P1

    # =========================================================================
    # LÓGICA FÍSICA 
    # O 'except *' avisa ao Cython para repassar erros pro Python
    # =========================================================================
    cdef double sound_velocity_c(self, double P, double h) except *:
        self._cp_update(HmassP_INPUTS, h, P)
        cdef double a, Q, T, rho, a_liq, rho_liq, cp_liq, beta_liq, v_liq, T_sat_v, T_sat_l
        cdef double a_vap, rho_vap, cp_vap, beta_vap, v_vap, e_vap, e_liq
        cdef double CP_liq, CP_vap, tau_liq, tau_vap, a_w

        if self._cp_phase() != iphase_twophase:
            a = self._cp_speed_sound()
        else:
            Q = self._cp_Q()
            T = self._cp_T()
            rho = self._cp_rhomass()

            
            self._cp_specify_phase(iphase_liquid)
            self._cp_update(PQ_INPUTS, P, 0.0)
            T_sat_l = self._cp_T()
            self._cp_update(PT_INPUTS, P, T_sat_l-1e-4) #pra não dar erro de ele pegar propriedades do lado errado do domo. Só é necessário com BICUBIC
            a_liq = self._cp_speed_sound()
            rho_liq = self._cp_rhomass()
            cp_liq = self._cp_cpmass()
            beta_liq = self._cp_isobaric_expansion_coefficient()
            v_liq = 1.0 / rho_liq

            self._cp_specify_phase(iphase_gas)
            self._cp_update(PQ_INPUTS, P, 1.0)
            T_sat_v = self._cp_T()
            self._cp_update(PT_INPUTS, P, T_sat_v+1e-4)
            a_vap = self._cp_speed_sound()
            rho_vap = self._cp_rhomass()
            cp_vap = self._cp_cpmass()
            beta_vap = self._cp_isobaric_expansion_coefficient()
            v_vap = 1.0 / rho_vap

            self._cp_unspecify_phase()

            e_vap = Q * rho_liq / (Q * rho_liq + (1.0 - Q) * rho_vap)
            e_liq = 1.0 - e_vap

            CP_liq = rho_liq * e_liq * cp_liq
            CP_vap = rho_vap * e_vap * cp_vap

            tau_liq = T * beta_liq * v_liq / cp_liq
            tau_vap = T * beta_vap * v_vap / cp_vap

            a_w = rho * (e_vap / (rho_vap * (a_vap**2)) + e_liq / (rho_liq * (a_liq**2)))
            a = (a_w + (rho / T) * (CP_vap * CP_liq * ((tau_liq - tau_vap)**2)) / (CP_vap + CP_liq))**(-0.5)
        return a

    cdef double nozzle_p_throat_c(self, double P_t) except *:
        self.P_t = P_t
        self._cp_update(PSmass_INPUTS, self.P_t, self.s_p0)
        self.h_t_is = self._cp_hmass()
        self.h_t = self.h_p0 - self.eta_t * (self.h_p0 - self.h_t_is) 
        self.V_t = sqrt(2.0 * (self.h_p0 - self.h_t))
        self.a_t = self.sound_velocity_c(self.P_t, self.h_t)
        self._cp_update(HmassP_INPUTS, self.h_t, self.P_t)
        self.rho_t = self._cp_rhomass()
        self.s_t = self._cp_smass()
        self.m_p = self.rho_t * self.V_t * self.A_t
        self.deviation_t = self.V_t - self.a_t
        return self.deviation_t
    
    cdef double nozzle_p_exit_c(self, double P_p1) except *:
        self.P_p1 = P_p1
        self.s_p1 = self.s_t
        self._cp_update(PSmass_INPUTS, self.P_p1, self.s_p1)
        self.h_p1 = self._cp_hmass()
        self.rho_p1 = self._cp_rhomass()
        self.V_p1 = self.m_p / (self.rho_p1 * self.A_p1)
        self.V_p1_line = sqrt(2.0 * (self.h_t - self.h_p1) + self.V_t**2)
        self.deviation_p1 = self.V_p1 - self.V_p1_line
        return self.deviation_p1
    
    cdef double aerodynamic_throat_c(self, double P_const) except *:
        self.P_const = P_const
        self._cp_update(PSmass_INPUTS, self.P_const, self.s_p1)
        self.h_p2_is = self._cp_hmass()
        self._cp_update(PSmass_INPUTS, self.P_const, self.s_s0)
        self.h_s2 = self._cp_hmass()
        self.rho_s2 = self._cp_rhomass()
        self.h_p2 = self.h_p1 - self.eta_m * (self.h_p1 - self.h_p2_is)
        self._cp_update(HmassP_INPUTS, self.h_p2, self.P_const)
        self.rho_p2 = self._cp_rhomass()
        self.V_p2 = self.m_p / (self.rho_p2 * self.A_p2)
        self.V_s2 = self.sound_velocity_c(self.P_const, self.h_s2)
        self.m_s = self.rho_s2 * self.V_s2 * self.A_s2
        
        self.deviation_ms = (self.m_p * (self.h_p1 + self.V_p1**2 / 2.0) + self.m_s * self.h_s0 - 
                             self.m_p * (self.h_p2 + self.V_p2**2 / 2.0) - 
                             self.m_s * (self.h_s2 + self.V_s2**2 / 2.0)) / \
                            (self.m_p * (self.h_p1 + self.V_p1**2 / 2.0) + self.m_s * self.h_s0)
        return self.deviation_ms

    cdef void mixing_chamber_c(self) except *:
        self.P_3 = self.P_const
        self.P_2 = self.P_const
        self.V_3 = self.phi_m * (self.m_p * self.V_p2 + self.m_s * self.V_s2) / (self.m_p + self.m_s)
        self.m = self.m_p + self.m_s
        self.h_3 = (self.m_p / self.m) * (self.h_p2 + self.V_p2**2 / 2.0) + \
                   (self.m_s / self.m) * (self.h_s2 + self.V_s2**2 / 2.0) - (self.V_3**2 / 2.0)
        self._cp_update(HmassP_INPUTS, self.h_3, self.P_3)
        self.rho_3 = self._cp_rhomass()
        self.s_3 = self._cp_smass()

    cdef double shock_c(self, double rho_4) except *:
        self.rho_4 = rho_4
        self.P_4 = self.V_3**2 * (self.rho_4 - self.rho_3) * self.rho_3 / self.rho_4 + self.P_3
        self.h_4 = (self.P_4 - self.P_3) / 2.0 * ((self.rho_3 + self.rho_4) / (self.rho_3 * self.rho_4)) + self.h_3
        self.V_4 = sqrt(2.0 * (self.h_3 - self.h_4) + self.V_3**2)
        self._cp_update(HmassP_INPUTS, self.h_4, self.P_4)
        self.rho_4_line = self._cp_rhomass()
        self.s_4 = self._cp_smass()
        self.deviation_rho_4 = self.rho_4_line - self.rho_4
        return self.deviation_rho_4

    cdef void diffuser_c(self) except *:
        self.h_5 = self.h_4 + self.V_4**2 / 2.0
        self.h_5_is = (self.h_5 - self.h_4) / self.eta_d + self.h_4
        cdef double P_5_guess = self.P_4 * 1.2
        self.P_5 = self._secant_P_5(P_5_guess)

    cdef double diffuser_pressure_c(self, double P_5) except *:
        self._cp_update(PSmass_INPUTS, P_5, self.s_4)
        cdef double h_guess = self._cp_hmass()
        return h_guess - self.h_5_is  

    # =========================================================================
    # MÉTODO PRINCIPAL (Acessível via Python)
    # =========================================================================
    cpdef void calculate(self, double P_t=0.0, double P_p1=0.0, double P_const=0.0, double rho_4=0.0) except *:
        self.P_p0 = self.primary_inlet_stream.p
        self.s_p0 = self.primary_inlet_stream.s
        self.h_p0 = self.primary_inlet_stream.h
        self.rho_p0 = self.primary_inlet_stream.rho

        self.P_s0 = self.secondary_inlet_stream.p
        self.s_s0 = self.secondary_inlet_stream.s
        self.h_s0 = self.secondary_inlet_stream.h
        self.rho_s0 = self.secondary_inlet_stream.rho

        cdef double P_t_guess = P_t if P_t != 0.0 else self.P_p0 * 0.58
        self.P_t = self._secant_P_t(P_t_guess)

        cdef double P_p1_guess = P_p1 if P_p1 != 0.0 else self.P_t * 0.1
        self.P_p1 = self._secant_P_p1(P_p1_guess)

        self._cp_update(PSmass_INPUTS, self.P_s0, self.s_p1)
        self.h_2h = self._cp_hmass()
        self.V_2h = sqrt(2.0 * (self.h_p1 - self.h_2h) + self.V_p1**2)
        self.rho_2h = self._cp_rhomass()
        self.A_2h = self.m_p / (self.rho_2h * self.V_2h)
        self.A_p2 = self.A_2h / (self.psi**2)
        self.A_s2 = self.A_const - self.A_p2

        cdef double P_const_guess = P_const if P_const != 0.0 else self.P_s0 * 0.65
        self.P_const = self._secant_P_const(P_const_guess)

        self.mixing_chamber_c()

        cdef double rho_4_guess = rho_4 if rho_4 != 0.0 else self.rho_3 * 3.5
        self.rho_4 = self._secant_rho_4(rho_4_guess)

        self.diffuser_c()

        self.validate()

        self.entrainment_ratio = self.m_s / self.m_p
        self.Pd_crit = self.P_5
        self.P_lift_ratio = self.Pd_crit / self.P_s0

        self.converged = (self.converged_Pt and self.converged_Pp1 and
                           self.converged_Pconst and self.converged_rho4 and self.converged_P5)
        self.max_residual = self.residual_Pt
        if self.residual_Pp1 > self.max_residual: self.max_residual = self.residual_Pp1
        if self.residual_Pconst > self.max_residual: self.max_residual = self.residual_Pconst
        if self.residual_rho4 > self.max_residual: self.max_residual = self.residual_rho4
        if self.residual_P5 > self.max_residual: self.max_residual = self.residual_P5

    cpdef void validate(self):
        conditions = [
            (self.P_s0 < self.P_p0, "P_s0 > P_p0"),
            (self.P_t < self.P_p0,  "P_t > P_p0"),
            (self.P_p1 < self.P_t,  "P_p1 > P_t"),
            (self.P_2 < self.P_s0,  "P_2 > P_s0"),
            (self.m_p > 0,          "m_p <= 0"),
            (self.m_s > 0,          "m_s <= 0"),
            (self.rho_4 > self.rho_3, "rho_4 <= rho_3")
        ]
        for cond, erro_msg in conditions:
            if not cond:
                raise ValueError(f'Invalid Result: {erro_msg}')

    # Wrappers para uso manual caso o usuário queira chamar do Python
    def nozzle_p_throat(self, P_t): return self.nozzle_p_throat_c(P_t)
    def nozzle_p_exit(self, P_p1): return self.nozzle_p_exit_c(P_p1)
    def aerodynamic_throat(self, P_const): return self.aerodynamic_throat_c(P_const)
    def shock(self, rho_4): return self.shock_c(rho_4)
    def speed_sound(self, P, h): return self.sound_velocity_c(P, h)

