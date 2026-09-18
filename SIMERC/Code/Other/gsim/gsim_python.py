import CoolProp.CoolProp as CP
from CoolProp import AbstractState
import numpy as np
import scipy.optimize as opt

"""
This was the first version of gsim, created without Cython. 
The gsim.pyx file was created from this code. In my tests, the Cython 
version is about 2 to 3 times faster than this code. However, I decided 
to leave this one available as well because it's easier to understand 
how it works, especially the ejector.
"""

class MaterialStream:
    def __init__(self, fluid: str| AbstractState, bicubic: bool = False):
        
        self.fluid = fluid
        self.h = None
        self.s = None
        self.T = None
        self.p = None
        self.Q = None
        self.rho = None 
        self.mass_flow = 1

        if isinstance(self.fluid, str):
            try:
                if bicubic == True:
                    self.fluid = AbstractState("BICUBIC&HEOS", self.fluid)
                else:
                    self.fluid = AbstractState("HEOS", self.fluid)
            except:
                raise ValueError(f"Error: {self.fluid} is not a valid fluid name in CoolProp.")
        elif isinstance(self.fluid, AbstractState):
            pass
        else:
            raise ValueError("Fluid must be a string or an AbstractState object.") 
            
    def calculate(self, flash_type: int, property_1: float, property_2: float):
        """
        Calculate the thermodynamic state using two properties and the CoolProp constant

        Example:
        my_stream = MaterialStream("Water")\n
        my_stream.calculate(CP.PQ_INPUTS, 101325, 1)\n
        print(my_stream.T)\n
        #output: 373.15 (K)
        """
        self.fluid.update(flash_type, property_1, property_2)
        self.h = self.fluid.hmass()
        self.s = self.fluid.smass()
        self.T = self.fluid.T()
        self.p = self.fluid.p()
        self.Q = self.fluid.Q()
        self.rho = self.fluid.rhomass()

class Pump:
    def __init__(self, inlet_stream: MaterialStream, outlet_stream: MaterialStream, fluid):
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.h_1 = self.inlet_stream.h
        self.s_1 = self.inlet_stream.s
        self.mass_flow = self.inlet_stream.mass_flow 
        self.isentropic_efficiency = 0.75
    def set_isentropic_efficiency(self, isentropic_efficiency, percent=False):
        if percent == True or percent == '%':
            self.isentropic_efficiency = isentropic_efficiency/100
        else:
            self.isentropic_efficiency = isentropic_efficiency
    def calculate(self, outlet_pressure: float):
        self.h_1 = self.inlet_stream.h
        self.s_1 = self.inlet_stream.s
        self.mass_flow = self.inlet_stream.mass_flow 
        self.fluid.update(CP.PSmass_INPUTS, outlet_pressure, self.s_1)
        self.h_2_is = self.fluid.hmass()
        self.h_2 = self.h_1 + (self.h_2_is - self.h_1)/self.isentropic_efficiency
        self.fluid.update(CP.HmassP_INPUTS, self.h_2, outlet_pressure)
        self.s_2 = self.fluid.smass()
        self.work = self.mass_flow*(self.h_2 - self.h_1)

        self.outlet_stream.mass_flow = self.mass_flow
        self.outlet_stream.p = outlet_pressure
        self.outlet_stream.T = self.fluid.T()
        self.outlet_stream.h = self.h_2
        self.outlet_stream.s = self.s_2
        self.outlet_stream.Q = self.fluid.Q()
        self.outlet_stream.rho = self.fluid.rhomass()

class Valve:
    def __init__(self, inlet_stream:MaterialStream, outlet_stream:MaterialStream, fluid):
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.h_1 = self.inlet_stream.h
        self.P_1 = self.inlet_stream.p
    
    def calculate(self, outlet_pressure:float):
        self.h_1 = self.inlet_stream.h
        self.P_1 = self.inlet_stream.p

        self.P_2 = outlet_pressure
        self.fluid.update(CP.HmassP_INPUTS, self.h_1, self.P_2)

        self.outlet_stream.mass_flow = self.inlet_stream.mass_flow
        self.outlet_stream.p = outlet_pressure
        self.outlet_stream.T = self.fluid.T()
        self.outlet_stream.h = self.fluid.hmass()
        self.outlet_stream.s = self.fluid.smass()
        self.outlet_stream.Q = self.fluid.Q()
        self.outlet_stream.rho = self.fluid.rhomass()
        
class HeaterCooler:
    def __init__(self, inlet_stream:MaterialStream, outlet_stream:MaterialStream, fluid, efficiency:float = 1):
        self.efficiency = efficiency
        self.inlet_stream = inlet_stream
        self.outlet_stream = outlet_stream
        self.fluid = fluid
        self.mass_flow = inlet_stream.mass_flow
        self.h_1 = inlet_stream.h
        self.P_1 = inlet_stream.p
        self.efficiency = efficiency

    def set_efficiency(self, efficiency, percent=False):
        if percent == True or percent == '%':
            self.efficiency = efficiency/100
        else:
            self.efficiency = efficiency
    
    def calculate(self, flash_type:int, value_1: float, value_2: float, efficiency:float=1):
        """
        Calculate the HeatCooler with the given parameters.\n
        Example:\n
        A heater receives 2 kg/s of water at 25 °C (1 atm) and produces saturated steam. Considering an efficiency of 90%, what will be the energy (heat) required and what will be the outlet temperature?\n
        water = AbstractState("HEOS", "Water")\n
        inlet_stream = MaterialStream(water)\n
        inlet_stream.mass_flow = 2\n
        outlet_stream = MaterialStream(water)\n
        inlet_stream.calculate(CP.PT_INPUTS, 101325, 298.15)\n
        heater = HeaterCooler(inlet_stream, outlet_stream, water)\n
        heater.set_efficiency(90,'%') \n
        heater.calculate(CP.PQ_INPUTS, 101325, 1)\n
        print(f'{heater.heat_duty/1000} kW') #output: 5712.46 kW\n
        print(f'{outlet_stream.T-273.15} °C') #output: 99.97 °C\n
        \n
        Another example:
        A cooler receives 10 kg/s of saturated benzene vapor at 5 bar and must cool it to 25 °C. The outlet pressure is 1 bar (pressure drop = 4 bar). How much heat must be removed, considering that the refrigerator operates with an efficiency of 96%?\n
        benzene = AbstractState("HEOS", "Benzene")\n
        inlet_stream = MaterialStream(benzene)\n
        inlet_stream.mass_flow = 10\n
        outlet_stream = MaterialStream(benzene)\n
        inlet_stream.calculate(CP.PQ_INPUTS, 500_000, 1)\n
        heater = HeaterCooler(inlet_stream, outlet_stream, benzene)\n
        heater.set_efficiency(96,'%')\n
        heater.calculate(CP.PT_INPUTS, 100_000, 298.15)\n
        print(f'{heater.heat_duty/1000} kW') #output: 5932.85 kW\n
        """
        self.mass_flow = self.inlet_stream.mass_flow
        self.fluid.update(flash_type, value_1, value_2)
        self.heat_duty = self.mass_flow*(self.fluid.hmass() - self.inlet_stream.h)/self.efficiency
        self.outlet_stream.mass_flow = self.mass_flow
        self.outlet_stream.p = self.fluid.p()
        self.outlet_stream.T = self.fluid.T()
        self.outlet_stream.h = self.fluid.hmass()
        self.outlet_stream.s = self.fluid.smass()
        self.outlet_stream.Q = self.fluid.Q()
        self.outlet_stream.rho = self.fluid.rhomass()

class Ejector:
    def __init__(self, primary_inlet_stream, secondary_inlet_stream, fluid):
        """
        A CPM ejector model. It uses the modeling developed by Cardemil and Cole (2012).
        The user must provide the inlet streams (primary and secondary) and the fluid.
        The fluid must be a CoolProp AbstractState object.
        Use the **set_dimensions** method to specify the ejector geometry.
        Use the **set_efficiencies** method to specify the ejector efficiencies and other experimental parameters.\n
        What the ejector calculates:\n
        * Primary and secondary mass_flows (**m_p** and **m_s**)
        * Entrainment ratio (**entrainmente_ratio**)
        * Critical backpressure (**Pd_crit**)
        * Area Ratio (**Ar**, wich is A_3/A_t)
        * Pressure Ratio (**Pr**, wich is P_s0/P_p0)\n
        notice that A_3 = A_2 = A_4 = A_const\n

        Nomenclature used for the ejector sections:\n
        ![Ejector Sections Nomenclature](https://github.com/GPCTM-BR/GSIM/blob/main/docstrings/images/ejector_docs.webp?raw=true)
        
        Reference:
        CARDEMIL, J. M.; COLLE, S. A general model for evaluation of vapor ejectors performance for application in refrigeration. Energy Conservation and Management, v. 64, p. 79-86, set. 2012. DOI: https://doi.org/10.1016/j.enconman.2012.05.009. Available at: https://www.sciencedirect.com/science/article/abs/pii/S0196890412002154. 
        """
        #====== TOLERANCE OF EACH ITERATIVE LOOP ======
        self.Pt_tol = 1.5e-8
        self.Pp1_tol = 1.5e-8
        self.Pconst_tol = 1.5e-8 #<---sometimes I need to set this one at 1.5e-3 to converge
        self.rho4_tol = 1.5e-8
        self.P5_tol = 1.5e-8
        #==============================================

        #========INLET STREAMS PROPERTIES==============
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
        #==============================================

        #======Initial experimental parameters=========
        self.eta_t = 0.95
        self.eta_m = 0.95
        self.eta_d = 0.95
        self.phi_m = 0.87
        self.psi = 0.88

        #Initial dimensions
        self.d_t = 2.64/1000
        self.d_p1 = 4.5/1000
        self.d_const = 6.7/1000
        self.A_t = np.pi * (self.d_t / 2)**2
        self.A_p1 = np.pi * (self.d_p1 / 2)**2
        self.A_const = np.pi * (self.d_const / 2)**2
        self.Ar = self.A_const/self.A_t
        #==============================================
    
    def set_dimensions(self, dimension: str, value_t: float, value_p1: float, value_const: float):
        """
        Set the dimensions of the ejector, where dimension can be 'd' for diameter or 'a' for area.\n
        value_t: diameter or area of the primary nozzle throat\n
        value_p1: diameter or area of the exit of the primary nozzle\n
        value_const: diameter or area of the constant area section\n

        See the image below for reference\n
        ![Ejector Dimensions](https://github.com/GPCTM-BR/GSIM/blob/main/docstrings/images/set_dimensions_docs.webp?raw=true)
        """
        self.dimension = dimension
        if not all(isinstance(value, (int, float)) for value in [value_t, value_p1, value_const]):
            raise ValueError("Dimensions must be numeric values.")
        if any(value <= 0 for value in [value_t, value_p1,value_const]):
            raise ValueError("Dimensions must be positive values.")
        if self.dimension.strip().lower()[0] == "d":
            self.d_t = value_t
            self.d_p1 = value_p1
            self.d_const = value_const
            self.A_t = np.pi * (self.d_t / 2)**2
            self.A_p1 = np.pi * (self.d_p1 / 2)**2
            self.A_const = np.pi * (self.d_const / 2)**2
            self.Ar = self.A_const/self.A_t
        elif self.dimension.strip().lower()[0] == "a":
            self.A_t = value_t
            self.A_p1 = value_p1
            self.A_const = value_const
            self.Ar = self.A_const/self.A_t
            self.d_t = np.sqrt(self.A_t * 4 / np.pi)
            self.d_p1 = np.sqrt(self.A_p1 * 4 / np.pi)
            self.d_const = np.sqrt(self.A_const * 4 / np.pi)
        else:
            raise ValueError("Invalid dimension type. Use 'd' for diameter or 'a' for area.")

    def set_efficiencies(self, eta_t: float, eta_m: float, eta_d: float, phi_m: float, psi: float):
        """
        Set the efficiencies of the ejector, where:\n
        eta_t (η_t): primary nozzle efficiency\n
        eta_m (η_m): mixing chamber efficiency\n
        eta_d (η_d): diffuser efficiency\n
        phi_m (φ_m): mixing loss factor\n
        psi (ψ): coefficient that accounts for the reduction of the aerodynamic throat area\n
        """
        if not all(isinstance(value, (int, float)) for value in [eta_t, eta_m, eta_d, phi_m, psi]):
            raise ValueError("Efficiencies must be numeric values.")
        if not all(0 < value <= 1 for value in [eta_t, eta_m, eta_d, phi_m, psi]):
            raise ValueError("Efficiencies must be between 0 and 1.")
        self.eta_t = eta_t
        self.eta_m = eta_m
        self.eta_d = eta_d
        self.phi_m = phi_m
        self.psi = psi

    def calculate(self, P_t: float = None, P_p1: float = None, P_const: float = None, rho_4: float = None):
        """
        Resolve the ejector model using the input parameters and dimensions. The user can choose to provide one or more of the following initial guesses (but is not necessary):\n
        P_t: Pressure at the throat (Pa)\n
        P_p1: Pressure at the exit of the primary nozzle (Pa)\n
        P_const: Pressure at the constant area section (Pa)\n
        rho_4: Density after the shock (kg/m³)\n
        """
        self.P_p0 = self.primary_inlet_stream.p
        self.s_p0 = self.primary_inlet_stream.s
        self.h_p0 = self.primary_inlet_stream.h
        self.rho_p0 = self.primary_inlet_stream.rho

        self.P_s0 = self.secondary_inlet_stream.p
        self.s_s0 = self.secondary_inlet_stream.s
        self.h_s0 = self.secondary_inlet_stream.h
        self.rho_s0 = self.secondary_inlet_stream.rho

        self.P_t = P_t
        self.P_p1 = P_p1
        self.P_const = P_const
        self.rho_4 = rho_4

        if self.P_t is None:
            self.P_t = self.P_p0*0.58
        result_t, infodict, ier, message = opt.fsolve(self.nozzle_p_throat, self.P_t, full_output=True, xtol= self.Pt_tol)
        if ier !=1:
            raise ValueError(f'Error in nozzle_p_throat: {message}')
        self.P_t = result_t[0]

        if self.P_p1 is None:
            self.P_p1 = self.P_t*0.1
        result_p1, infodict, ier, message = opt.fsolve(self.nozzle_p_exit, self.P_p1, full_output=True, xtol= self.Pp1_tol)
        if ier !=1:
            raise ValueError(f'Error in nozzle_p_exit: {message}')
        self.P_p1 = result_p1[0]

        self.fluid.update(CP.PSmass_INPUTS, self.P_s0, self.s_p1)
        self.h_2h = self.fluid.hmass()
        self.V_2h = np.sqrt(2*(self.h_p1 - self.h_2h) + self.V_p1**2)
        self.rho_2h = self.fluid.rhomass()
        self.A_2h = self.m_p/(self.rho_2h*self.V_2h)
        self.A_p2 = self.A_2h/(self.psi**2)
        self.A_s2 = self.A_const - self.A_p2

        if self.P_const is None:
            self.P_const = self.P_s0*0.65
        result_at, infodict, ier, message = opt.fsolve(self.aerodynamic_throat, self.P_const, full_output=True, xtol= self.Pconst_tol)
        if ier !=1:
            raise ValueError(f'Error in aerodynamic_throat: {message} Residue: {infodict['fvec']:.4e}')
        self.P_const = result_at[0]

        self.mixing_chamber()

        if self.rho_4 is None:
            self.rho_4 = self.rho_3*3.5
        result_shock, infodict, ier, message = opt.fsolve(self.shock, self.rho_4, full_output=True, xtol= self.rho4_tol)
        if ier !=1:
            raise ValueError(f'Error in shock: {message}')
        self.rho_4 = result_shock[0]

        self.diffuser()

        self.validate()

        #======== SALVANDO OS RESULTADOS ================
        self.entrainment_ratio = self.m_s/self.m_p
        self.Pd_crit = self.P_5
        self.P_lift_ratio = self.Pd_crit/self.P_s0

    def nozzle_p_throat(self, P_t):
        self.P_t = P_t[0] if isinstance(P_t, np.ndarray) else P_t
        self.fluid.update(CP.PSmass_INPUTS, self.P_t, self.s_p0)
        self.h_t_is = self.fluid.hmass()
        self.h_t = self.h_p0 - self.eta_t*(self.h_p0 - self.h_t_is) 
        self.V_t = np.sqrt(2*(self.h_p0 - self.h_t))
        self.a_t = self.sound_velocity(self.P_t, self.h_t, self.fluid)
        self.fluid.update(CP.HmassP_INPUTS, self.h_t, self.P_t)
        self.rho_t = self.fluid.rhomass()
        self.s_t = self.fluid.smass()
        self.m_p = self.rho_t*self.V_t*self.A_t
        self.deviation_t = self.V_t - self.a_t
        return self.deviation_t
    
    def nozzle_p_exit(self, P_p1):
        self.P_p1 = P_p1[0] if isinstance(P_p1, np.ndarray) else P_p1
        self.s_p1 = self.s_t
        self.fluid.update(CP.PSmass_INPUTS, self.P_p1, self.s_p1)
        self.h_p1 = self.fluid.hmass()
        self.rho_p1 = self.fluid.rhomass()
        self.V_p1 = self.m_p/(self.rho_p1*self.A_p1)
        self.V_p1_line = np.sqrt(2*(self.h_t - self.h_p1) + self.V_t**2)
        self.deviation_p1 = self.V_p1 - self.V_p1_line
        return self.deviation_p1
    
    def aerodynamic_throat(self, P_const):
        self.P_const = P_const[0] if isinstance(P_const, np.ndarray) else P_const
        self.fluid.update(CP.PSmass_INPUTS, self.P_const, self.s_p1)
        self.h_p2_is = self.fluid.hmass()
        self.fluid.update(CP.PSmass_INPUTS, self.P_const, self.s_s0)
        self.h_s2 = self.fluid.hmass()
        self.rho_s2 = self.fluid.rhomass()
        self.h_p2 = self.h_p1 - self.eta_m*(self.h_p1 - self.h_p2_is)
        self.fluid.update(CP.HmassP_INPUTS, self.h_p2, self.P_const)
        self.rho_p2 = self.fluid.rhomass()
        self.V_p2 = self.m_p/(self.rho_p2*self.A_p2)
        self.V_s2 = self.sound_velocity(self.P_const, self.h_s2, self.fluid)
        self.m_s = self.rho_s2*self.V_s2*self.A_s2
        #self.m_s_line = (self.m_p*(self.h_p2 + self.V_p2**2/2 - self.h_p1 - self.V_p1**2/2))/(self.h_s0 - self.h_s2 - self.V_s2**2/2)
        #self.deviation_ms = (self.m_s - self.m_s_line)/self.m_s
        self.deviation_ms = (self.m_p*(self.h_p1 + self.V_p1**2/2) + self.m_s*self.h_s0 - self.m_p*(self.h_p2 + self.V_p2**2/2) - self.m_s*(self.h_s2 + self.V_s2**2/2))/(self.m_p*(self.h_p1 + self.V_p1**2/2) + self.m_s*self.h_s0)
        return self.deviation_ms

    def mixing_chamber(self):
        self.P_3 = self.P_const
        self.P_2 = self.P_const
        self.V_3 = self.phi_m*(self.m_p * self.V_p2 + self.m_s * self.V_s2)/(self.m_p + self.m_s)
        self.m = self.m_p + self.m_s
        self.h_3 = self.m_p/self.m*(self.h_p2 + self.V_p2**2/2) + self.m_s/self.m*(self.h_s2 + self.V_s2**2/2) - self.V_3**2/2
        self.fluid.update(CP.HmassP_INPUTS, self.h_3, self.P_3)
        self.rho_3 = self.fluid.rhomass()
        self.s_3 = self.fluid.smass()

    def shock(self, rho_4):
        self.rho_4 = rho_4[0] if isinstance(rho_4, np.ndarray) else rho_4
        self.P_4 = self.V_3**2*(self.rho_4 - self.rho_3)*self.rho_3/self.rho_4+self.P_3
        self.h_4 = (self.P_4 - self.P_3)/2*((self.rho_3 + self.rho_4)/(self.rho_3 * self.rho_4)) + self.h_3
        self.V_4 = np.sqrt(2*(self.h_3 - self.h_4) + self.V_3**2)
        self.fluid.update(CP.HmassP_INPUTS, self.h_4, self.P_4)
        self.rho_4_line = self.fluid.rhomass()
        self.s_4 = self.fluid.smass()
        self.deviation_rho_4 = self.rho_4_line - self.rho_4
        return self.deviation_rho_4

    def diffuser(self):
        self.h_5 = self.h_4 + self.V_4**2/2
        self.h_5_is = (self.h_5 - self.h_4)/self.eta_d + self.h_4
        #self.fluid.update(CP.HmassSmass_INPUTS, self.h_5_is, self.s_4)
        P_5_guess = self.P_4*1.2
        self.P_5 = opt.fsolve(self.diffuser_pressure, P_5_guess, xtol= self.P5_tol)[0]
        #self.P_5 = self.fluid.p()

    def diffuser_pressure(self, P_5):
        P_5 = P_5[0] if isinstance(P_5, np.ndarray) else P_5
        self.fluid.update(CP.PSmass_INPUTS, P_5, self.s_4)
        h_guess = self.fluid.hmass()
        return h_guess - self.h_5_is  

    def sound_velocity(self, P, h, fluid):
        "TÁ ERRADO, OLHA O 'Notas.txt' !"
        fluid.update(CP.HmassP_INPUTS, h, P)
        try:
            a = fluid.speed_sound()
        except:
            Q = fluid.Q()
            T = fluid.T()
            rho = fluid.rhomass()

            fluid.update(CP.PQ_INPUTS, P, 0)
            a_liq = fluid.speed_sound()
            rho_liq = fluid.rhomass()
            cp_liq = fluid.cpmass()
            beta_liq = fluid.isobaric_expansion_coefficient()
            v_liq = 1 / rho_liq

            fluid.update(CP.PQ_INPUTS, P, 1)
            a_vap = fluid.speed_sound()
            rho_vap = fluid.rhomass()
            cp_vap = fluid.cpmass()
            beta_vap = fluid.isobaric_expansion_coefficient()
            v_vap = 1 / rho_vap

            e_vap = Q*rho_liq/(Q*rho_liq + (1-Q)*rho_vap)
            e_liq = 1 - e_vap

            CP_liq = rho_liq*e_liq*cp_liq
            CP_vap = rho_vap*e_vap*cp_vap

            tau_liq = T*beta_liq*v_liq/cp_liq
            tau_vap = T*beta_vap*v_vap/cp_vap

            a_w = rho*(e_vap/(rho_vap*(a_vap**2)) + e_liq/(rho_liq*(a_liq**2)))
            a = (a_w + (rho/T)*(CP_vap*CP_liq*((tau_liq-tau_vap)**2))/(CP_vap+CP_liq))**(-0.5)
        return a

    def validate(self):
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

