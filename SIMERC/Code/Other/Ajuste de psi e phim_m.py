import pandas as pd
import numpy as np
import scipy.optimize as opt
from CoolProp.CoolProp import AbstractState
import CoolProp.CoolProp as CP
import gsim as gs

#Pega a planilha com os valores
df_total = pd.read_excel("Fluidos dados.xlsx", sheet_name="R134a") #troca pra "R1234ze(E)" depois

TAMANHO_CALIBRACAO = 29 #depende da quantidade dos dados, pro R1234ze(E) eu usei 
SEED = 42

# Faz a separação em conjunto de calibração e validação
df_calib = df_total.sample(n=TAMANHO_CALIBRACAO, random_state=SEED)
df_valid = df_total.drop(df_calib.index).reset_index(drop=True)
df_calib = df_calib.reset_index(drop=True)

fluido = AbstractState("BICUBIC&HEOS", "R134a")


# Função auxiliar para processar os blocos de dados no CoolProp
def pre_processar_dados(df_subset):
    N = len(df_subset)
    h_s0, s_s0 = np.zeros(N), np.zeros(N)
    h_p0, s_p0 = np.zeros(N), np.zeros(N)
    Ar, Pr = np.zeros(N), np.zeros(N)
    P_p0 = df_subset['Pp0'].values
    P_s0 = df_subset['Ps0'].values
    
    stream_p_temp = gs.MaterialStream("R134a")
    stream_s_temp = gs.MaterialStream("R134a")
    ejetor_temp = gs.Ejector(stream_p_temp, stream_s_temp, fluido)
    
    for i in range(N):
        # Evaporador
        fluido.update(CP.QT_INPUTS, 1, df_subset.loc[i, 'Ts0'])
        h_s0[i] = fluido.hmass()
        s_s0[i] = fluido.smass()

        # Gerador
        fluido.update(CP.QT_INPUTS, 1, df_subset.loc[i, 'Tp0'])
        h_p0[i] = fluido.hmass()
        s_p0[i] = fluido.smass()

        # Geometria (passando de mm para m)
        ejetor_temp.set_dimensions('d', df_subset.loc[i, 'dt']/1000, 
                                        df_subset.loc[i, 'dp1']/1000, 
                                        df_subset.loc[i, 'd3']/1000)
        Ar[i] = ejetor_temp.A_const / ejetor_temp.A_t
        Pr[i] = P_s0[i] / P_p0[i]
        
    return P_s0, h_s0, s_s0, P_p0, h_p0, s_p0, Ar, Pr, df_subset['w'].values, df_subset['Pd*'].values

# Pré-processa os dois grupos separadamente
print("Processando dados do CoolProp...")
P_s0_c, h_s0_c, s_s0_c, P_p0_c, h_p0_c, s_p0_c, Ar_c, Pr_c, omega_exp_c, Pc_exp_c = pre_processar_dados(df_calib)
P_s0_v, h_s0_v, s_s0_v, P_p0_v, h_p0_v, s_p0_v, Ar_v, Pr_v, omega_exp_v, Pc_exp_v = pre_processar_dados(df_valid)

# =================================================================
# OTIMIZAÇÃO (Apenas com os pontos de calibração)
# =================================================================
stream_p_main = gs.MaterialStream(fluido)
stream_s_main = gs.MaterialStream(fluido)
ejetor = gs.Ejector(stream_p_main, stream_s_main, fluido)

def f_obj(param):
    a, b, c, d = param
    erro_tot = 0.0
    for i in range(TAMANHO_CALIBRACAO):
        stream_s_main.p, stream_s_main.h, stream_s_main.s = P_s0_c[i], h_s0_c[i], s_s0_c[i]
        stream_p_main.p, stream_p_main.h, stream_p_main.s = P_p0_c[i], h_p0_c[i], s_p0_c[i]

        ejetor.set_dimensions('d', df_calib.loc[i, 'dt']/1000, df_calib.loc[i, 'dp1']/1000, df_calib.loc[i, 'd3']/1000)
        
        phi_m = a - b * Ar_c[i]
        psi = c / (Pr_c[i] * Ar_c[i]) + d
        ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)
        
        try:
            ejetor.calculate()
            erro_omega = ((omega_exp_c[i] - ejetor.entrainment_ratio) / omega_exp_c[i])**2
            erro_pc = ((Pc_exp_c[i] - (ejetor.Pd_crit)) / Pc_exp_c[i])**2
            erro_tot += (erro_omega + erro_pc)
        except:
            erro_tot += 10.0 # Penalidade caso o a conta falhe
            
    return erro_tot

def restricoes_ineq(param):
    a, b, c, d = param
    phi_m_vetor = a - b * Ar_c
    psi_vetor = c / (Pr_c * Ar_c) + d
    return np.concatenate((phi_m_vetor - 0.1, 0.99 - phi_m_vetor, psi_vetor - 0.1, 0.99 - psi_vetor))

chute = [0.9788, 0.0073, 0.046, 0.764]
restricoes = {'type': 'ineq', 'fun': restricoes_ineq}

print("Iniciando Ajuste dos Parâmetros...")
resultado = opt.minimize(f_obj, chute, method='SLSQP', constraints=restricoes, options={'maxiter': 300})

# =================================================================
# AVALIAÇÃO DOS ERROS (Calibração e Validação)
# =================================================================
if resultado.success:
    a_opt, b_opt, c_opt, d_opt = resultado.x
    print(f"\nAjuste concluído! Coeficientes Ótimos: \na={a_opt:.6f}, \nb={b_opt:.6f}, \nc={c_opt:.6f}, \nd={d_opt:.6f}")
    
    # -----------------------------------------------------------------
    # ERROS DO CONJUNTO DE CALIBRAÇÃO
    # -----------------------------------------------------------------
    print(f"\n--- AVALIANDO OS {TAMANHO_CALIBRACAO} PONTOS DE CALIBRAÇÃO ---")
    erros_omega_calib = []
    erros_pc_calib = []
    
    for i in range(TAMANHO_CALIBRACAO):
        stream_s_main.p, stream_s_main.h, stream_s_main.s = P_s0_c[i], h_s0_c[i], s_s0_c[i]
        stream_p_main.p, stream_p_main.h, stream_p_main.s = P_p0_c[i], h_p0_c[i], s_p0_c[i]

        ejetor.set_dimensions('d', df_calib.loc[i, 'dt']/1000, df_calib.loc[i, 'dp1']/1000, df_calib.loc[i, 'd3']/1000)
        
        phi_m = a_opt - b_opt * Ar_c[i]
        psi = c_opt / (Pr_c[i] * Ar_c[i]) + d_opt
        ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)
        
        ejetor.calculate()
        
        err_omega = abs(omega_exp_c[i] - ejetor.entrainment_ratio) / omega_exp_c[i] * 100
        err_pc = abs(Pc_exp_c[i] - (ejetor.Pd_crit)) / Pc_exp_c[i] * 100
        
        erros_omega_calib.append(err_omega)
        erros_pc_calib.append(err_pc)
        
        print(f"Calibração {i+1:02d}: Erro Omega = {err_omega:.2f}% | Erro Pc* = {err_pc:.2f}%")

    # -----------------------------------------------------------------
    # ERROS DO CONJUNTO DE VALIDAÇÃO (Pontos ocultos)
    # -----------------------------------------------------------------
    TAMANHO_VALIDACAO = len(df_valid)
    print(f"\n--- INICIANDO VALIDAÇÃO COM OS {TAMANHO_VALIDACAO} PONTOS OCULTOS ---")
    
    erros_omega_valid = []
    erros_pc_valid = []
    
    for i in range(TAMANHO_VALIDACAO):
        stream_s_main.p, stream_s_main.h, stream_s_main.s = P_s0_v[i], h_s0_v[i], s_s0_v[i]
        stream_p_main.p, stream_p_main.h, stream_p_main.s = P_p0_v[i], h_p0_v[i], s_p0_v[i]

        ejetor.set_dimensions('d', df_valid.loc[i, 'dt']/1000, df_valid.loc[i, 'dp1']/1000, df_valid.loc[i, 'd3']/1000)
        
        phi_m = a_opt - b_opt * Ar_v[i]
        psi = c_opt / (Pr_v[i] * Ar_v[i]) + d_opt
        ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)
        
        ejetor.calculate()
        
        err_omega = abs(omega_exp_v[i] - ejetor.entrainment_ratio) / omega_exp_v[i] * 100
        err_pc = abs(Pc_exp_v[i] - (ejetor.Pd_crit)) / Pc_exp_v[i] * 100
        
        erros_omega_valid.append(err_omega)
        erros_pc_valid.append(err_pc)
        
        print(f"Validação {i+1:02d} : Erro Omega = {err_omega:.2f}% | Erro Pc* = {err_pc:.2f}%")
        
    # -----------------------------------------------------------------
    # RESUMO FINAL
    # -----------------------------------------------------------------
    print("\n============================================================")
    print("                      RESUMO DOS ERROS                      ")
    print("============================================================")
    print(f"-> CALIBRAÇÃO ({TAMANHO_CALIBRACAO} pontos):")
    print(f"   Média Erro Omega: {np.mean(erros_omega_calib):.2f}%")
    print(f"   Média Erro Pc*:   {np.mean(erros_pc_calib):.2f}%")
    print(f"\n-> VALIDAÇÃO ({TAMANHO_VALIDACAO} pontos ocultos):")
    print(f"   Média Erro Omega: {np.mean(erros_omega_valid):.2f}%")
    print(f"   Média Erro Pc*:   {np.mean(erros_pc_valid):.2f}%")
    print("============================================================")
    print("\n============================================================")
    print("                   FORMA FINAL DE ψ e φₘ                     ")
    print("============================================================")
    print(f"\n   ψ = {c_opt}/(Pr * Ar) + {d_opt}")
    print(f"   φₘ = {a_opt} - {b_opt} * Ar\n")
    print(f"   Backend: {fluido.backend_name()}")
    print("============================================================\n")
else:
    print("O otimizador falhou:", resultado.message)