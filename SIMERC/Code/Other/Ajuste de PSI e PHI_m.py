import pandas as pd
import numpy as np
import scipy.optimize as opt
import CoolProp.CoolProp as CP
from CoolProp import AbstractState
import matplotlib.pyplot as plt
import gsim

fluido = "R134a"
backend = "BICUBIC&HEOS"

df = pd.read_excel("Fluidos dados.xlsx", sheet_name = fluido)
NUMERO_DE_PONTOS = len(df['Tp0'])
FRACAO_CALIBRACAO = 1# 0.7
SEED = 7 #numero inteiro usado para gerar a semente do random_state, garantindo que a amostra seja sempre a mesma (para que voce, pessoa que esta lendo, consiga reproduzir os resultados)

caso = 1 #caso 1: ajusta a equação pros dados experimentais
         #caso 2: acha o valor de psi e phi que minimizam o erro entre os dados experimentais e os calculados (para analizar como se comportam frente a Pr e Ar)


NOMES_PARAM = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']  # nomes dos parâmetros, na ordem usada por param/chute

def phi_m_func(Ar, Pr, param):
    a, b, c, d, e, f, g, h, i, j = param
    #return a + b * Ar  #forma antiga
    return a + b * Ar + c * Pr + d * (Ar*Pr) #+ e * (Ar)**2  #forma nova
    

def psi_func(Ar, Pr, param):
    a, b, c, d, e, f, g, h, i, j = param
    #return c / (Pr * Ar) + d  #forma antiga
    return f + g * Ar + h * Pr + i * (Ar*Pr) #+ j * (Ar)**2  #forma nova

def formula_str(param):
    a, b, c, d, e, f, g, h, i, j = param
    txt_psi = f"ψ =  {f} + {g} * Ar + {h} * Pr + {i} * (Ar * Pr)"
    txt_phi = f"φₘ = {a} + {b} * Ar + {c} * Pr + {d} * (Ar * Pr)"
    return txt_psi, txt_phi



def pre_processar_dados(df_subset):
    #essa função cria vetores com os valores de P, h, s, Ar e Pr para cada ponto da planilha, usando o CoolProp
    #isso vai ser usado mais tarde como entrada para a função de ajuste de PSI e PHI
    N = len(df_subset)
    h_s0, s_s0 = np.zeros(N), np.zeros(N)
    h_p0, s_p0 = np.zeros(N), np.zeros(N)
    Ar, Pr = np.zeros(N), np.zeros(N)
    P_p0 = df_subset['Pp0'].values
    P_s0 = df_subset['Ps0'].values
 
    stream_p_temp = gsim.MaterialStream(fluido)
    stream_s_temp = gsim.MaterialStream(fluido)
    ejetor_temp = gsim.Ejector(stream_p_temp, stream_s_temp, fluido)

    #nem todos os dados estão na saturação, então nao dá pra fazr a avaliação do estado como f(T, Q=1)
    #Porem, fazer a avaliação f(T,P) para dados que ESTÃO na saturação costuma a dar erro. Essa função corrige isso
    TOL_SATURACAO = 1e-3  
    def calcular_h_s(T, P):
        fluido.update(CP.QT_INPUTS, 1, T)
        P_sat = fluido.p()
        if abs(P - P_sat) / P_sat < TOL_SATURACAO:
            return fluido.hmass(), fluido.smass()   
        else:
            fluido.update(CP.PT_INPUTS, P, T)
            return fluido.hmass(), fluido.smass()
 
    for i in range(N):
        # Evaporador
        h_s0[i], s_s0[i] = calcular_h_s(df_subset.loc[i, 'Ts0'], P_s0[i])
 
        # Gerador
        h_p0[i], s_p0[i] = calcular_h_s(df_subset.loc[i, 'Tp0'], P_p0[i])
 
        # Geometria (passando de mm para m)
        ejetor_temp.set_dimensions('d', df_subset.loc[i, 'dt']/1000,
                                        df_subset.loc[i, 'dp1']/1000,
                                        df_subset.loc[i, 'd3']/1000)
        Ar[i] = ejetor_temp.A_const / ejetor_temp.A_t
        Pr[i] = P_s0[i] / P_p0[i]
 
    return P_s0, h_s0, s_s0, P_p0, h_p0, s_p0, Ar, Pr, df_subset['w'].values, df_subset['Pd*'].values

if caso == 1:
    TAMANHO_CALIBRACAO = round(NUMERO_DE_PONTOS*FRACAO_CALIBRACAO) 
    TAMANHO_VALIDACAO = NUMERO_DE_PONTOS - TAMANHO_CALIBRACAO 
    

    df_calib = df.sample(n=TAMANHO_CALIBRACAO, random_state=SEED)
    df_valid = df.drop(df_calib.index).reset_index(drop=True)
    df_calib = df_calib.reset_index(drop=True)

    fluido = AbstractState(backend, fluido)

    if df_valid.empty:
        print("Não há dados de validação. Todos os dados serão usados para calibração.")
        P_s0_c, h_s0_c, s_s0_c, P_p0_c, h_p0_c, s_p0_c, Ar_c, Pr_c, w_exp_c, Pd_exp_c = pre_processar_dados(df_calib)
    else:
        P_s0_c, h_s0_c, s_s0_c, P_p0_c, h_p0_c, s_p0_c, Ar_c, Pr_c, w_exp_c, Pd_exp_c = pre_processar_dados(df_calib)
        P_s0_v, h_s0_v, s_s0_v, P_p0_v, h_p0_v, s_p0_v, Ar_v, Pr_v, w_exp_v, Pd_exp_v = pre_processar_dados(df_valid)

    corrente_principal = gsim.MaterialStream(fluido)
    corrente_secundaria = gsim.MaterialStream(fluido)
    ejetor = gsim.Ejector(corrente_principal, corrente_secundaria, fluido)

    def funcao_objetivo(param):
        erro_total = 0.0
        for i in range(TAMANHO_CALIBRACAO):
            corrente_secundaria.p, corrente_secundaria.h, corrente_secundaria.s = P_s0_c[i], h_s0_c[i], s_s0_c[i]
            corrente_principal.p, corrente_principal.h, corrente_principal.s = P_p0_c[i], h_p0_c[i], s_p0_c[i]
            ejetor.set_dimensions('d', df_calib.loc[i, 'dt']/1000, 
                                        df_calib.loc[i, 'dp1']/1000,
                                        df_calib.loc[i, 'd3']/1000)

            phi_m = phi_m_func(Ar_c[i], Pr_c[i], param)
            psi = psi_func(Ar_c[i], Pr_c[i], param)
            

            try:
                ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)
                ejetor.calculate()
                erro_w = ((ejetor.entrainment_ratio - w_exp_c[i])/ w_exp_c[i])**2
                erro_Pd = ((ejetor.Pd_crit - Pd_exp_c[i])/ Pd_exp_c[i])**2
                erro_total += (erro_w + erro_Pd)
            except:
                erro_total += 10  # Penalidade caso a conta falhe

        return erro_total

    def funcao_restricoes(param):
        phi_m = phi_m_func(Ar_c, Pr_c, param)
        psi = psi_func(Ar_c, Pr_c, param)
        return np.concatenate((phi_m - 0.1, 0.99 - phi_m, psi - 0.1, 0.99 - psi))  # Restrições: 0.1 < phi_m < 0.99 e 0.1 < psi < 0.99

    chute = [0.9788, 0.0073, 0.046, 0.75, 0.75, 0.5, 0.25, 0.25, 0.25, 0.25]  
    #chute = [5, -6, 0.7, -0.5, 0.05, -7.5, 4, -0.5, -5.5, 1] 
    #chute = [0.8, 0.01, 1, -0.15, 1, 2, -0.2, -5.5, 1, 1]
    #chute = [0.8, 0.01, 0.01, -0.15, 1, 2, -0.2, -5.5, 1, 1]
    restricoes = {'type': 'ineq', 'fun': funcao_restricoes}

    #INICIO DO AJUSTE

    print("Iniciando Ajuste dos Parâmetros...")
    #o método SLSQP é usado para otimização com restrições, SLSQP significa Sequential Least Squares Programming
    resultado = opt.minimize(funcao_objetivo, chute, method='SLSQP', constraints=restricoes, options={'maxiter': 300})

    if resultado.success:
        param_opt = resultado.x
        print("\nAjuste concluído! Coeficientes Ótimos:")
        for nome, valor in zip(NOMES_PARAM, param_opt):
            print(f"{nome}={valor:.6f}")

        #Erros para os dados de calibração
        print("\nCalculando erros para os dados de calibração...")
        erro_w_calib = np.zeros(TAMANHO_CALIBRACAO)
        erro_Pd_calib = np.zeros(TAMANHO_CALIBRACAO)
        w_calc_calib = np.zeros(TAMANHO_CALIBRACAO)
        Pd_calc_calib = np.zeros(TAMANHO_CALIBRACAO)
        for i in range(TAMANHO_CALIBRACAO):
            corrente_secundaria.p, corrente_secundaria.h, corrente_secundaria.s = P_s0_c[i], h_s0_c[i], s_s0_c[i]
            corrente_principal.p, corrente_principal.h, corrente_principal.s = P_p0_c[i], h_p0_c[i], s_p0_c[i]
            ejetor.set_dimensions('d', df_calib.loc[i, 'dt']/1000, 
                                        df_calib.loc[i, 'dp1']/1000,
                                        df_calib.loc[i, 'd3']/1000)

            phi_m = phi_m_func(Ar_c[i], Pr_c[i], param_opt)
            psi = psi_func(Ar_c[i], Pr_c[i], param_opt)
            ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)

            ejetor.calculate()

            w_calc_calib[i] = ejetor.entrainment_ratio
            Pd_calc_calib[i] = ejetor.Pd_crit

            erro_w_calib[i] = (abs(ejetor.entrainment_ratio - w_exp_c[i])/ w_exp_c[i]) * 100
            erro_Pd_calib[i] = (abs(ejetor.Pd_crit - Pd_exp_c[i])/ Pd_exp_c[i]) * 100

            #print(f"Calibração {i+1:02d}: Erro Omega = {erro_w_calib[i]:.2f}% | Erro Pc* = {erro_Pd_calib[i]:.2f}%")
            print(f"Tp0: {df_calib.loc[i, 'Tp0']:.2f} | Pp0: {df_calib.loc[i, 'Pp0']:.2f} | Ts0: {df_calib.loc[i, 'Ts0']:.2f} | Ps0: {df_calib.loc[i, 'Ps0']:.4f} | Ar: {Ar_c[i]:.4f} | Pr: {Pr_c[i]:.4f} | phi_m: {phi_m:.4f} | psi: {psi:.4f} | w_calc: {w_calc_calib[i]:.4f} | Pd_calc: {Pd_calc_calib[i]:.4f}")


        #Erros para os dados de validação, se houver
        if not df_valid.empty:
            print("\nCalculando erros para os dados de validação...")
            erro_w_valid = np.zeros(TAMANHO_VALIDACAO)
            erro_Pd_valid = np.zeros(TAMANHO_VALIDACAO)
            w_calc_valid = np.zeros(TAMANHO_VALIDACAO)
            Pd_calc_valid = np.zeros(TAMANHO_VALIDACAO)
            for i in range(TAMANHO_VALIDACAO):
                corrente_secundaria.p, corrente_secundaria.h, corrente_secundaria.s = P_s0_v[i], h_s0_v[i], s_s0_v[i]
                corrente_principal.p, corrente_principal.h, corrente_principal.s = P_p0_v[i], h_p0_v[i], s_p0_v[i]
                ejetor.set_dimensions('d', df_valid.loc[i, 'dt']/1000, 
                                            df_valid.loc[i, 'dp1']/1000,
                                            df_valid.loc[i, 'd3']/1000)

                phi_m = phi_m_func(Ar_v[i], Pr_v[i], param_opt)
                psi = psi_func(Ar_v[i], Pr_v[i], param_opt)
                ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)

                ejetor.calculate()

                w_calc_valid[i] = ejetor.entrainment_ratio
                Pd_calc_valid[i] = ejetor.Pd_crit

                erro_w_valid[i] = (abs(ejetor.entrainment_ratio - w_exp_v[i])/ w_exp_v[i]) * 100
                erro_Pd_valid[i] = (abs(ejetor.Pd_crit - Pd_exp_v[i])/ Pd_exp_v[i]) * 100

                print(f"Validação {i+1:02d}: Erro Omega = {erro_w_valid[i]:.2f}% | Erro Pc* = {erro_Pd_valid[i]:.2f}%")


        print("\n============================================================")
        print("                      RESUMO DOS ERROS                      ")
        print("============================================================")
        print(f"-> CALIBRAÇÃO ({TAMANHO_CALIBRACAO} pontos):")
        print("    Para a Razão de Arraste:")
        print(f"   \tErro Mínimo: {np.min(erro_w_calib):.2f}%")
        print(f"   \tErro Máximo: {np.max(erro_w_calib):.2f}%")
        print(f"   \tErro Médio: {np.mean(erro_w_calib):.2f}%")
        print("    Para a Pressão Crítica:")
        print(f"   \tErro Mínimo: {np.min(erro_Pd_calib):.2f}%")
        print(f"   \tErro Máximo: {np.max(erro_Pd_calib):.2f}%")
        print(f"   \tErro Médio: {np.mean(erro_Pd_calib):.2f}%")
        if not df_valid.empty:
            print(f"\n-> VALIDAÇÃO ({TAMANHO_VALIDACAO} pontos ocultos):")
            print("    Para a Razão de Arraste:")
            print(f"   \tErro Mínimo: {np.min(erro_w_valid):.2f}%")
            print(f"   \tErro Máximo: {np.max(erro_w_valid):.2f}%")
            print(f"   \tErro Médio: {np.mean(erro_w_valid):.2f}%")
            print("    Para a Pressão Crítica:")
            print(f"   \tErro Mínimo: {np.min(erro_Pd_valid):.2f}%")
            print(f"   \tErro Máximo: {np.max(erro_Pd_valid):.2f}%")
            print(f"   \tErro Médio: {np.mean(erro_Pd_valid):.2f}%")
        
        print("============================================================")
        print("\n============================================================")
        print("                   FORMA FINAL DE ψ e φₘ                     ")
        print("============================================================")
        txt_psi, txt_phi = formula_str(param_opt)
        print(f"\n   {txt_psi}")
        print(f"   {txt_phi}\n")
        print(f"   Backend: {fluido.backend_name()}")
        print("============================================================\n")

    else:
        print("O ajuste não foi bem-sucedido. Tente outro chute inicial ou verifique os dados.")

elif caso == 2:
    #esse caso dois é pra achar "psi e phi_m experimental". Achar o psi e phi perfeitos que levam ao dado experimental/CDF. Acho que foi assim que o Huang observou a relação de psi e phi_m com Ar e Pr. 
    N = len(df)
    fluido = AbstractState(backend, fluido)
    P_s0, h_s0, s_s0, P_p0, h_p0, s_p0, Ar, Pr, w_exp, Pd_exp = pre_processar_dados(df)
    

    corrente_principal = gsim.MaterialStream(fluido)
    corrente_secundaria = gsim.MaterialStream(fluido)
    ejetor = gsim.Ejector(corrente_principal, corrente_secundaria, fluido)

    def funcao_objetivo(param, i):
        phi_m, psi = param

        corrente_secundaria.p, corrente_secundaria.h, corrente_secundaria.s = P_s0[i], h_s0[i], s_s0[i]
        corrente_principal.p, corrente_principal.h, corrente_principal.s = P_p0[i], h_p0[i], s_p0[i]
        ejetor.set_dimensions('d', df.loc[i, 'dt']/1000,
                                    df.loc[i, 'dp1']/1000,
                                    df.loc[i, 'd3']/1000)
        ejetor.set_efficiencies(0.95, 0.95, 0.95, phi_m, psi)

        try:
            ejetor.calculate()
            erro_w = ((ejetor.entrainment_ratio - w_exp[i]) / w_exp[i])**2
            erro_Pd = ((ejetor.Pd_crit - Pd_exp[i]) / Pd_exp[i])**2
            return erro_w + erro_Pd #quem eu quero zerar/minimizar
        except Exception:
            return 10.0  # mesmo esquema de antes, se a conta der errado o erro sovbe

    def funcao_restricoes(param):
        phi_m, psi = param
        return np.array([phi_m - 0.1, 0.99 - phi_m, psi - 0.1, 0.99 - psi])

    chute = [0.8, 0.8]
    restricoes = {'type': 'ineq', 'fun': funcao_restricoes}

    phi_m_opt = np.full(N, np.nan)
    psi_opt = np.full(N, np.nan)
    convergiu = np.zeros(N, dtype=bool)

    print("Otimizando phi_m e psi ponto a ponto...")
    for i in range(N):
        resultado = opt.minimize(funcao_objetivo, chute, args=(i,), method='SLSQP',
                                  constraints=restricoes, options={'maxiter': 300})
        if resultado.success:
            phi_m_opt[i], psi_opt[i] = resultado.x
            convergiu[i] = True
            print(f"Ponto {i+1:03d}: phi_m={phi_m_opt[i]:.4f} | psi={psi_opt[i]:.4f} | Erro={resultado.fun:.2e}")
        else:
            print(f"Ponto {i+1:03d}: falhou -> {resultado.message}")

    
    df_resultado = df.copy()
    df_resultado['Ar'] = Ar
    df_resultado['Pr'] = Pr
    df_resultado['phi_m_opt'] = phi_m_opt
    df_resultado['psi_opt'] = psi_opt

    ok = convergiu  # máscara só com pontos que convergiram

    print("\n============================================================")
    print("Média de phi_m ótimo: {:.4f}".format(np.mean(phi_m_opt)))
    print("Média de psi ótimo: {:.4f}".format(np.mean(psi_opt)))

    fig, axs = plt.subplots(2, 3, figsize=(12, 8))

    axs[0, 0].scatter(Ar[ok], phi_m_opt[ok])
    axs[0, 0].set_xlabel('Ar'); axs[0, 0].set_ylabel(r'$\phi_m$')

    axs[0, 1].scatter(Pr[ok], phi_m_opt[ok])
    axs[0, 1].set_xlabel('Pr'); axs[0, 1].set_ylabel(r'$\phi_m$')

    axs[0, 2].scatter(Ar[ok] * Pr[ok], phi_m_opt[ok])
    axs[0, 2].set_xlabel('Ar * Pr'); axs[0, 2].set_ylabel(r'$\phi_m$')

    axs[1, 0].scatter(Ar[ok], psi_opt[ok])
    axs[1, 0].set_xlabel('Ar'); axs[1, 0].set_ylabel(r'$\psi$')

    axs[1, 1].scatter(Pr[ok], psi_opt[ok])
    axs[1, 1].set_xlabel('Pr'); axs[1, 1].set_ylabel(r'$\psi$')

    axs[1, 2].scatter(Ar[ok]* Pr[ok] , psi_opt[ok])
    axs[1, 2].set_xlabel('Ar * Pr'); axs[1, 2].set_ylabel(r'$\psi$')

    plt.tight_layout()
    plt.show()