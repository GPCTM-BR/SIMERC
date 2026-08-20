import flet as ft
import CoolProp.CoolProp as CP
from CoolProp import AbstractState
import os
import sys
import multiprocessing
import asyncio
from concurrent.futures import ProcessPoolExecutor
import uuid
from gsim import Ejector, MaterialStream, Pump, HeaterCooler, Valve
import numpy as np
import itertools
import csv
import flet_datatable2 as fdt
import xlsxwriter
import pyarrow as pa
import pyarrow.parquet as pq

def obter_caminho_recurso(caminho_relativo):
    try:
        # O PyInstaller cria uma pasta temporária e armazena o caminho em _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, caminho_relativo)

def _write_csv_file(file_path, columns, rows):
    with open(file_path, "w", newline="", encoding="utf-8-sig") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(columns)
        writer.writerows(rows)

def _write_xlsx_file(file_path, columns, rows):
    workbook = xlsxwriter.Workbook(file_path, {"constant_memory": True})
    worksheet = workbook.add_worksheet("Batch Results")
 
    header_format = workbook.add_format({"bold": True, "bg_color": "#D9D9D9"})
    worksheet.write_row(0, 0, columns, header_format)
 
    for row_i, linha in enumerate(rows, start=1):
        worksheet.write_row(row_i, 0, [None if v != v else v for v in linha])
 
    for col_idx, col_name in enumerate(columns):
        worksheet.set_column(col_idx, col_idx, len(str(col_name)) + 2)
 
    workbook.close()

def _write_parquet_file(file_path, columns, rows, float_dtype=None):
    if float_dtype is None:
        float_dtype = pa.float64()
    
    n = len(rows)
    arrays = []
    for col_idx, nome_coluna in enumerate(columns):
        if nome_coluna == "Status":
            valores = [row[col_idx] for row in rows]
            arrays.append(pa.array(valores, type=pa.string()))
        else:
            buf = np.fromiter((row[col_idx] for row in rows), dtype=np.float64, count=n)
            arrays.append(pa.array(buf, type=float_dtype))
 
    tabela = pa.Table.from_arrays(arrays, names=columns)
    pq.write_table(tabela, file_path, compression="zstd", compression_level=1)

def ERC(fluid,
        p_eta,
        g_flash, g_par1, g_par2, g_eta, 
        c_flash, c_par1, c_par2, c_eta,
        e_flash, e_par1, e_par2, e_eta,
        dt, dp1, dconst, eta_t, eta_m, eta_d, phi_m, psi,
        batch=False, max_iter = [], tol = [], initial_guess = [], *batch_param):

    if g_flash != "PQ":
        g_par1, g_par2 = g_par2, g_par1

    
    if c_flash != "PQ":
        c_par1, c_par2 = c_par2, c_par1

    if e_flash != "PQ":
        e_par1, e_par2 = e_par2, e_par1


    flash_inputs ={
        "TP": CP.PT_INPUTS,
        "TQ": CP.QT_INPUTS,
        "PQ": CP.PQ_INPUTS
    }

    #definindo as correntes
    stream1 = MaterialStream(fluid)
    stream2 = MaterialStream(fluid)
    stream3 = MaterialStream(fluid)
    stream4 = MaterialStream(fluid)
    stream5 = MaterialStream(fluid)
    stream6 = MaterialStream(fluid)
    stream7 = MaterialStream(fluid)
    stream8 = MaterialStream(fluid)

    #definindo os componentes
    ejector = Ejector(stream1, stream2, fluid)
    if max_iter != []:
        ejector.max_iter_Pt, ejector.max_iter_Pp1, ejector.max_iter_Pconst, ejector.max_iter_rho4, ejector.max_iter_P5 = max_iter
    if tol != []:
        ejector.Pt_tol, ejector.Pp1_tol, ejector.Pconst_tol, ejector.rho4_tol, ejector.P5_tol = tol
    if initial_guess != []:
        Pt_guess, Pp1_guess, Pconst_guess, rho4_guess = initial_guess
    generator = HeaterCooler(stream7, stream1, fluid)
    condenser = HeaterCooler(stream3, stream4, fluid)
    evaporator = HeaterCooler(stream8, stream2, fluid)
    pump = Pump(stream5, stream7, fluid)
    expansion_valve = Valve(stream6, stream8, fluid)

    #calculando estado na saída do gerador, condensador e evaporador
    stream1.calculate(flash_inputs[g_flash], g_par1, g_par2)
    stream4.calculate(flash_inputs[c_flash], c_par1, c_par2)
    stream2.calculate(flash_inputs[e_flash], e_par1, e_par2)


    #calculando Ar e Pr (caso os coeficientes do ejetor use) e eval no phi_m e psi
    Ar = (dconst/dt)**2
    Pr = stream2.p/stream1.p
    
    try:
        phi_m = float(eval(phi_m))
    except:
        raise ValueError("The Mixing Loss Factor is invalid.")
    try:
        psi = float(eval(psi))
    except:
        raise ValueError("The Expansion Coefficient is invalid.")

    if phi_m>1 or phi_m<0:
        raise ValueError(f"The Mixing Loss Factor must be a number between 0 and 1. Current value: {phi_m:.4f}")
    if psi>1 or psi<=0:
        raise ValueError(f"The Expansion Coefficient must be a number between 0 and 1. Current value: {psi:.4f}")

    #calculando o ejetor
    ejector.set_dimensions("d", dt, dp1, dconst)
    ejector.set_efficiencies(eta_t, eta_m, eta_d, phi_m, psi)
    if batch:
        ejector.calculate(Pt_guess, Pp1_guess, Pconst_guess, rho4_guess)
    else:
        try:
            ejector.calculate()
        except:
            ejector.Pconst_tol = 1.5e-3
            ejector.calculate()
    
    #salvando mp, ms, entrainment ratio, Pd crit
    mp = ejector.m_p
    ms = ejector.m_s
    m = ejector.m
    entrainment_ratio= ejector.entrainment_ratio
    Pd_crit = ejector.Pd_crit
    
    P_lift_ratio = stream4.p/stream2.p

    #vverificando se está na região crítica
    if Pd_crit < stream4.p:
        raise ValueError("The condenser pressure is greater than the critical backpressure of the ejector. This ejector model is not valid for the subcritical region.")
    
    #definindo o fluxo das outras correntes
    stream1.mass_flow = mp
    stream2.mass_flow = ms
    stream3.mass_flow = m
    stream4.mass_flow = m
    stream5.mass_flow = mp
    stream6.mass_flow = ms
    stream7.mass_flow = mp
    stream8.mass_flow = ms
    #calculando stream3
    stream3.calculate(CP.HmassP_INPUTS, ejector.h_5, stream4.p)

    #calculando o condensador
    condenser.set_efficiency(c_eta, True)
    condenser.calculate(flash_inputs[c_flash], c_par1, c_par2) #isso definiu stream4
    condenser_heat = condenser.heat_duty

    #calculando stream 5 e 6
    stream5.calculate(CP.HmassP_INPUTS, stream4.h, stream4.p)
    stream6.calculate(CP.HmassP_INPUTS, stream4.h, stream4.p)

    #calculando a bomba
    pump.set_isentropic_efficiency(p_eta, True)
    pump.calculate(stream1.p) #isso definiiu stream7
    pump_work = pump.work

    #calculando a válvula
    expansion_valve.calculate(stream2.p) #isso definiu stream8

    #calculando o gerador
    generator.set_efficiency(g_eta, True)
    generator.calculate(flash_inputs[g_flash], g_par1, g_par2) #isso definiu stream 1
    generator_heat = generator.heat_duty

    #calculando o evaporador
    evaporator.set_efficiency(e_eta, True)
    evaporator.calculate(flash_inputs[e_flash], e_par1, e_par2) #isso definiu stream 2
    evaporator_heat = evaporator.heat_duty

    #verificando se Tger > Tcond > Tevap
    if not (stream1.T > stream4.T):
        raise ValueError("For the cycle to make sense, GENERATOR: Outlet Temperature must be greater than CONDENSER: Outlet Temperature")
    if not (stream4.T > stream2.T):
        raise ValueError("For the cycle to make sense, CONDENSER: Outlet Temperature must be greater than EVAPORATOR: Outlet Temperature")
    
    #calculando o COP
    COP = abs(evaporator_heat)/(abs(generator_heat) + abs(pump_work))

    #retorno das respostas
    if batch:
        batch_param_list = list(batch_param)
        near_result = [ejector.P_t, ejector.P_p1, ejector.P_const, ejector.rho_4]
        param_dict = {
            "COP": COP,
            "EJECTOR: Entrainment Ratio": entrainment_ratio,
            "EJECTOR: Pressure Lift Ratio": P_lift_ratio,
            "EJECTOR: Critical Backpressure": Pd_crit,
            "EJECTOR: Primary Inlet Mass Flow": stream1.mass_flow,
            "EJECTOR: Secondary Inlet Mass Flow": stream2.mass_flow,
            "EJECTOR: Outlet Mass Flow": stream3.mass_flow,
            "EJECTOR: Primary Inlet Temperature": stream1.T,
            "EJECTOR: Primary Inlet Pressure": stream1.p,
            "EJECTOR: Primary Inlet Specific Enthalpy": stream1.h,
            "EJECTOR: Primary Inlet Specific Entropy": stream1.s,
            "EJECTOR: Primary Inlet Density": stream1.rho,
            "EJECTOR: Secondary Inlet Temperature": stream2.T,
            "EJECTOR: Secondary Inlet Pressure": stream2.p,
            "EJECTOR: Secondary Inlet Specific Enthalpy": stream2.h,
            "EJECTOR: Secondary Inlet Specific Entropy": stream2.s,
            "EJECTOR: Secondary Inlet Density": stream2.rho,
            "EJECTOR: Outlet Temperature": stream3.T,
            "EJECTOR: Outlet Pressure": stream3.p,
            "EJECTOR: Outlet Specific Enthalpy": stream3.h,
            "EJECTOR: Outlet Specific Entropy": stream3.s,
            "EJECTOR: Outlet Density": stream3.rho,
            "EVAPORATOR: Heaty Duty": evaporator_heat,
            "EVAPORATOR: Mass Flow": stream8.mass_flow,
            "EVAPORATOR: Inlet Temperature": stream8.T,
            "EVAPORATOR: Inlet Pressure": stream8.p,
            "EVAPORATOR: Inlet Specific Enthalpy": stream8.h,
            "EVAPORATOR: Inlet Specific Entropy": stream8.s,
            "EVAPORATOR: Inlet Density": stream8.rho,
            "EVAPORATOR: Outlet Temperature": stream2.T,
            "EVAPORATOR: Outlet Pressure": stream2.p,
            "EVAPORATOR: Outlet Specific Enthalpy": stream2.h,
            "EVAPORATOR: Outlet Specific Entropy": stream2.s,
            "EVAPORATOR: Outlet Density": stream2.rho,
            "GENERATOR: Heaty Duty": generator_heat,
            "GENERATOR: Mass Flow": stream7.mass_flow,
            "GENERATOR: Inlet Temperature": stream7.T,
            "GENERATOR: Inlet Pressure": stream7.p,
            "GENERATOR: Inlet Specific Enthalpy": stream7.h,
            "GENERATOR: Inlet Specific Entropy": stream7.s,
            "GENERATOR: Inlet Density": stream7.rho,
            "GENERATOR: Outlet Temperature": stream1.T,
            "GENERATOR: Outlet Pressure": stream1.p,
            "GENERATOR: Outlet Specific Enthalpy": stream1.h,
            "GENERATOR: Outlet Specific Entropy": stream1.s,
            "GENERATOR: Outlet Density": stream1.rho,
            "CONDENSER: Heaty Duty": condenser_heat,
            "CONDENSER: Mass Flow": stream3.mass_flow,
            "CONDENSER: Inlet Temperature": stream3.T,
            "CONDENSER: Inlet Pressure": stream3.p,
            "CONDENSER: Inlet Specific Enthalpy": stream3.h,
            "CONDENSER: Inlet Specific Entropy": stream3.s,
            "CONDENSER: Inlet Density": stream3.rho,
            "CONDENSER: Outlet Temperature": stream4.T,
            "CONDENSER: Outlet Pressure": stream4.p,
            "CONDENSER: Outlet Specific Enthalpy": stream4.h,
            "CONDENSER: Outlet Specific Entropy": stream4.s,
            "CONDENSER: Outlet Density": stream4.rho,
            "EXPANSION VALVE: Mass Flow": stream6.mass_flow,
            "EXPANSION VALVE: Inlet Temperature": stream6.T,
            "EXPANSION VALVE: Inlet Pressure": stream6.p,
            "EXPANSION VALVE: Inlet Specific Enthalpy": stream6.h,
            "EXPANSION VALVE: Inlet Specific Entropy": stream6.s,
            "EXPANSION VALVE: Inlet Density": stream6.rho,
            "EXPANSION VALVE: Outlet Temperature": stream8.T,
            "EXPANSION VALVE: Outlet Pressure": stream8.p,
            "EXPANSION VALVE: Outlet Specific Enthalpy": stream8.h,
            "EXPANSION VALVE: Outlet Specific Entropy": stream8.s,
            "EXPANSION VALVE: Outlet Density": stream8.rho,
            "PUMP: Work": pump_work,
            "PUMP: Mass Flow": stream5.mass_flow,
            "PUMP: Inlet Temperature": stream5.T,
            "PUMP: Inlet Pressure": stream5.p,
            "PUMP: Inlet Specific Enthalpy": stream5.h,
            "PUMP: Inlet Specific Entropy": stream5.s,
            "PUMP: Inlet Density": stream5.rho,
            "PUMP: Outlet Temperature": stream7.T,
            "PUMP: Outlet Pressure": stream7.p,
            "PUMP: Outlet Specific Enthalpy": stream7.h,
            "PUMP: Outlet Specific Entropy": stream7.s,
            "PUMP: Outlet Density": stream7.rho,
            "SOLVER: Converged (all)": float(ejector.converged),
            "SOLVER: Max Residual": ejector.max_residual,
            "SOLVER: Residual P_t": ejector.residual_Pt,
            "SOLVER: Residual P_p1": ejector.residual_Pp1,
            "SOLVER: Residual P_const": ejector.residual_Pconst,
            "SOLVER: Residual rho_4": ejector.residual_rho4,
            "SOLVER: Residual P_5": ejector.residual_P5,
            "SOLVER: Iterations P_t": float(ejector.iter_Pt),
            "SOLVER: Iterations P_p1": float(ejector.iter_Pp1),
            "SOLVER: Iterations P_const": float(ejector.iter_Pconst),
            "SOLVER: Iterations rho_4": float(ejector.iter_rho4),
            "SOLVER: Iterations P_5": float(ejector.iter_P5),
        }
        solution_erc = {}
        for param in batch_param_list:
            solution_erc[param] = param_dict[param]
        return solution_erc, near_result
    else:
        solution_erc = {
            "stream1": stream1,
            "stream2": stream2,
            "stream3": stream3,
            "stream4": stream4,
            "stream5": stream5,
            "stream6": stream6,
            "stream7": stream7,
            "stream8": stream8,
            "entrainment_ratio": entrainment_ratio,
            "Pd_crit": Pd_crit,
            "P_lift_ratio": P_lift_ratio,
            "condenser_heat": condenser_heat,
            "generator_heat": generator_heat,
            "evaporator_heat": evaporator_heat,
            "pump_work": pump_work,
            "COP": COP
        }
        return solution_erc

def _run_batch_chunk(rows, backend_name, caminho_tabelas, max_iter, tol):
    try:
        CP.set_config_string(CP.ALTERNATIVE_TABLES_DIRECTORY, caminho_tabelas + os.sep)
    except Exception:
        pass
    fluid = None
    fluid_name_atual = None
    resultados_chunk = []
    chute = [0.0,0.0,0.0,0.0] #se for zero ele usa o chute interno
    for row in rows:
        try:
            fluid_name = row[0]
            if fluid is None or fluid_name != fluid_name_atual:
                fluid = AbstractState(backend_name, fluid_name)
                fluid_name_atual = fluid_name
            rest = row[1:22]
            dependent_params = row[23:]
            try:
                resultado, near_result = ERC(fluid, *rest, True, max_iter, tol, list(chute), *dependent_params)
                chute = near_result
            except Exception:
                chute = [0.0,0.0,0.0,0.0]
                resultado, near_result = ERC(fluid, *rest, True, max_iter, tol, list(chute), *dependent_params)
            resultados_chunk.append(resultado)
        except Exception as exc:
            resultados_chunk.append({"error": str(exc)})
    return resultados_chunk

def fluid_list():
        string_fluids = CP.get_global_param_string('FluidsList')
        list_fluids = string_fluids.split(',')
        dropdown_list =[]
        for fluid in list_fluids:
             dropdown_list.append(ft.DropdownOption(key=fluid, text=fluid, style=ft.ButtonStyle(color="#B7C3C1")))
        return dropdown_list

async def main(page: ft.Page):

    VERSAO = "1.0" #só pra fazer um charme kkkkkk
    page.window.maximized = True
    page.bgcolor = "#011314"
    page.padding = ft.Padding(10,10,10,0)
    page.spacing= 0
    page.fonts = {"Open Sans Regular": "/fonts/OpenSans-Regular.ttf",
                  "Open Sans Light": "/fonts/OpenSans-Light.ttf"}
    page.theme = ft.Theme(font_family='Open Sans Light')
    page.theme_mode = ft.ThemeMode.DARK
    page.title = "SIMERC"
    page.window.icon ="icon.ico"

    painel_aberto = False
    LARGURA_PAINEL = 500
    SPACING = 5
    SIMULATED = False #vc nao pode abrir as aba Batch sim  se nao tivwer simulado antes
    BATCHED = False #vc nao pode abrir a aba Batch sim table se nao tivwer feito um batch
    batch_cancel_requested = False #caso o cara clique em cancelar, isso vira True
    batch_results = None
    
    batch_independent_labels = None
    batch_independent_values = None
    batch_dependent_params = None
    batch_table_columns = None
    batch_table_rows = None
    current_page = 0
    batch_independent_si_columns = None    # lista de colunas (uma lista por eixo), cada uma com o valor SI de cada linha
    batch_independent_units = None         # unidade de exibição de cada eixo (lista, mesma ordem das colunas independentes)
    PAGE_SIZE = 100
    
    base = os.path.dirname(sys.executable)
    caminho_tabelas = os.path.join(base, "CoolPropTables")
    os.makedirs(caminho_tabelas, exist_ok=True)
    CP.set_config_string(CP.ALTERNATIVE_TABLES_DIRECTORY, caminho_tabelas + os.sep)

    def tooltip_msg(msg: str, duration = 1500):
        return ft.Tooltip(msg, wait_duration=duration, vertical_offset=30, decoration=ft.BoxDecoration(border = ft.Border.all(1,"#50B7C3C1"),border_radius=2, bgcolor="#011314"), text_style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.2, size=11))

    def pegar_tol_iter():
        tol = []
        iter =[]
        for row in tol_iter_column.controls:
            if "Tolerance" in row.controls[1].value:
                try:
                    n = float(eval(row.controls[2].value))
                    tol.append(n)
                except:
                    row.controls[2].value = "1.5e-08"
                    n = float(eval(row.controls[2].value))
                    tol.append(n)
            else:
                try:
                    n = float(eval(row.controls[2].value))
                    if not n.is_integer():
                        raise ValueError()
                    iter.append(n)
                except:
                    row.controls[2].value = "60"
                    n = int(eval(row.controls[2].value))
                    iter.append(n)
        return tol, iter

    async def batch_sim(e):
        if not SIMULATED:
            await mostrar_erro("It is necessary to simulate a base case first.")
        else:
            if painel_aberto:
                animar_janela(e)
            
            if coluna.controls[1].content.key == "batch":
                coluna.controls[1].content = ft.Container(
                        key ="erc",
                        content=erc,
                        alignment=ft.Alignment.CENTER,
                        expand=True,
                        data='out_click',
                        on_click=gerenciar_clique,
                        image= ft.DecorationImage(src="/images/pattern3.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=20, opacity=0.5)
                    )
            else:
                #desceleciona todo mundo (desceleciona ou deseleciona? aff preciso durmir)
                for i in range(1, 7):
                    erc.controls[i].gradient.colors = ['#00B3F364', '#00FFFFFF']
                    erc.controls[i].border = ft.Border.all(0.5, "#00B6FF57")
                nonlocal last_clicked
                last_clicked = ''
                #troca o conteudo pelo conteudo de batch_sim
                coluna.controls[1].content=ft.Container(
                    key = "batch",
                    expand=True,
                    content=batch_sim_content,
                    image= ft.DecorationImage(src="/images/pattern4.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=50, opacity=.2)
                )
            page.update()
    
    async def results(e):
        if not SIMULATED:
            await mostrar_erro("It is necessary to simulate a base case first.")
        elif not BATCHED:
            await mostrar_erro("There are no Batch Sim results to show.")
        else:
            if painel_aberto:
                animar_janela(e)
            
            if coluna.controls[1].content.key == "results":
                coluna.controls[1].content = ft.Container(
                        key ="erc",
                        content=erc,
                        alignment=ft.Alignment.CENTER,
                        expand=True,
                        data='out_click',
                        on_click=gerenciar_clique,
                        image= ft.DecorationImage(src="/images/pattern3.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=20, opacity=0.5)
                    )
            else:
                #desceleciona todo mundo (desceleciona ou deseleciona? aff preciso durmir)
                for i in range(1, 7):
                    erc.controls[i].gradient.colors = ['#00B3F364', '#00FFFFFF']
                    erc.controls[i].border = ft.Border.all(0.5, "#00B6FF57")
                nonlocal last_clicked
                last_clicked = ''
                coluna.controls[1].content=ft.Container(
                    key = "results",
                    expand=True,
                    content=results_content,
                    image= ft.DecorationImage(src="/images/pattern4.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=50, opacity=.2)
                )
            page.update()

    async def run(e):
        if coluna.controls[1].content.key == "batch":
            await run_batch_ERC(e)
            return
        else:
            await run_single_ERC(e)
            return
 
    async def fechar_settings(e):
        overlay_settings.opacity = 0
        page.update()
        await asyncio.sleep(0.3)
        overlay_settings.visible = False
        page.update()
    
    async def abrir_settings(e):
        overlay_settings.visible = True
        overlay_settings.opacity = 0
        page.update()
        await asyncio.sleep(0.05)
        overlay_settings.opacity = 1
        page.update()
    
    close_settings_button = ft.IconButton(
        expand = False,
        tooltip=tooltip_msg(f"Close Settings", 100),
        icon=ft.Icon(icon= ft.Icons.CLOSE_ROUNDED, color="#FF5C5C"),
        height=48,   # mesma altura do dropdown
        width=48,

        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(
                radius=ft.BorderRadius(top_left=10, bottom_left=0, top_right=0, bottom_right=10)
            ),
            
        ),
        on_click=fechar_settings
    )
    
    settings_title = ft.Container(expand = True, content = ft.Row(expand = False, controls=[
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2, margin=ft.Margin.only(left=0)),
                        ft.Text("SIMULATION SETTINGS", style=ft.TextStyle(font_family="Open Sans Light", size=20, letter_spacing=2, color="#B7C3C1")),
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2, margin=ft.Margin.only(right=20)),
                    ]))
    
    batch_parallel_checkbox = ft.Checkbox("Run batch in parallel", value=False,
                                          tooltip=tooltip_msg("Defines whether Batch Sim will use parallel processing. Useful when there are many simulations."),
                                          label_style= ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1"),
                                          active_color= "#B6FF57",
                                          check_color="#011314"
                                          )
    
    def reset_settings(e):
        batch_parallel_checkbox.value= False
        for row in tol_iter_column.controls:
            if "Tolerance" in row.controls[1].value:
                row.controls[2].value = "1.5e-08"
            else:
                row.controls[2].value = "60"

    reset_settings_button = ft.Button(content= "Reset Settings", on_click=reset_settings,
                                      bgcolor="#B6FF57", color="#011314",
                                      style=ft.ButtonStyle(text_style=ft.TextStyle(font_family="Open Sans Regular", letter_spacing=1.5, size=16),
                                                           shape=ft.RoundedRectangleBorder(radius=5)))

    def check_config(e):
        for row in tol_iter_column.controls:
            if "Tolerance" in row.controls[1].value:
                try:
                    n = float(eval(row.controls[2].value))
                    if n == 0:
                        raise ValueError()
                except:
                    row.controls[2].value = "1.5e-08"
            else:
                try:
                    n = float(eval(row.controls[2].value))
                    if not n.is_integer():
                        raise ValueError()
                except:
                    row.controls[2].value = "60"

    def tol_iter_row(msg, valor):
        return ft.Row(
            controls=[
                    ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                    ft.Text(f"{msg}", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),
                    ft.TextField(text_align=ft.TextAlign.RIGHT,value= f"{valor}",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5, input_filter=ft.InputFilter(allow=True, regex_string=r"^[+-]?(\d+(\.\d*)?|\.\d*)?([eE][+-]?\d*)?$"), on_blur=check_config)
            ]
        )
    
    tol_iter_column = ft.Column(
        margin = ft.Margin.only(left=10,right=10),
        spacing=10,
        scroll=ft.ScrollMode.AUTO,
        expand=True,
        controls=[
            tol_iter_row("Primary Nozzle Throat: Loop Tolerance", 1.5e-8),
            tol_iter_row("Primary Nozzle Exit: Loop Tolerance", 1.5e-8),
            tol_iter_row("Constant Area Section: Loop Tolerance", 1.5e-8),
            tol_iter_row("Shock: Loop Tolerance", 1.5e-8),
            tol_iter_row("Diffuser: Loop Tolerance", 1.5e-8),
            tol_iter_row("Primary Nozzle Throat: Max Iterations", 60),
            tol_iter_row("Primary Nozzle Exit: Max Iterations", 60),
            tol_iter_row("Constant Area Section: Max Iterations", 60),
            tol_iter_row("Shock: Max Iterations", 60),
            tol_iter_row("Diffuser: Max Iterations", 60),
        ]
    )
    
    async def scroll_baixo(e):
        await tol_iter_column.scroll_to(delta=100, duration=500)

    overlay_settings = ft.Container(
        bgcolor="#50000000",
        expand=True,
        visible=False,
        blur=10,
        opacity=0,
        alignment=ft.Alignment.CENTER,
        animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT),
        content=ft.Container(
            alignment=ft.Alignment.CENTER,
            height=600/830.4*page.window.height,
            width=800/1550.4*page.window.width,
            bgcolor="#0A2022",
            border_radius=10,
            border=ft.Border.all(1, "#50B7C3C1"),
            padding=0,
            content=ft.Column(
                alignment=ft.MainAxisAlignment.START,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=20,
                expand = True,
                controls=[
                    ft.Row(expand=False, controls=[close_settings_button, settings_title], alignment=ft.MainAxisAlignment.START),
                    ft.Column(
                        margin = 10,
                        expand =True,
                        spacing = 30,
                        controls = [
                            ft.Row(margin= ft.Margin.only(right=10),expand = False, controls = [ft.Container(expand = True, content = batch_parallel_checkbox, alignment=ft.Alignment.CENTER_LEFT), ft.Container(expand = True, content = reset_settings_button, alignment=ft.Alignment.CENTER_RIGHT)]),
                            tol_iter_column,
                            ft.Row(controls=[ft.IconButton(icon = ft.Icons.ARROW_DROP_DOWN_ROUNDED, icon_color="#B6FF57", on_click= scroll_baixo)], alignment=ft.MainAxisAlignment.CENTER)
                        ]
                        
                    )
                ],
            ),
        ),
    )

    async def run_single_ERC(e):
        #COMENTÁRIO PÓS ESCRITA: como é que eu tive paciencia pra escrever isso eu não sei
        page.update()
        #verifica se o fluido é válido
        if dropdown.value not in CP.get_global_param_string('FluidsList').split(","):
            await mostrar_erro('Please select a valid fluid')
            return
        #verifica se o backend foi selecionado
        if backend.value != "HEOS" and backend.value != "BICUBIC&HEOS":
            await mostrar_erro('Please select a backend')
            return
        #se BICUBIC, faz o bake (caaso já não tenha feito)
        if backend.value == 'BICUBIC&HEOS':
            await bake_fluid(e, run = True)   


        #checando input da bomba ===============================
        p_eta = pump_column_controls[2].content.controls[1].content.value
        try:
            p_eta = float(eval(p_eta))
        except:
            await mostrar_erro("The isentropic efficiency of the pump is invalid.")
            return
        if p_eta > 100 or p_eta <= 0:
            await mostrar_erro("The isentropic efficiency of the pump should be a number between 0% and 100%.")
            return
        #=======================================================

        #checando input do gerador =============================
        g_eta = generator_column_controls[2].controls[1].controls[3].content.value
        g_flash = g_flash_type_dropdown.value
        g_par1 = generator_column_controls[2].controls[1].controls[1].initial_value
        try:
            g_eta = float(eval(g_eta))
        except:
            await mostrar_erro("The efficiency of the generator is invalid.")
            return
        if g_eta > 100 or g_eta <= 0:
            await mostrar_erro("The efficiency of the generator should be a number between 0% and 100%.")
            return
        try:
            g_par1 = float(eval(g_par1))
        except:
            await mostrar_erro(f'Generator {heater_cooler_flash[g_flash][0]} is invalid')
            return
        if (g_par1 < 0 and g_flash == "PQ") or (g_par1 < 0 and generator_column_controls[2].controls[1].controls[1].controls[1].value == "K"):
            await mostrar_erro(f"Generator {heater_cooler_flash[g_flash][0]} can't be negative (absolute {heater_cooler_flash[g_flash][0].split()[1].lower()})")
            return
        if g_flash == "TP":
            g_par2 = generator_column_controls[2].controls[1].controls[2].initial_value
        else:
            g_par2 = generator_column_controls[2].controls[1].controls[2].content.value
        try:
            g_par2 = float(eval(g_par2))
        except:
            await mostrar_erro(f'Generator {heater_cooler_flash[g_flash][1]} is invalid')
            return
        if g_par2 < 0 and g_flash == "TP":
            await mostrar_erro(f"Generator {heater_cooler_flash[g_flash][1]} can't be negative (absolute {heater_cooler_flash[g_flash][1].split()[1].lower()})")
            return
        if (g_flash != "TP") and (g_par2>1 or g_par2<0):
            await mostrar_erro("Generator Outlet Vapor Fraction must be a number between 0 and 1.")
            return
        #=======================================================
        
        #checando input do condensador =============================
        c_eta = condenser_column_controls[2].controls[1].controls[3].content.value
        c_flash = c_flash_type_dropdown.value
        c_par1 = condenser_column_controls[2].controls[1].controls[1].initial_value
        try:
            c_eta = float(eval(c_eta))
        except:
            await mostrar_erro("The efficiency of the condenser is invalid.")
            return
        if c_eta > 100 or c_eta <= 0:
            await mostrar_erro("The efficiency of the condenser should be a number between 0% and 100%.")
            return
        try:
            c_par1 = float(eval(c_par1))
        except:
            await mostrar_erro(f'Condenser {heater_cooler_flash[c_flash][0]} is invalid')
            return
        if (c_par1 < 0 and c_flash == "PQ") or (c_par1 < 0 and condenser_column_controls[2].controls[1].controls[1].controls[1].value == "K"):
            await mostrar_erro(f"Condenser {heater_cooler_flash[c_flash][0]} can't be negative (absolute {heater_cooler_flash[c_flash][0].split()[1].lower()})")
            return
        if c_flash == "TP":
            c_par2 = condenser_column_controls[2].controls[1].controls[2].initial_value
        else:
            c_par2 = condenser_column_controls[2].controls[1].controls[2].content.value
        try:
            c_par2 = float(eval(c_par2))
        except:
            await mostrar_erro(f'Condenser {heater_cooler_flash[c_flash][1]} is invalid')
            return
        if c_par2 < 0 and c_flash == "TP":
            await mostrar_erro(f"Condenser {heater_cooler_flash[c_flash][1]} can't be negative (absolute {heater_cooler_flash[c_flash][1].split()[1].lower()})")
            return
        if (c_flash != "TP") and (c_par2>1 or c_par2<0):
            await mostrar_erro("Condenser Outlet Vapor Fraction must be a number between 0 and 1.")
            return
        #=======================================================

        #checando input do evaporador =============================
        e_eta = evaporator_column_controls[2].controls[1].controls[3].content.value
        e_flash = e_flash_type_dropdown.value
        e_par1 = evaporator_column_controls[2].controls[1].controls[1].initial_value
        try:
            e_eta = float(eval(e_eta))
        except:
            await mostrar_erro("The efficiency of the evaporator is invalid.")
            return
        if e_eta > 100 or e_eta <= 0:
            await mostrar_erro("The efficiency of the evaporator should be a number between 0% and 100%.")
            return
        try:
            e_par1 = float(eval(e_par1))
        except:
            await mostrar_erro(f'Evaporator {heater_cooler_flash[e_flash][0]} is invalid')
            return
        if (e_par1 < 0 and e_flash == "PQ") or (e_par1 < 0 and evaporator_column_controls[2].controls[1].controls[1].controls[1].value == "K"):
            await mostrar_erro(f"Evaporator {heater_cooler_flash[e_flash][0]} can't be negative (absolute {heater_cooler_flash[e_flash][0].split()[1].lower()})")
            return
        if e_flash == "TP":
            e_par2 = evaporator_column_controls[2].controls[1].controls[2].initial_value
        else:
            e_par2 = evaporator_column_controls[2].controls[1].controls[2].content.value
        try:
            e_par2 = float(eval(e_par2))
        except:
            await mostrar_erro(f'Evaporator {heater_cooler_flash[e_flash][1]} is invalid')
            return
        if e_par2 < 0 and e_flash == "TP":
            await mostrar_erro(f"Evaporator {heater_cooler_flash[e_flash][1]} can't be negative (absolute {heater_cooler_flash[e_flash][1].split()[1].lower()})")
            return
        if (e_flash != "TP") and (e_par2>1 or e_par2<0):
            await mostrar_erro("Evaporator Outlet Vapor Fraction must be a number between 0 and 1.")
            return
        #=======================================================

        #checando input do ejetor ==============================
        dt = ejector_column_controls[2].controls[1].controls[0].initial_value
        dp1 = ejector_column_controls[2].controls[1].controls[1].initial_value
        dconst = ejector_column_controls[2].controls[1].controls[2].initial_value
        eta_t = ejector_column_controls[2].controls[1].controls[3].content.value
        eta_m = ejector_column_controls[2].controls[1].controls[4].content.value
        eta_d = ejector_column_controls[2].controls[1].controls[5].content.value
        phi_m = ejector_column_controls[2].controls[1].controls[6].content.value
        psi = ejector_column_controls[2].controls[1].controls[7].content.value
        try:
            dt = float(eval(dt))
        except:
            await mostrar_erro("Primary Nozzle Throat Diameter is invalid")
            return
        if dt <= 0:
            await mostrar_erro("Primary Nozzle Throat Diameter must be greater than zero.")
            return
        try:
            dp1 = float(eval(dp1))
        except:
            await mostrar_erro("Primary Nozzle Exit Diameter is invalid")
            return
        if dp1 <= 0:
            await mostrar_erro("Primary Nozzle Exit Diameter must be greater than zero.")
            return
        try:
            dconst = float(eval(dconst))
        except:
            await mostrar_erro("Constant Area Section Diameter is invalid")
            return
        if dconst <= 0:
            await mostrar_erro("Constant Area Section Diameter must be greater than zero.")
            return
        try:
            eta_t = float(eval(eta_t))
        except:
            await mostrar_erro("The Throat Isentropic Efficiency is invalid.")
            return
        if eta_t > 100 or eta_t <= 0:
            await mostrar_erro("The Throat Isentropic Efficiency should be a number between 0 and 1.")
            return
        try:
            eta_m = float(eval(eta_m))
        except:
            await mostrar_erro("The Aerodynamic Throat Isentropic Efficiency is invalid.")
            return
        if eta_m > 100 or eta_m <= 0:
            await mostrar_erro("The Aerodynamic Throat Isentropic Efficiency should be a number between 0 and 1.")
            return
        try:
            eta_d = float(eval(eta_d))
        except:
            await mostrar_erro("The Diffuser Isentropic Efficiency is invalid.")
            return
        if eta_d > 100 or eta_d <= 0:
            await mostrar_erro("The Diffuser Isentropic Efficiency should be a number between 0 and 1.")
            return
        #=======================================================

        #iniciando a simulação
        fluido = AbstractState(backend.value, dropdown.value)

        g_par1 = return_unit_to_SI(g_par1, generator_column_controls[2].controls[1].controls[1].unit)
        if g_flash =="TP":
            g_par2 = return_unit_to_SI(g_par2, generator_column_controls[2].controls[1].controls[2].unit)

        c_par1 = return_unit_to_SI(c_par1, condenser_column_controls[2].controls[1].controls[1].unit)
        if c_flash =="TP":
            c_par2 = return_unit_to_SI(c_par2, condenser_column_controls[2].controls[1].controls[2].unit)
        
        e_par1 = return_unit_to_SI(e_par1, evaporator_column_controls[2].controls[1].controls[1].unit)
        if e_flash =="TP":
            e_par2 = return_unit_to_SI(e_par2, evaporator_column_controls[2].controls[1].controls[2].unit)

        dt = return_unit_to_SI(dt, ejector_column_controls[2].controls[1].controls[0].unit)
        dp1 = return_unit_to_SI(dp1, ejector_column_controls[2].controls[1].controls[1].unit)
        dconst = return_unit_to_SI(dconst, ejector_column_controls[2].controls[1].controls[2].unit)

        try:
            tol, max_iter = pegar_tol_iter()
            results = ERC(fluido, p_eta, 
                g_flash, g_par1, g_par2, g_eta, 
                c_flash, c_par1, c_par2,c_eta, 
                e_flash, e_par1, e_par2, e_eta, 
                dt, dp1, dconst, eta_t, eta_m, eta_d, phi_m, psi,
                False, max_iter, tol)
            await mostrar_erro(f"Simulation completed successfully!", icone= ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size = 50))
            
            #Anotando os resultados nos paineis esquerdos
            #Bomba================
            pump_column_controls[4].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 21.5, controls=[
                    texto_na_esquerda("Work"),
                    texto_na_esquerda("Mass Flow"),
                    ft.Text("INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Inlet Temperature"),
                    texto_na_esquerda("Inlet Pressure"),
                    texto_na_esquerda("Inlet Specific Enthalpy"),
                    texto_na_esquerda("Inlet Specific Entropy"),
                    texto_na_esquerda("Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy"),
                    texto_na_esquerda("Outlet Specific Entropy"),
                    texto_na_esquerda("Outlet Density"),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['pump_work']:.4f}","W", "kW", "MW", "hp"),
                    ValorUnidade2(f"{results['stream5'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream5'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream5'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream5'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream5'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream5'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream7'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream7'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream7'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream7'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream7'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])
                
            #Válvula de Expansão================
            expansion_valve_column_controls[2].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 21.5, controls=[
                    texto_na_esquerda("Mass Flow"),
                    ft.Text("INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Inlet Temperature"),
                    texto_na_esquerda("Inlet Pressure"),
                    texto_na_esquerda("Inlet Specific Enthalpy"),
                    texto_na_esquerda("Inlet Specific Entropy"),
                    texto_na_esquerda("Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy"),
                    texto_na_esquerda("Outlet Specific Entropy"),
                    texto_na_esquerda("Outlet Density"),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['stream6'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream6'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream6'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream6'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream6'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream6'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream8'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream8'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream8'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream8'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream8'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])
                
            #Gerador================
            generator_column_controls[4].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 21.5, controls=[
                    texto_na_esquerda("Heat Duty"),
                    texto_na_esquerda("Mass Flow"),
                    ft.Text("INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Inlet Temperature"),
                    texto_na_esquerda("Inlet Pressure"),
                    texto_na_esquerda("Inlet Specific Enthalpy"),
                    texto_na_esquerda("Inlet Specific Entropy"),
                    texto_na_esquerda("Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy"),
                    texto_na_esquerda("Outlet Specific Entropy"),
                    texto_na_esquerda("Outlet Density"),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['generator_heat']:.4f}","W", "kW", "MW", "hp"),
                    ValorUnidade2(f"{results['stream7'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream7'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream7'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream7'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream7'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream7'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream1'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream1'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream1'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream1'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream1'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])

            #Condensador================
            condenser_column_controls[4].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 21.5, controls=[
                    texto_na_esquerda("Heat Duty"),
                    texto_na_esquerda("Mass Flow"),
                    ft.Text("INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Inlet Temperature"),
                    texto_na_esquerda("Inlet Pressure"),
                    texto_na_esquerda("Inlet Specific Enthalpy"),
                    texto_na_esquerda("Inlet Specific Entropy"),
                    texto_na_esquerda("Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy"),
                    texto_na_esquerda("Outlet Specific Entropy"),
                    texto_na_esquerda("Outlet Density"),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['condenser_heat']:.4f}","W", "kW", "MW", "hp"),
                    ValorUnidade2(f"{results['stream3'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream3'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream3'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream3'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream3'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream3'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream4'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream4'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream4'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream4'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream4'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])

            #evaporador================
            evaporator_column_controls[4].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 21.5, controls=[
                    texto_na_esquerda("Heat Duty"),
                    texto_na_esquerda("Mass Flow"),
                    ft.Text("INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Inlet Temperature"),
                    texto_na_esquerda("Inlet Pressure"),
                    texto_na_esquerda("Inlet Specific Enthalpy"),
                    texto_na_esquerda("Inlet Specific Entropy"),
                    texto_na_esquerda("Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy"),
                    texto_na_esquerda("Outlet Specific Entropy"),
                    texto_na_esquerda("Outlet Density"),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['evaporator_heat']:.4f}","W", "kW", "MW", "hp"),
                    ValorUnidade2(f"{results['stream8'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream8'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream8'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream8'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream8'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream8'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=15)),
                    ValorUnidade2(f"{results['stream2'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream2'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream2'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream2'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream2'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])
            
            #ejetor================
            ejector_column_controls[4].content = ft.Row(expand=True, controls=[
                ft.Column(expand=True,spacing = 20, controls=[
                    texto_na_esquerda("Primary Inlet Mass Flow", margin=ft.Margin.only(top=0)),
                    texto_na_esquerda("Secondary Inlet Mass Flow"),
                    texto_na_esquerda("Outlet Mass Flow"),
                    texto_na_esquerda("Entrainment Ratio", margin=ft.Margin.only(top=10)),
                    texto_na_esquerda("Pressure Lift Ratio", margin=ft.Margin.only(top=5)),
                    texto_na_esquerda("Critical Backpressure", margin=ft.Margin.only(top=5)),
                    ft.Text("PRIMARY INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=6)),
                    texto_na_esquerda("Primary Inlet Temperature"),
                    texto_na_esquerda("Primary Inlet Pressure"),
                    texto_na_esquerda("Primary Inlet Specific Enthalpy"),
                    texto_na_esquerda("Primary Inlet Specific Entropy"),
                    texto_na_esquerda("Primary Inlet Density"),
                    ft.Text("SECONDARY INLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=10)),
                    texto_na_esquerda("Secondary Inlet Temperature"),
                    texto_na_esquerda("Secondary Inlet Pressure"),
                    texto_na_esquerda("Secondary Inlet Specific Enthalpy"),
                    texto_na_esquerda("Secondary Inlet Specific Entropy"),
                    texto_na_esquerda("Secondary Inlet Density"),
                    ft.Text("OUTLET RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=10)),
                    texto_na_esquerda("Outlet Temperature"),
                    texto_na_esquerda("Outlet Pressure"),
                    texto_na_esquerda("Outlet Specific Enthalpy",margin=ft.Margin.only(top=5)),
                    texto_na_esquerda("Outlet Specific Entropy",margin=ft.Margin.only(top=5)),
                    texto_na_esquerda("Outlet Density", margin=ft.Margin.only(top=5)),
                    ]),
                ft.Column(expand=True, controls=[
                    ValorUnidade2(f"{results['stream1'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ValorUnidade2(f"{results['stream2'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ValorUnidade2(f"{results['stream3'].mass_flow:.4f}","kg/s","kg/h", "g/s", "lb/s", "lb/h"),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= f"{results["entrainment_ratio"]:.4f}",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5, read_only=True),),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= f"{results["P_lift_ratio"]:.4f}",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5, read_only=True),),
                    ValorUnidade2(f"{results['Pd_crit']:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=22)),
                    ValorUnidade2(f"{results['stream1'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream1'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream1'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream1'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream1'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=20)),
                    ValorUnidade2(f"{results['stream2'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream2'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream2'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream2'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream2'].rho:.4f}","kg/m³", "g/cm³"),
                    ft.Text("", size=10, margin=ft.Margin.only(top=20)),
                    ValorUnidade2(f"{results['stream3'].T:.4f}","K", "°C", "°F"),
                    ValorUnidade2(f"{results['stream3'].p:.4f}","Pa", "kPa", "MPa", "bar", "atm", "psi"),
                    ValorUnidade2(f"{results['stream3'].h:.4f}","J/kg", "kJ/kg", "MJ/kg", "Btu/lb"),
                    ValorUnidade2(f"{results['stream3'].s:.4f}","J/kg/K", "kJ/kg/K", "MJ/kg/K"),
                    ValorUnidade2(f"{results['stream3'].rho:.4f}","kg/m³", "g/cm³"),
                    ])
                ])
            nonlocal SIMULATED
            SIMULATED = True
            return
        except Exception as error:
            await mostrar_erro(f"{error}")
            return
    
    async def run_batch_ERC(e):
        try:
            #verifica se o fluido é válido
            if dropdown.value not in CP.get_global_param_string('FluidsList').split(","):
                await mostrar_erro('Please select a valid fluid')
                return
            #verifica se o backend foi selecionado
            if backend.value != "HEOS" and backend.value != "BICUBIC&HEOS":
                await mostrar_erro('Please select a backend')
                return
            #se BICUBIC, faz o bake (caaso já não tenha feito)
            if backend.value == 'BICUBIC&HEOS':
                await bake_fluid(e, run = True)  

            #pegando todas as entradas do single ERC
            p_eta = pump_column_controls[2].content.controls[1].content.value
            p_eta = float(eval(p_eta))

            g_eta = generator_column_controls[2].controls[1].controls[3].content.value
            g_flash = g_flash_type_dropdown.value
            g_par1 = generator_column_controls[2].controls[1].controls[1].initial_value
            g_eta = float(eval(g_eta))
            g_par1 = float(eval(g_par1))
            if g_flash == "TP":
                g_par2 = generator_column_controls[2].controls[1].controls[2].initial_value
            else:
                g_par2 = generator_column_controls[2].controls[1].controls[2].content.value
            g_par2 = float(eval(g_par2))

            c_eta = condenser_column_controls[2].controls[1].controls[3].content.value
            c_flash = c_flash_type_dropdown.value
            c_par1 = condenser_column_controls[2].controls[1].controls[1].initial_value
            c_eta = float(eval(c_eta))
            c_par1 = float(eval(c_par1))
            if c_flash == "TP":
                c_par2 = condenser_column_controls[2].controls[1].controls[2].initial_value
            else:
                c_par2 = condenser_column_controls[2].controls[1].controls[2].content.value
            c_par2 = float(eval(c_par2))

            e_eta = evaporator_column_controls[2].controls[1].controls[3].content.value
            e_flash = e_flash_type_dropdown.value
            e_par1 = evaporator_column_controls[2].controls[1].controls[1].initial_value
            e_eta = float(eval(e_eta))
            e_par1 = float(eval(e_par1))
            if e_flash == "TP":
                e_par2 = evaporator_column_controls[2].controls[1].controls[2].initial_value
            else:
                e_par2 = evaporator_column_controls[2].controls[1].controls[2].content.value
            e_par2 = float(eval(e_par2))

            dt = ejector_column_controls[2].controls[1].controls[0].initial_value
            dp1 = ejector_column_controls[2].controls[1].controls[1].initial_value
            dconst = ejector_column_controls[2].controls[1].controls[2].initial_value
            eta_t = ejector_column_controls[2].controls[1].controls[3].content.value
            eta_m = ejector_column_controls[2].controls[1].controls[4].content.value
            eta_d = ejector_column_controls[2].controls[1].controls[5].content.value
            phi_m = ejector_column_controls[2].controls[1].controls[6].content.value
            psi = ejector_column_controls[2].controls[1].controls[7].content.value
            dt = float(eval(dt))
            dp1 = float(eval(dp1))
            dconst = float(eval(dconst))
            eta_t = float(eval(eta_t))
            eta_m = float(eval(eta_m))
            eta_d = float(eval(eta_d))

            #retornando todo mundo pro SI
            g_par1 = return_unit_to_SI(g_par1, generator_column_controls[2].controls[1].controls[1].unit)
            if g_flash =="TP":
                g_par2 = return_unit_to_SI(g_par2, generator_column_controls[2].controls[1].controls[2].unit)

            c_par1 = return_unit_to_SI(c_par1, condenser_column_controls[2].controls[1].controls[1].unit)
            if c_flash =="TP":
                c_par2 = return_unit_to_SI(c_par2, condenser_column_controls[2].controls[1].controls[2].unit)
            
            e_par1 = return_unit_to_SI(e_par1, evaporator_column_controls[2].controls[1].controls[1].unit)
            if e_flash =="TP":
                e_par2 = return_unit_to_SI(e_par2, evaporator_column_controls[2].controls[1].controls[2].unit)

            dt = return_unit_to_SI(dt, ejector_column_controls[2].controls[1].controls[0].unit)
            dp1 = return_unit_to_SI(dp1, ejector_column_controls[2].controls[1].controls[1].unit)
            dconst = return_unit_to_SI(dconst, ejector_column_controls[2].controls[1].controls[2].unit)

            #pegando os outputs que o usuário escolheu
            dependent_params =[dependent_row.controls[0].value for dependent_row in dependent_rows_column.controls]
            if not dependent_params:
                await mostrar_erro("Select at least one dependent variable to track.")
                return
            

            #criando array de tamanho N
            caso_base = [
                dropdown.value, p_eta, g_flash, g_par1, g_par2, g_eta, 
                c_flash, c_par1, c_par2, c_eta,
                e_flash, e_par1, e_par2, e_eta,
                dt, dp1, dconst, eta_t, eta_m, eta_d, phi_m, psi,
                True, *dependent_params
            ]


            #substituição dos parametros do erc
            indice_escolha = {
                "Condenser : Outlet Temperature" : 7,
                "Condenser : Outlet Pressure" : 7,
                "Condenser : Outlet Vapor Fraction" : 8,
                "Condenser : Efficiency" : 9,
                "Evaporator : Outlet Temperature" : 11,
                "Evaporator : Outlet Pressure" : 11,
                "Evaporator : Outlet Vapor Fraction" : 12,
                "Evaporator : Efficiency" : 13,
                "Generator : Outlet Temperature" : 3,
                "Generator : Outlet Pressure" : 3,
                "Generator : Outlet Vapor Fraction" : 4,
                "Generator : Efficiency" : 5,
                "Pump : Isentropic Efficiency" : 1,
                "Ejector : Primary Nozzle Throat Diameter" : 14,
                "Ejector : Primary Nozzle Exit Diameter" : 15,
                "Ejector : Constant Area Section Diameter" : 16,
                "Ejector : Throat Isentropic Efficiency" : 17,
                "Ejector : Aerodynamic Throat Isentropic Efficiency" : 18,
                "Ejector : Diffuser Isentropic Efficiency" : 19,
                "Ejector : Mixing Loss Factor" : 20,
                "Ejector : Expansion Coefficient" : 21,
            }
            
            eixos_indices = []
            eixos_labels = []
            eixos_params = []
            eixos_pairs = []
            for independent_row in independent_rows_column.controls:
                equipamento = independent_row.controls[0].value
                param = independent_row.controls[1].value
                start_vu = independent_row.controls[2]
                end_vu = independent_row.controls[3]
                unit = start_vu.unit

                start_display = float(eval(start_vu.value))
                end_display = float(eval(end_vu.value))
                n = int(independent_row.controls[4].content.value)

                array_display = np.linspace(start_display, end_display, n)
                array_si = return_unit_to_SI(array_display, unit)

                eixos_indices.append(indice_escolha[f"{equipamento} : {param}"])
                eixos_labels.append(f"{equipamento}: {param}")
                eixos_params.append(param)
                eixos_pairs.append(list(zip(array_si, array_display)))

            total_estimado = 1
            for pares in eixos_pairs:
                total_estimado *= len(pares)

            batch_progress_bar.value = 0
            batch_progress_text.value = f"Preparing {total_estimado:,} scenarios..."
            nonlocal batch_cancel_requested
            batch_cancel_requested = False
            await abre_overlay_batch(e)
            await asyncio.sleep(0)

            def _gerar_cenarios():
                cenarios_ = []
                combinacoes_display_ = []
                eixos_si_colunas_ = [[] for _ in eixos_pairs]

                for combo in itertools.product(*eixos_pairs):
                    si_values = [par[0] for par in combo]
                    display_values = [par[1] for par in combo]

                    linha = list(caso_base)
                    for eixo_idx, si_valor in zip(eixos_indices, si_values):
                        linha[eixo_idx] = str(si_valor) if eixo_idx in (20, 21) else si_valor
                    cenarios_.append(linha)
                    combinacoes_display_.append(display_values)

                    for j, si_valor in enumerate(si_values):
                        eixos_si_colunas_[j].append(si_valor)
                return cenarios_, combinacoes_display_, eixos_si_colunas_

            cenarios, combinacoes_display, eixos_si_colunas = await asyncio.get_event_loop().run_in_executor(
                None, _gerar_cenarios
            )

            total = len(cenarios)
            resultados = [None] * total

            n_workers = os.cpu_count()
            CHUNK_SIZE = max(1, total // (4 * n_workers))
            chunks = [cenarios[i:i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)]
            chunk_starts = list(range(0, total, CHUNK_SIZE))
            tol, max_iter = pegar_tol_iter()
            batch_progress_text.value = f"0 / {total:,} simulations completed"
            page.update()

            if batch_parallel_checkbox.value:
                cancelled = False
                with ProcessPoolExecutor(max_workers=n_workers-1) as executor:
                    futures = [
                        executor.submit(_run_batch_chunk, chunk, backend.value, caminho_tabelas, max_iter, tol)
                        for chunk in chunks
                    ]

                    restantes = set(range(len(futures)))
                    concluidos_linhas = 0
                    while restantes:
                        if batch_cancel_requested:
                            cancelled = True
                            for i in restantes:
                                futures[i].cancel()
                            executor.shutdown(wait=False, cancel_futures=True)
                            break

                        await asyncio.sleep(0.1)
                        concluidos_agora = [i for i in restantes if futures[i].done()]
                        for i in concluidos_agora:
                            restantes.remove(i)
                            start = chunk_starts[i]
                            try:
                                chunk_result = futures[i].result()
                            except Exception as exc:
                                chunk_result = [{"error": str(exc)}] * len(chunks[i])
                            for offset, r in enumerate(chunk_result):
                                resultados[start + offset] = r
                            concluidos_linhas += len(chunk_result)

                        if concluidos_agora:
                            batch_progress_bar.value = concluidos_linhas / total
                            batch_progress_text.value = f"{concluidos_linhas:,} / {total:,} simulations completed"
                            page.update()

                if cancelled:
                    await fecha_overlay_batch(e)
                    await mostrar_erro("Batch simulation cancelled.")
                    return
            else:
                fluid = AbstractState(backend.value, dropdown.value)
                concluidos_linhas = 0
                cancelled = False
                for i, chunk in enumerate(chunks):
                    if batch_cancel_requested:
                        cancelled = True
                        break

                    chunk_result = []
                    chute = [0.0,0.0,0.0,0.0]
                    for row in chunk:
                        try:
                            rest = row[1:22]
                            dependent_params = row[23:]
                            try:
                                resultado, near_result = ERC(fluid, *rest, True, max_iter, tol, list(chute), *dependent_params)
                                chute = near_result
                            except Exception: 
                                chute = [0.0,0.0,0.0,0.0]
                                resultado, near_result = ERC(fluid, *rest, True, max_iter, tol, list(chute), *dependent_params)
                            chunk_result.append(resultado)
                        except Exception as exc:
                            chunk_result.append({"error": str(exc)})

                    start = chunk_starts[i]
                    for offset, r in enumerate(chunk_result):
                        resultados[start + offset] = r
                    concluidos_linhas += len(chunk_result)

                    batch_progress_bar.value = concluidos_linhas / total
                    batch_progress_text.value = f"{concluidos_linhas:,} / {total:,} simulations completed"
                    page.update()
                    await asyncio.sleep(0)

                if cancelled:
                    await fecha_overlay_batch(e)
                    await mostrar_erro("Batch simulation cancelled.")
                    return

            await fecha_overlay_batch(e)

            nonlocal BATCHED
            BATCHED = True

            nonlocal batch_results, batch_independent_labels, batch_independent_values, batch_dependent_params
            nonlocal batch_table_columns, batch_table_rows, current_page
            nonlocal batch_independent_si_columns, batch_independent_units

            batch_results = resultados
            batch_independent_labels = eixos_labels
            batch_independent_values = combinacoes_display
            batch_dependent_params = dependent_params
            batch_independent_si_columns = eixos_si_colunas
            batch_independent_units = [si_unit_label(p) for p in eixos_params]
            batch_table_columns, batch_table_rows = await asyncio.get_event_loop().run_in_executor(
                None, montar_linhas_tabela
            )
            current_page = 0
            atualizar_tabela_pagina(0)

            await mostrar_erro(f"Batch simulation completed: {total} runs.", icone=ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size=50))
            
            
        except Exception as erro:
            await fecha_overlay_batch(e)
            await mostrar_erro(f'{erro}')
            return

    def si_unit_label(param_name):
        if "Ratio" in param_name or param_name == "COP":
            return ""
        if "Temperature" in param_name:
            return "K"
        if "Specific Enthalpy" in param_name:
            return "J/kg"
        if "Specific Entropy" in param_name:
            return "J/kg/K"
        if "Density" in param_name:
            return "kg/m³"
        if "Mass Flow" in param_name:
            return "kg/s"
        if "Heat" in param_name or "Work" in param_name:
            return "W"
        if "Pressure" in param_name or "Backpressure" in param_name:
            return "Pa"
        if "Diameter" in param_name:
            return "m"
        return ""  # Vapor Fraction, Efficiency, Mixing Loss Factor, Expansion Coefficient, etc.

    def return_unit_to_SI(valor, unidade):
        if unidade == "°C":
            valor = valor + 273.15
        if unidade == "K":
            valor = valor
        if unidade == "°F":
            valor = (valor - 32)*5/9 + 273.15
        if unidade == "Pa":
            valor = valor
        if unidade == "kPa":
            valor = valor*1000
        if unidade == "MPa":
            valor = valor*1000*1000
        if unidade == "bar":
            valor = valor*1e5
        if unidade == "atm":
            valor = valor*101325
        if unidade == "psi":
            valor = valor*6894.7572931783
        if unidade == "mm":
            valor = valor/1000
        if unidade == "m":
            valor = valor
        if unidade == "cm":
            valor = valor/100
        if unidade == "W":
            valor = valor
        if unidade == "kW":
            valor = valor*1000
        if unidade == "MW":
            valor = valor*1000*1000
        if unidade == "hp":
            valor = valor*745.699872
        if unidade == "kg/s":
            valor = valor
        if unidade == "kg/h":
            valor = valor/3600
        if unidade == "g/s":
            valor = valor/1000
        if unidade == "lb/s":
            valor = valor*0.45359237
        if unidade == "lb/h":
            valor = valor*0.45359237/3600
        if unidade == "J/kg":
            valor = valor
        if unidade == "kJ/kg":
            valor = valor*1000
        if unidade == "MJ/kg":
            valor = valor*1000*1000
        if unidade == "Btu/lb":
            valor = valor/4.29922613938997e-4
        if unidade == "J/kg/K":
            valor = valor
        if unidade == "kJ/kg/K":
            valor = valor*1000
        if unidade == "MJ/kg/K":
            valor = valor*1000*1000
        if unidade == "kg/m³":
            valor = valor
        if unidade == "g/cm³":
            valor = valor*1000
        return valor

    def return_SI_to_unit(valor, unidade):
        if unidade == "°C":
            valor = valor - 273.15
        if unidade == "K":
            valor = valor
        if unidade == "°F":
            valor = (valor - 273.15)*9/5 + 32
        if unidade == "Pa":
            valor = valor
        if unidade == "kPa":
            valor = valor/1000
        if unidade == "MPa":
            valor = valor/1000/1000
        if unidade == "bar":
            valor = valor/1e5
        if unidade == "atm":
            valor = valor/101325
        if unidade == "psi":
            valor = valor/6894.7572931783
        if unidade == "mm":
            valor = valor*1000
        if unidade == "m":
            valor = valor
        if unidade == "cm":
            valor = valor*100
        if unidade == "W":
            valor = valor
        if unidade == "kW":
            valor = valor/1000
        if unidade == "MW":
            valor = valor/1000/1000
        if unidade == "hp":
            valor = valor/745.699872
        if unidade == "kg/s":
            valor = valor
        if unidade == "kg/h":
            valor = valor*3600
        if unidade == "g/s":
            valor = valor*1000
        if unidade == "lb/s":
            valor = valor/0.45359237
        if unidade == "lb/h":
            valor = valor/0.45359237*3600
        if unidade == "J/kg":
            valor = valor
        if unidade == "kJ/kg":
            valor = valor/1000
        if unidade == "MJ/kg":
            valor = valor/1000/1000
        if unidade == "Btu/lb":
            valor = valor*4.29922613938997e-4
        if unidade == "J/kg/K":
            valor = valor
        if unidade == "kJ/kg/K":
            valor = valor/1000
        if unidade == "MJ/kg/K":
            valor = valor/1000/1000
        if unidade == "kg/m³":
            valor = valor
        if unidade == "g/cm³":
            valor = valor/1000
        return valor

    class ValorUnidade2(ft.Row):
        def __init__(self, valor, *unidades):
            super().__init__()
            self.initial_value = valor
            self.unidades = list(unidades)
            self.unit = self.unidades[0]
            self.expand = True
            self.spacing= 0
            self.scale=ft.Scale(scale_y=1)
            self.controls = [
                ft.TextField(value= self.initial_value, expand=4, max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.only(top_left=10, bottom_left=10), border_width=1, focused_border_width=1.5, text_align=ft.TextAlign.RIGHT, hint_text="-",hint_style=ft.TextStyle(color="#677E80"), read_only=True),
                ft.Dropdown(value=self.unidades[0], options= [ft.DropdownOption(key=unidade, text=unidade, style=ft.ButtonStyle(color="#B7C3C1")) for unidade in self.unidades],focused_border_color="#B6FF57", focused_border_width=1.5, color="#B6FF57", text_style=ft.TextStyle(font_family= "Open Sans Regular", italic=False), fill_color="#354E42", filled = True, border_radius=ft.BorderRadius.only(top_right=10, bottom_right=10), border_color="#1F3835",bgcolor="#011314",
                            trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED', size=12),
                            selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED', size=12),expand = 3,
                            on_select=self.update_value)
            ]

        def update_value(self, e):
            valor_SI = return_unit_to_SI(float(eval(self.initial_value)), self.unit) # pega o valor atual e a unidade e leva pro SI
            self.unit = self.controls[1].value #agora a unidade é a nova
            self.initial_value = return_SI_to_unit(valor_SI, self.unit)
            valor_novo = f"{self.initial_value:.4f}"
            #verifica se as tres primeiras casas decimais são zero. Se for, usa notação científica
            if valor_novo.startswith("0.00") or valor_novo.startswith("-0.00"):
                valor_novo = f"{self.initial_value:.4e}" #o valor do textfield é o valor em SI convertido pra nova unidade
            #verifica se o valor é maior do que 10000 ou menor do que -10000. Se for, usa notação científica
            if float(valor_novo) > 10000 or float(valor_novo) < -10000:
                valor_novo = f"{self.initial_value:.4e}" #o valor do textfield é o valor em SI convertido pra nova unidade
            self.initial_value = str(return_SI_to_unit(valor_SI, self.unit))
            self.controls[0].value = f"{valor_novo}" #o valor do textfield é o valor em SI convertido pra nova unidade
            #isso deu mais trabalho do que eu esperava, se voce ta lendo isso, eu fiquei uns 15-30 min só tentando fazer isso
            #mas a disgrama do valor não convertia direito por causa do arredondamento, AHHHHHH

    class ValorUnidade(ft.Row):

        def __init__(self, valor, *unidades):
            super().__init__()
            self.initial_value = valor
            self.unidades = list(unidades)
            self.unit = self.unidades[0]
            self.expand = True
            self.spacing= 0
            self.scale=ft.Scale(scale_y=1)
            self.controls = [
                ft.TextField(value= self.initial_value, expand=4, max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.only(top_left=10, bottom_left=10), border_width=1, focused_border_width=1.5, text_align=ft.TextAlign.RIGHT, hint_text="-",hint_style=ft.TextStyle(color="#677E80"), on_change=self.update_valor_inicial),
                ft.Dropdown(value=self.unidades[0], options= [ft.DropdownOption(key=unidade, text=unidade, style=ft.ButtonStyle(color="#B7C3C1")) for unidade in self.unidades],focused_border_color="#B6FF57", focused_border_width=1.5, color="#B6FF57", text_style=ft.TextStyle(font_family= "Open Sans Regular", italic=False), fill_color="#354E42", filled = True, border_radius=ft.BorderRadius.only(top_right=10, bottom_right=10), border_color="#1F3835",bgcolor="#011314",
                            trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED', size=12),
                            selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED', size=12),expand = 3)
            ]

        def update_valor_inicial(self, e):
            self.initial_value = self.controls[0].value
            self.unit = self.controls[1].value

        #acabou q eu nem usei esse método né, mas fazer o que, itis de laifi
        def return_to_SI(self): 
            if self.controls[1].value == "°C":
                self.value = float(eval(self.controls[0].value) +273.15)
            if self.controls[1].value == "K":
                self.value = float(eval(self.controls[0].value))
            if self.controls[1].value == "°F":
                self.value = float((eval(self.controls[0].value)-32)*5/9 +273.15)
            if self.controls[1].value == "Pa":
                self.value = float(eval(self.controls[0].value))
            if self.controls[1].value == "kPa":
                self.value = float(eval(self.controls[0].value)*1000)
            if self.controls[1].value == "MPa":
                self.value = float(eval(self.controls[0].value)*1000*1000)
            if self.controls[1].value == "bar":
                self.value = float(eval(self.controls[0].value)*1e5)
            if self.controls[1].value == "atm":
                self.value = float(eval(self.controls[0].value)*101325)
            if self.controls[1].value == "psi":
                self.value = float(eval(self.controls[0].value)*6894.7572931783)

    class ValorUnidade3(ft.Row):
        def __init__(self, valor, hint, *unidades, param_name=None):
            super().__init__()
            self.initial_value = valor
            self.hint = hint
            self.param_name = param_name or ""

            # Flags de validação
            self.is_temperature = "Temperature" in self.param_name
            self.is_pressure = "Pressure" in self.param_name
            self.is_vapor_fraction = "Vapor Fraction" in self.param_name
            self.is_efficiency = "Efficiency" in self.param_name
            self.is_dimension = "Diameter" in self.param_name
            self.is_unit_ratio = self.param_name in ("Mixing Loss Factor", "Expansion Coefficient")

            self.text_field = ft.TextField(
                value=self.initial_value,
                expand=12 if unidades else 4,
                max_lines=1, multiline=False,
                text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),
                bgcolor="#172B2B",
                focused_border_color="#B6FF57",
                border_color="#1F3835",
                border_radius=ft.BorderRadius.only(top_left=0, bottom_left=0),
                border_width=1,
                focused_border_width=1.5,
                text_align=ft.TextAlign.RIGHT,
                hint_text=self.hint,
                hint_style=ft.TextStyle(color="#677E80", italic=True, size=14, letter_spacing=1.5),
                on_blur=self.validate_value,
            )

            if unidades:
                self.unidades = list(unidades)
                self.unit = self.unidades[0]
                self.expand = True
                self.spacing = 0
                self.scale = ft.Scale(scale_y=1)
                self.unit_dropdown = ft.Dropdown(
                    value=self.unidades[0],
                    options=[ft.DropdownOption(key=u, text=u, style=ft.ButtonStyle(color="#B7C3C1")) for u in self.unidades],
                    focused_border_color="#B6FF57",
                    focused_border_width=1.5,
                    color="#B6FF57",
                    text_style=ft.TextStyle(font_family="Open Sans Regular", italic=False),
                    fill_color="#354E42",
                    filled=True,
                    border_radius=ft.BorderRadius.only(top_right=0, bottom_right=0),
                    border_color="#1F3835",
                    bgcolor="#011314",
                    trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED', size=12),
                    selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED', size=12),
                    expand=11,
                    on_select=self.on_unit_change,
                )
                self.controls = [self.text_field, self.unit_dropdown]
            else:
                self.unidades = []
                self.unit = None
                self.expand = True
                self.spacing = 0
                self.scale = ft.Scale(scale_y=1)
                self.controls = [self.text_field]

        def on_unit_change(self, e):
            # Troca de unidade pode invalidar um valor já digitado (ex: negativo virando Kelvin)
            self.unit = self.unit_dropdown.value
            self.validate_value(None)

        def validate_value(self, e):
            text = self.text_field.value

            if not text:
                self.text_field.error_text = None
                self.text_field.update()
                return

            try:
                value = float(text)
            except ValueError:
                self.text_field.value = ""
                self.text_field.error_text = "Must be a number"
                self.text_field.update()
                return

            if self.is_temperature and self.unit == "K" and value < 0:
                self.text_field.value = ""
                self.text_field.error_text = "Kelvin cannot be negative"
                self.text_field.update()
                return

            if self.is_pressure and value < 0:
                self.text_field.value = ""
                self.text_field.error_text = "Pressure cannot be negative"
                self.text_field.update()
                return

            if self.is_vapor_fraction and not (0 <= value <= 1):
                self.text_field.value = ""
                self.text_field.error_text = "Must be between 0 and 1"
                self.text_field.update()
                return

            if self.is_efficiency and value <= 0:
                self.text_field.value = ""
                self.text_field.error_text = "Cannot be negative"
                self.text_field.update()
                return

            if self.is_dimension and value <= 0:
                self.text_field.value = ""
                self.text_field.error_text = "Dimension cannot be negative"
                self.text_field.update()
                return

            if self.is_unit_ratio and not (0 < value <= 1):
                self.text_field.value = ""
                self.text_field.error_text = "Must be strictly between 0 and 1"
                self.text_field.update()
                return

            self.text_field.error_text = None
            self.text_field.update()

        @property
        def value(self):
            return self.text_field.value

    def gerar_tabela(fluido):
        #eu me sinto meio burro de ter criado uma função só pra isso, mas agora to com preguiça de mudar
        return AbstractState("BICUBIC&HEOS", fluido) 
    
    def tabela_existe(fluido: str) -> bool:
        try:
            arquivos = os.listdir(caminho_tabelas)
            return any(fluido in arquivo for arquivo in arquivos)
        except FileNotFoundError:
            return False
        
    async def mostrar_erro(mensagem: str, icone: ft.Icon = ft.Icon(icon=ft.Icons.ERROR_OUTLINE, size = 50, color = "#FF5C5C")):
        overlay_erro.content.content.controls[0] = icone
        overlay_erro.content.content.controls[1].value = mensagem
        overlay_erro.visible = True
        overlay_erro.opacity = 0
        page.update()
        await asyncio.sleep(0.05) #parece que fica meio bugado sem isso
        overlay_erro.opacity = 1
        page.update()

    async def fechar_erro(e):
        overlay_erro.opacity = 0
        page.update()
        await asyncio.sleep(0.3)
        overlay_erro.visible = False
        page.update()

    overlay_erro = ft.Container(
        bgcolor="#50000000",
        on_click=fechar_erro,
        expand=True,
        visible=False,
        blur=10,
        opacity=0,
        alignment=ft.Alignment.CENTER,
        animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT),
        content=ft.Container(
            alignment=ft.Alignment.CENTER,
            height=160,
            width=320,
            bgcolor="#0A2022",
            border_radius=10,
            border=ft.Border.all(1, "#50B7C3C1"),
            content=ft.Column(
                alignment=ft.MainAxisAlignment.CENTER,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=12,
                controls=[
                    ft.Icon(ft.Icons.ERROR_OUTLINE, color="#FF5C5C", size=32),
                    ft.Text(
                        "mensagem aqui",   # atualizado dinamicamente
                        style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5),
                        text_align=ft.TextAlign.CENTER,
                    ),
                    ft.Text(
                        "(click anywhere to close)",
                        style=ft.TextStyle(color="#677E80", size=12, letter_spacing=1.5),
                    ),
                ],
            ),
        ),
    )

    async def abre_overlay(e):
        overlay_bake_fluid.visible = True
        overlay_bake_fluid.opacity = 0
        page.update()
        await asyncio.sleep(0.05)
        overlay_bake_fluid.opacity = 1
        page.update()

    async def fecha_overlay(e):
        overlay_bake_fluid.opacity = 0
        page.update()
        await asyncio.sleep(0.3)
        overlay_bake_fluid.visible=False
        page.update()

    async def bake_fluid(e, run=False):
        if run:
            if not tabela_existe(dropdown.value):
                await abre_overlay(e)
                await asyncio.sleep(0.05)
                await asyncio.get_event_loop().run_in_executor(None, lambda: gerar_tabela(dropdown.value))
                await asyncio.sleep(0.05)
                await fecha_overlay(e)
        else:
            try:
                AbstractState("HEOS", dropdown.value)
            except Exception:
                await mostrar_erro('Please select a valid fluid')
                return
            if tabela_existe(dropdown.value):
                await mostrar_erro('Fluid already baked!', icone= ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size = 50))
            else:
                await abre_overlay(e)
                await asyncio.sleep(0.05)
                await asyncio.get_event_loop().run_in_executor(None, lambda: gerar_tabela(dropdown.value))
                await asyncio.sleep(0.05)
                await fecha_overlay(e)

    def animar_janela(e):
        nonlocal painel_aberto
        if painel_aberto:
            # Desloca o painel esquerdo totalmente para fora da tela (em pixels)
            painel_esquerdo.left = -(LARGURA_PAINEL + 20)
            # Zera a margem do direito para ele expandir e ocupar a tela toda
            painel_direito.margin = ft.Margin.only(left=0)
        else:
            # Traz o painel esquerdo de volta para a tela
            painel_esquerdo.left = 0
            # Restaura a margem para empurrar o painel direito
            painel_direito.margin = ft.Margin.only(left=LARGURA_PAINEL + SPACING)
        
        painel_aberto = not painel_aberto
        page.update()

    last_clicked = ''
    def gerenciar_clique(e):
        nonlocal last_clicked
        
        # Mudar os values pelo conteudo em si
        dicionario_nomes = {
            'pump': bomba_painel_esquerdo,
            'generator': gerador_painel_esquerdo,
            'condenser': condensador_painel_esquerdo,
            'ejector': ejetor_painel_esquerdo,
            'evaporator': evaporador_painel_esquerdo,
            'expansion_valve': valvula_de_expansao_painel_esquerdo
        }

        # Se clicou FORA do seletor
        if e.control.data == 'out_click':
            # Apaga todo mundo
            for i in range(1, 7):
                erc.controls[i].gradient.colors = ['#00B3F364', '#00FFFFFF']
                erc.controls[i].border = ft.Border.all(0.5, "#00B6FF57")
            
            last_clicked = ''
            
            # Fecha o painel se ele estiver aberto
            if painel_aberto == True:
                animar_janela(e)
            
            e.page.update()
            return 

        #Se clicou em um componente E o painel está fechado, abre o painel
        if painel_aberto == False:
            animar_janela(e)

        # Atualiza o texto do painel
        if e.control.data in dicionario_nomes:
            animador_de_telas.content = dicionario_nomes[e.control.data]

        #RESET GERAL: Apaga o brilho de todos os componentes primeiro
        for i in range(1, 7):
            erc.controls[i].gradient.colors = ['#00B3F364', '#00FFFFFF']
            erc.controls[i].border = ft.Border.all(0.5, "#00B6FF57")

        #ACENDE APENAS O CLICADO
        e.control.gradient.colors = ['#80B6FF57', '#00FFFFFF']
        e.control.border = ft.Border.all(0.5, "#B6FF57")
        
        last_clicked = e.control.data

        e.page.update()

#=================== INICIO:  Botões Superiores ==============
    dropdown = ft.Dropdown(
        expand=True,
        hint_text='Fluid',
        hint_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),
        options=fluid_list(),
        color="#E3EDED",
        bgcolor="#011314",
        menu_height=400,
        fill_color="#172B2B",
        filled=True,   
        text_style=ft.TextStyle(font_family="Open Sans Light",size=14, letter_spacing=1.5),
        dense=True,
        editable=True,
        enable_filter=True,
        enable_search=True,
        border_radius=ft.BorderRadius(top_left=10, top_right=0, bottom_left=0, bottom_right=0),
        trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
        selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
        width=150,
        focused_border_color="#B6FF57",
        focused_border_width=1.5,
        border_color="#1F3835",
        tooltip = tooltip_msg("Select the working fluid for the ERC"),
    )

    bake_button = ft.OutlinedButton(
        expand = True,
        tooltip=tooltip_msg(f"Create the bicubic interpolation tables for the fluid and save them in \n{caminho_tabelas}\nYou only need to bake the fluid once. You also need to choose the BICUBIC backend to use it.", 800),
        content=ft.Text("Bake Fluid", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5)),
        height=48,   # mesma altura do dropdown
        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(radius=0),
        ),
        on_click=bake_fluid
    )

    backend = ft.Dropdown(
        expand=True,
        hint_text='Backend',
        hint_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),
        options=[ft.DropdownOption(key='HEOS', text='HEOS', style=ft.ButtonStyle(color="#B7C3C1")), ft.DropdownOption(key='BICUBIC&HEOS', text='BICUBIC', style=ft.ButtonStyle(color="#B7C3C1"))],
        color="#E3EDED",
        bgcolor="#011314",
        fill_color="#172B2B",
        filled=True,   
        text_style=ft.TextStyle(font_family="Open Sans Light",size=14, letter_spacing=1.5),
        dense=True,
        border_radius=0,
        trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
        selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
        focused_border_color="#B6FF57",
        focused_border_width=1.5,
        border_color="#1F3835",
        tooltip = tooltip_msg("Select the backend:\n\nHEOS (Helmholtz Equation of State) is the standard and most robust backend of the CoolProp.\n\nBICUBIC uses pre-calculated tables (based on HEOS) in a regular grid of states and uses bicubic interpolation.\n             It is faster than HEOS, but requires that the tables for each fluid be generated once (Bake Fluid). Useful when there are many simulations."),
    )

    run_simulation_button = ft.IconButton(
        expand = True,
        tooltip=tooltip_msg(f"RUN SIMULATION, RUN!", 100),
        icon=ft.Icon(icon= ft.Icons.PLAY_ARROW_ROUNDED, color="#B6FF57"),
        height=48,   # mesma altura do dropdown
        width=48,

        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(
                radius=ft.BorderRadius(top_left=0, bottom_left=0, top_right=0, bottom_right=0)
            ),
            
        ),
        on_click=run
    )
    
    batch_sim_button = ft.OutlinedButton(
        expand = True,
        tooltip=tooltip_msg(f"Run batch simulations or parametric studies", 800),
        content=ft.Text("Batch Sim", style=ft.TextStyle(font_family="Open Sans Light", size=14, letter_spacing=1.5)),
        height=48,   # mesma altura do dropdown
        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(radius=0),
        ),
        on_click= batch_sim
    )

    results_button = ft.OutlinedButton(
        expand = True,
        tooltip=tooltip_msg(f"View batch simulation results", 800),
        content=ft.Text("Batch Sim Table", style=ft.TextStyle(font_family="Open Sans Light", size=14, letter_spacing=1.5)),
        height=48,   # mesma altura do dropdown
        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(radius=0),
        ),
        on_click= results
    )

    settings_button = ft.IconButton(
        expand = False,
        tooltip=tooltip_msg(f"Simulation Settings", 100),
        icon=ft.Icon(icon= ft.Icons.SETTINGS, color="#B7C3C1"),
        height=48,   # mesma altura do dropdown
        width=48,

        style=ft.ButtonStyle(
            side={
                ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#1F3835"),
                ft.ControlState.HOVERED: ft.BorderSide(1.5, "#B6FF57"),
            },
            color={
                ft.ControlState.DEFAULT: "#B7C3C1",
                ft.ControlState.HOVERED: "#B6FF57",
            },
            bgcolor="#172B2B",
            shape=ft.RoundedRectangleBorder(
                radius=ft.BorderRadius(top_left=0, bottom_left=0, top_right=10, bottom_right=0)
            ),
            
        ),
        on_click=abrir_settings
    )
#=================== FIM:  Botões Superiores =================
    
    def texto_com_destaque(texto, titulo):
        texto = ft.Column(controls=[ft.Container(content = ft.Row(
            spacing=12,
            alignment=ft.MainAxisAlignment.START,
            controls=[
                ft.Container(
                    width=4,
                    height=46,
                    bgcolor="#B6FF57", 
                    border_radius=2,
                )
                ,
                ft.Column(
                    spacing=2,
                    controls=[
                        ft.Text(titulo, size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular')),
                        ft.Text(texto, size=22, color="#E3EDED", style=ft.TextStyle(letter_spacing=2))
                    ]
                )
            ]
        ))])
        return texto

    def nome_equipamento(nome):
        texto = ft.Row(
            spacing=12,
            alignment=ft.MainAxisAlignment.START,
            controls=[
                ft.Container(
                    width=4,
                    height=46,
                    bgcolor="#B6FF57", 
                    border_radius=2,
                )
                ,
                ft.Column(
                    spacing=2,
                    controls=[
                        ft.Text("SYSTEM COMPONENT", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular')),
                        ft.Text(nome, size=22, color="#E3EDED", style=ft.TextStyle(letter_spacing=2))
                    ]
                )
            ]
        )
        return texto

    #======= BOMBA PAINEL ESQUERDO ==============
    pump_description_text = ft.Text(spans=[ft.TextSpan('A device that takes the liquid coming out of the condenser to a high pressure. ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://github.com/GPCTM-BR/SIMERC/blob/f70d33508b99be33ee1ecaa04b4c895f418e602d/docstrings/images/System%20Components%20Cards.pdf", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    pump_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        pump_description_text
    ])
    pump_column_controls = [
        pump_description,
        ft.Text("CALCULATION PARAMETERS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand = True, content=
            ft.Row(
                expand = True,
                controls=[
                    ft.Container(expand = True, content=
                        ft.Row(controls= [
                            ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                            ft.Text("Isentropic Efficiency (0-100%)", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),    
                            ],
                            
                        )
                    ),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "75",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
                ]
            )
        ),
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    #Para adicionar mais coisas, adicionar no pump_column_controls
    #Depois de simular, trocarr o content do container (pump_column_controls[4].content = os resultados (ou só pump_column_controls[4], sei lá))
    bomba_painel_esquerdo = ft.Column(
        key='tela_bomba',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(expand=3, content=ft.Image("/images/pump.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER),
            nome_equipamento("PUMP"),
            ft.Column(controls=pump_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )
    
    #======= GERADOR PAINEL ESQUERDO ==============
    def update_generator_parameters(e):
        generator_column_controls[2].controls[0].controls[1] = texto_na_esquerda(heater_cooler_flash[g_flash_type_dropdown.value][0], size = 14,tooltip_str=heater_cooler_flash[g_flash_type_dropdown.value][2])
        generator_column_controls[2].controls[0].controls[2] = texto_na_esquerda(heater_cooler_flash[g_flash_type_dropdown.value][1], size = 14,tooltip_str=heater_cooler_flash[g_flash_type_dropdown.value][3])
        if g_flash_type_dropdown.value == "TP":
            generator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            generator_column_controls[2].controls[1].controls[2] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
        if g_flash_type_dropdown.value == "TQ":
            generator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            generator_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
        if g_flash_type_dropdown.value == "PQ":
            generator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
            generator_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)

    def texto_na_esquerda(texto, size=14, latex_symbol = None, tooltip_str = None, tooltip_time = 300, margin = 0):
        return ft.Container(expand = True, margin = margin, content=
                        ft.Row(controls= [
                            ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                            ft.Text(texto, size=size, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),    
                            ft.Markdown(value= latex_symbol, latex_style=ft.TextStyle(font_family="Open Sans Light", size=20, color="#B7C3C1"))    
                            ],
                            
                        ), tooltip=tooltip_msg(tooltip_str,tooltip_time)
                    )
    
    heater_cooler_flash ={
        "TP": ["Outlet Temperature", "Outlet Pressure","Define the outlet temperature", "Define the outlet pressure"],
        "TQ": ["Outlet Temperature", "Outlet Vapor Fraction", "Define the outlet temperature","The vapor fraction must be a number between 0 and 1"],
        "PQ": ["Outlet Pressure", "Outlet Vapor Fraction","Define the outlet pressure", "The vapor fraction must be a number between 0 and 1"],
    }
    g_flash_type_dropdown = ft.Dropdown(
        value="TP",
        options=[ft.DropdownOption(key='TP', text='Outlet Temperature and Pressure', style=ft.ButtonStyle(color="#B7C3C1")), 
                 ft.DropdownOption(key='TQ', text='Outlet Temperature and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1")),
                 ft.DropdownOption(key='PQ', text='Outlet Pressure and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1"))],
        color="#E3EDED",
        bgcolor="#011314",
        fill_color="#354E42",
        filled=True,   
        text_style=ft.TextStyle(font_family="Open Sans Light",size=10),
        dense=True,
        border_radius=10,
        trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
        selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
        focused_border_color="#B6FF57",
        focused_border_width=1.5,
        border_color="#1F3835",
        tooltip = tooltip_msg("Select the calculation method"),
        on_select= update_generator_parameters
    )
    generator_description_text = ft.Text(spans=[ft.TextSpan('A heater that produces high-pressure vapor that acts as a motive flow for the ejector. ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://github.com/GPCTM-BR/SIMERC/blob/f70d33508b99be33ee1ecaa04b4c895f418e602d/docstrings/images/System%20Components%20Cards.pdf", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    generator_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        generator_description_text
    ])
    generator_column_controls =[
        generator_description,
        ft.Text("CALCULATION PARAMETERS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Row(
            expand = True,
            controls =[
                ft.Column(spacing=25, expand = True,controls=[
                    texto_na_esquerda("Flash Type", size = 14),
                    texto_na_esquerda(heater_cooler_flash[g_flash_type_dropdown.value][0], size = 14, tooltip_str=heater_cooler_flash[g_flash_type_dropdown.value][2]),
                    texto_na_esquerda(heater_cooler_flash[g_flash_type_dropdown.value][1], size = 14, tooltip_str=heater_cooler_flash[g_flash_type_dropdown.value][3]),
                    ft.Container(expand = True, content=
                        ft.Row(controls= [
                            ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                            ft.Text("Efficiency (0-100%)", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),    
                            ],
                            
                        )
                    ),
                    
                ]),
                ft.Column(expand = True,controls=[
                    g_flash_type_dropdown,
                    ValorUnidade("", "°C","K","°F"),
                    ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi"),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "100",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
                ])
            ]
        ),
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    gerador_painel_esquerdo = ft.Column(
        key = 'tela_gerador',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(expand=3, content=ft.Image("/images/generator.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER),
            nome_equipamento("GENERATOR"),
            ft.Column(controls=generator_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )

    #======= EJETOR PAINEL ESQUERDO ==============
    ejector_description_text = ft.Text(spans=[ft.TextSpan('A CPM ejector model. It uses the modeling developed by Cardemil and Colle (2012). ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://doi.org/10.1016/j.enconman.2012.05.009", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    ejector_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ejector_description_text
    ])
    ejector_column_controls =[
        ejector_description,
        ft.Text("CALCULATION PARAMETERS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Row(
            expand = True,
            controls = [
                #PRIMEIRA COLUNA
                ft.Column(
                    spacing=16, expand = True,
                    controls=[
                        texto_na_esquerda("Primary Nozzle Throat Diameter", latex_symbol=r"$d_t$"),
                        texto_na_esquerda("Primary Nozzle Exit Diameter", latex_symbol=r"$d_{p1}$"),
                        texto_na_esquerda("Constant Area Section Diameter", latex_symbol=r"$d_{const}$"),
                        texto_na_esquerda("Throat Isentropic Efficiency (0-1)", latex_symbol=r"$\eta_t$"),
                        texto_na_esquerda("Aerodynamic Throat Isentropic Efficiency (0-1)", latex_symbol=r"$\eta_m$"),
                        texto_na_esquerda("Diffuser Isentropic Efficiency (0-1)", latex_symbol=r"$\eta_d$"),
                        texto_na_esquerda("Mixing Loss Factor", latex_symbol=r"$\phi_m$", tooltip_str="It can be a fixed value or an expression in terms of Ar or Pr (or both).", tooltip_time=100),
                        texto_na_esquerda("Expansion Coefficient", latex_symbol=r"$\psi$",tooltip_str="It can be a fixed value or an expression in terms of Ar or Pr (or both).", tooltip_time=100)
                    ]
                ),
                #SEGUNDA COLUNA
                ft.Column(
                    expand = True,
                    controls=[
                        ValorUnidade("","mm","m","cm"),
                        ValorUnidade("","mm","m","cm"),
                        ValorUnidade("","mm","m","cm"),
                        ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "0.95",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5)),
                        ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "0.95",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5)),
                        ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "0.95",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5)),
                        ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5)),
                        ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5)),
                    ]

                )
            ]
        ),
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    ejetor_painel_esquerdo = ft.Column(
        key = 'tela_ejetor',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(expand=3, content=ft.Image("/images/ejector.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER),
            nome_equipamento("EJECTOR"),
            ft.Column(controls=ejector_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )

    #======= CONDENSADOR PAINEL ESQUERDO ==============
    def update_condenser_parameters(e):
        condenser_column_controls[2].controls[0].controls[1] = texto_na_esquerda(heater_cooler_flash[c_flash_type_dropdown.value][0], size = 14,tooltip_str=heater_cooler_flash[c_flash_type_dropdown.value][2])
        condenser_column_controls[2].controls[0].controls[2] = texto_na_esquerda(heater_cooler_flash[c_flash_type_dropdown.value][1], size = 14,tooltip_str=heater_cooler_flash[c_flash_type_dropdown.value][3])
        if c_flash_type_dropdown.value == "TP":
            condenser_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            condenser_column_controls[2].controls[1].controls[2] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
        if c_flash_type_dropdown.value == "TQ":
            condenser_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            condenser_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
        if c_flash_type_dropdown.value == "PQ":
            condenser_column_controls[2].controls[1].controls[1] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
            condenser_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
    
    c_flash_type_dropdown = ft.Dropdown(
        value="TP",
        options=[ft.DropdownOption(key='TP', text='Outlet Temperature and Pressure', style=ft.ButtonStyle(color="#B7C3C1")), 
                 ft.DropdownOption(key='TQ', text='Outlet Temperature and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1")),
                 ft.DropdownOption(key='PQ', text='Outlet Pressure and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1"))],
        color="#E3EDED",
        bgcolor="#011314",
        fill_color="#354E42",
        filled=True,   
        text_style=ft.TextStyle(font_family="Open Sans Light",size=10),
        dense=True,
        border_radius=10,
        trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
        selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
        focused_border_color="#B6FF57",
        focused_border_width=1.5,
        border_color="#1F3835",
        tooltip = tooltip_msg("Select the calculation method"),
        on_select= update_condenser_parameters
    )
    condenser_description_text = ft.Text(spans=[ft.TextSpan('A cooler that works by condensing the fluid that comes out of the ejector. ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://github.com/GPCTM-BR/SIMERC/blob/f70d33508b99be33ee1ecaa04b4c895f418e602d/docstrings/images/System%20Components%20Cards.pdf", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    condenser_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        condenser_description_text
    ])
    condenser_column_controls =[
        condenser_description,
        ft.Text("CALCULATION PARAMETERS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Row(
            expand = True,
            controls =[
                ft.Column(spacing=25, expand = True,controls=[
                    texto_na_esquerda("Flash Type", size = 14),
                    texto_na_esquerda(heater_cooler_flash[c_flash_type_dropdown.value][0], size = 14,tooltip_str=heater_cooler_flash[c_flash_type_dropdown.value][2]),
                    texto_na_esquerda(heater_cooler_flash[c_flash_type_dropdown.value][1], size = 14,tooltip_str=heater_cooler_flash[c_flash_type_dropdown.value][3]),
                    ft.Container(expand = True, content=
                        ft.Row(controls= [
                            ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                            ft.Text("Efficiency (0-100%)", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),    
                            ],
                            
                        )
                    ),
                    
                ]),
                ft.Column(expand = True,controls=[
                    c_flash_type_dropdown,
                    ValorUnidade("", "°C","K","°F"),
                    ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi"),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "100",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
                ])
            ]
        ),
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    condensador_painel_esquerdo = ft.Column(
        key = 'tela_condensador',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(expand=3, content=ft.Image("/images/condenser.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER),
            nome_equipamento("CONDENSER"),
            ft.Column(controls=condenser_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )

    #======= EVAPORADOR PAINEL ESQUERDO ==============
    def update_evaporator_parameters(e):
        evaporator_column_controls[2].controls[0].controls[1] = texto_na_esquerda(heater_cooler_flash[e_flash_type_dropdown.value][0], size = 14,tooltip_str=heater_cooler_flash[e_flash_type_dropdown.value][2])
        evaporator_column_controls[2].controls[0].controls[2] = texto_na_esquerda(heater_cooler_flash[e_flash_type_dropdown.value][1], size = 14,tooltip_str=heater_cooler_flash[e_flash_type_dropdown.value][3])
        if e_flash_type_dropdown.value == "TP":
            evaporator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            evaporator_column_controls[2].controls[1].controls[2] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
        if e_flash_type_dropdown.value == "TQ":
            evaporator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "°C", "K", "°F")
            evaporator_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
        if e_flash_type_dropdown.value == "PQ":
            evaporator_column_controls[2].controls[1].controls[1] = ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi")
            evaporator_column_controls[2].controls[1].controls[2] = ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
    
    e_flash_type_dropdown = ft.Dropdown(
        value="TP",
        options=[ft.DropdownOption(key='TP', text='Outlet Temperature and Pressure', style=ft.ButtonStyle(color="#B7C3C1")), 
                 ft.DropdownOption(key='TQ', text='Outlet Temperature and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1")),
                 ft.DropdownOption(key='PQ', text='Outlet Pressure and Vapor Fraction', style=ft.ButtonStyle(color="#B7C3C1"))],
        color="#E3EDED",
        bgcolor="#011314",
        fill_color="#354E42",
        filled=True,   
        text_style=ft.TextStyle(font_family="Open Sans Light",size=10),
        dense=True,
        border_radius=10,
        trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
        selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
        focused_border_color="#B6FF57",
        focused_border_width=1.5,
        border_color="#1F3835",
        tooltip = tooltip_msg("Select the calculation method"),
        on_select= update_evaporator_parameters
    )
    evaporator_description_text = ft.Text(spans=[ft.TextSpan('A heater that produces low-pressure vapor, which acts as entrainment flow for the ejector. ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://github.com/GPCTM-BR/SIMERC/blob/f70d33508b99be33ee1ecaa04b4c895f418e602d/docstrings/images/System%20Components%20Cards.pdf", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    evaporator_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        evaporator_description_text
    ])
    evaporator_column_controls =[
        evaporator_description,
        ft.Text("CALCULATION PARAMETERS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Row(
            expand = True,
            controls =[
                ft.Column(spacing=25, expand = True,controls=[
                    texto_na_esquerda("Flash Type", size = 14),
                    texto_na_esquerda(heater_cooler_flash[e_flash_type_dropdown.value][0], size = 14,tooltip_str=heater_cooler_flash[e_flash_type_dropdown.value][2]),
                    texto_na_esquerda(heater_cooler_flash[e_flash_type_dropdown.value][1], size = 14,tooltip_str=heater_cooler_flash[e_flash_type_dropdown.value][3]),
                    ft.Container(expand = True, content=
                        ft.Row(controls= [
                            ft.Container(height = 35,width=4,bgcolor="#3C544B", border_radius=2),
                            ft.Text("Efficiency (0-100%)", size=14, color="#B7C3C1", style=ft.TextStyle(font_family='Open Sans Light', letter_spacing=1.5), max_lines=10, expand = True),    
                            ],
                            
                        )
                    ),
                    
                ]),
                ft.Column(expand = True,controls=[
                    e_flash_type_dropdown,
                    ValorUnidade("", "°C","K","°F"),
                    ValorUnidade("", "kPa", "Pa", "MPa", "atm", "bar", "psi"),
                    ft.Container(expand=True, content=ft.TextField(text_align=ft.TextAlign.RIGHT,value= "100",hint_text="-",hint_style=ft.TextStyle(color="#677E80"), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.all(10), border_width=1, focused_border_width=1.5),)
                ])
            ]
        ),
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    evaporador_painel_esquerdo = ft.Column(
        key = 'tela_evaporador',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(content=ft.Image("/images/evaporator.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER,expand=3),
            nome_equipamento("EVAPORATOR"),
            ft.Column(controls=evaporator_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )

    #======= VÁLVULA PAINEL ESQUERDO ==============
    expansion_valve_description_text = ft.Text(spans=[ft.TextSpan('A throttling device that reduces the pressure and temperature of the working fluid. ',style=ft.TextStyle(italic=True, letter_spacing=1.5, size = 12, color = "#FFFFFF")),ft.TextSpan("[See here]",url ="https://github.com/GPCTM-BR/SIMERC/blob/f70d33508b99be33ee1ecaa04b4c895f418e602d/docstrings/images/System%20Components%20Cards.pdf", style=ft.TextStyle(size=12, italic=False, color="#B6FF57"))])
    expansion_valve_description = ft.Column(controls=[
        ft.Text("DESCRIPTION", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        expansion_valve_description_text
    ])
    expansion_valve_column_controls = [
        expansion_valve_description,
        ft.Text("SIMULATION RESULTS", size=10, color="#677E80", style=ft.TextStyle(font_family='Open Sans Regular'), margin=ft.Margin.only(top=15)),
        ft.Container(expand=True, content=
                     ft.Text("There are no simulation results yet.", style=ft.TextStyle(font_family="Opens Sans Light", letter_spacing=1.5, color="#B7C3C1", italic=True))
                     )
    ]
    #trocar o expansion_valve_column_controls[2]
    valvula_de_expansao_painel_esquerdo = ft.Column(
        key = 'tela_valvula_de_expansao',
        horizontal_alignment= ft.CrossAxisAlignment.START,
        controls=[
            ft.Container(content=ft.Image("/images/expansion_valve.webp", fit=ft.BoxFit.FIT_HEIGHT), padding=10, border_radius=10, gradient=ft.LinearGradient(colors=["#375043","#041A1C"], begin=ft.Alignment.BOTTOM_LEFT, end=ft.Alignment.TOP_RIGHT, stops=[0,0.5]), border=ft.Border.all(1,color="#504B5F5F"), alignment=ft.Alignment.CENTER, expand = 3),
            nome_equipamento("EXPANSION VALVE"),
            ft.Column(controls=expansion_valve_column_controls, scroll=ft.ScrollMode.AUTO, expand=5) #posso adicionar coisa aqui depois, e só vai scrollar esse conteúdo
        ]
    )
    #=============================================

    #======== Conteúdo das outras abas ===============
    def back_to_fluxogram(e):
        coluna.controls[1].content = ft.Container(
            key ="erc",
            content=erc,
            alignment=ft.Alignment.CENTER,
            expand=True,
            data='out_click',
            on_click=gerenciar_clique,
            image= ft.DecorationImage(src="/images/pattern3.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=20, opacity=0.5)
        )

    independent_rows_column = ft.Column(spacing=10, controls=[])
    dependent_rows_column = ft.Column(spacing=10, controls=[])

    def create_dependent_row():
        
        row_id = str(uuid.uuid4())
        def remove_row(e):
            for r in dependent_rows_column.controls:
                if r.data == row_id:
                    dependent_rows_column.controls.remove(r)
                    break
            dependent_rows_column.update()
        
        async def validate(e):
            current_param = e.control.value
            e.control.value = current_param
            e.control.update()
            if not current_param:
                return

            if current_param not in param_list:
                await mostrar_erro(f"'{current_param}' is not a valid dependent variable.")
                e.control.value = None
                e.control.update()  # por alguma razão não ta funcionando
                return

            for r in dependent_rows_column.controls:
                row_param = r.controls[0].value
                if r.data != row_id and row_param == current_param:
                    await mostrar_erro(f"'{current_param}' has already been added.")
                    e.control.value = None
                    e.control.update()  
                    return

        param_list =["COP",
                     "EJECTOR: Entrainment Ratio", "EJECTOR: Pressure Lift Ratio", "EJECTOR: Critical Backpressure", "EJECTOR: Primary Inlet Mass Flow",
                     "EJECTOR: Secondary Inlet Mass Flow", "EJECTOR: Outlet Mass Flow", "EJECTOR: Primary Inlet Temperature", "EJECTOR: Primary Inlet Pressure",
                     "EJECTOR: Primary Inlet Specific Enthalpy", "EJECTOR: Primary Inlet Specific Entropy", "EJECTOR: Primary Inlet Density",
                     "EJECTOR: Secondary Inlet Temperature", "EJECTOR: Secondary Inlet Pressure",
                     "EJECTOR: Secondary Inlet Specific Enthalpy", "EJECTOR: Secondary Inlet Specific Entropy", "EJECTOR: Secondary Inlet Density",
                     "EJECTOR: Outlet Temperature", "EJECTOR: Outlet Pressure",
                     "EJECTOR: Outlet Specific Enthalpy", "EJECTOR: Outlet Specific Entropy", "EJECTOR: Outlet Density",
                     "EVAPORATOR: Heaty Duty", "EVAPORATOR: Mass Flow", 
                     "EVAPORATOR: Inlet Temperature", "EVAPORATOR: Inlet Pressure",
                     "EVAPORATOR: Inlet Specific Enthalpy", "EVAPORATOR: Inlet Specific Entropy", "EVAPORATOR: Inlet Density",
                     "EVAPORATOR: Outlet Temperature", "EVAPORATOR: Outlet Pressure",
                     "EVAPORATOR: Outlet Specific Enthalpy", "EVAPORATOR: Outlet Specific Entropy", "EVAPORATOR: Outlet Density",
                     "GENERATOR: Heaty Duty", "GENERATOR: Mass Flow", 
                     "GENERATOR: Inlet Temperature", "GENERATOR: Inlet Pressure",
                     "GENERATOR: Inlet Specific Enthalpy", "GENERATOR: Inlet Specific Entropy", "GENERATOR: Inlet Density",
                     "GENERATOR: Outlet Temperature", "GENERATOR: Outlet Pressure",
                     "GENERATOR: Outlet Specific Enthalpy", "GENERATOR: Outlet Specific Entropy", "GENERATOR: Outlet Density",
                     "CONDENSER: Heaty Duty", "CONDENSER: Mass Flow", 
                     "CONDENSER: Inlet Temperature", "CONDENSER: Inlet Pressure",
                     "CONDENSER: Inlet Specific Enthalpy", "CONDENSER: Inlet Specific Entropy", "CONDENSER: Inlet Density",
                     "CONDENSER: Outlet Temperature", "CONDENSER: Outlet Pressure",
                     "CONDENSER: Outlet Specific Enthalpy", "CONDENSER: Outlet Specific Entropy", "CONDENSER: Outlet Density",
                     "EXPANSION VALVE: Mass Flow", 
                     "EXPANSION VALVE: Inlet Temperature", "EXPANSION VALVE: Inlet Pressure",
                     "EXPANSION VALVE: Inlet Specific Enthalpy", "EXPANSION VALVE: Inlet Specific Entropy", "EXPANSION VALVE: Inlet Density",
                     "EXPANSION VALVE: Outlet Temperature", "EXPANSION VALVE: Outlet Pressure",
                     "EXPANSION VALVE: Outlet Specific Enthalpy", "EXPANSION VALVE: Outlet Specific Entropy", "EXPANSION VALVE: Outlet Density",
                     "PUMP: Work", "PUMP: Mass Flow", 
                     "PUMP: Inlet Temperature", "PUMP: Inlet Pressure",
                     "PUMP: Inlet Specific Enthalpy", "PUMP: Inlet Specific Entropy", "PUMP: Inlet Density",
                     "PUMP: Outlet Temperature", "PUMP: Outlet Pressure",
                     "PUMP: Outlet Specific Enthalpy", "PUMP: Outlet Specific Entropy", "PUMP: Outlet Density",
                     "SOLVER: Converged (all)", "SOLVER: Max Residual",
                     "SOLVER: Residual P_t", "SOLVER: Residual P_p1", "SOLVER: Residual P_const",
                     "SOLVER: Residual rho_4", "SOLVER: Residual P_5",
                     "SOLVER: Iterations P_t", "SOLVER: Iterations P_p1", "SOLVER: Iterations P_const",
                     "SOLVER: Iterations rho_4", "SOLVER: Iterations P_5",
                     ]
        param_dropdown = ft.Dropdown(
            expand=2,
            menu_height=400,
            hint_text='Dependent Variable',
            hint_style=ft.TextStyle(color="#677E80", italic=True, size = 14, letter_spacing=1.5),
            options=[ft.DropdownOption(key=param, text=param, style=ft.ButtonStyle(color="#B7C3C1")) for param in param_list],
            color="#E3EDED",
            bgcolor="#011314",
            fill_color="#172B2B",
            filled=True,
            text_style=ft.TextStyle(font_family="Open Sans Light", size=14, letter_spacing=1.5),
            dense=True,
            editable=True,
            enable_filter=True,
            enable_search=True,
            border_radius=ft.BorderRadius(top_left=10, top_right=10, bottom_left=10, bottom_right=10),
            trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
            selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
            focused_border_color="#B6FF57",
            focused_border_width=1.5,
            border_color="#1F3835",
            tooltip=tooltip_msg("Select the parameter"),
            on_select= validate,
            on_blur=validate,
            on_text_change= validate
        )

        row = ft.Row(expand=True, controls=[], spacing=2)
        delete_button = ft.IconButton(
            icon=ft.Icons.DELETE_FOREVER_ROUNDED,
            icon_color="#FF5C5C",
            icon_size=30,
            on_click=remove_row,
            data=row,  # a linha vai carregar a referência a si mesma
        )
        
        row = ft.Row(expand=True, data=row_id, controls=[param_dropdown, delete_button], spacing=2)
        return row
    
    def calculate_combinations():
        total = 1

        for r in independent_rows_column.controls:
            equipament_dropdown = r.controls[0]
            param_dropdown = r.controls[1]
            start = r.controls[2]
            end = r.controls[3]
            n_points_field = r.controls[4].content  # TextField dentro do Container

            start_value = start.controls[0].value
            end_value = end.controls[0].value

            if (
                not equipament_dropdown.value
                or not param_dropdown.value
                or not start_value
                or not end_value
                or not n_points_field.value
            ):
                return "-"

            total *= int(n_points_field.value)

        return total
    
    combinations = None
    combinations_text = texto_com_destaque("-", "NUMBER OF SIMULATIONS")
    
    async def validate_n_points(e):

        text = e.control.value

        if not text:
            e.control.error_text = None
            e.control.update()
            return

        # Tenta converter para número
        try:
            value = float(text)
        except ValueError:
            e.control.value = ""
            e.control.update()
            return

        # Verifica se é um número inteiro 
        if not value.is_integer():
            e.control.value = ""
            e.control.update()
            return

        # Verifica se é maior que zero
        if int(value) <= 0:
            e.control.value = ""
            e.control.update()
            return

        # Válido
        e.control.error_text = None
        e.control.value = str(int(value))
        e.control.update()
        nonlocal combinations
        combinations = calculate_combinations()
        try:
            combinations_text.controls[0].content.controls[1].controls[1].value = f"{combinations:,}"
        except:
            combinations_text.controls[0].content.controls[1].controls[1].value = f"{combinations}"
        page.update()

    def add_dependent_row(e):
        new_row = create_dependent_row()
        dependent_rows_column.controls.append(new_row)
        dependent_rows_column.update()

    def create_independent_row():
        row_id = str(uuid.uuid4())
        def remove_row(e):
            for r in independent_rows_column.controls:
                if r.data == row_id:
                    independent_rows_column.controls.remove(r)
                    break
            independent_rows_column.update()

        def update_param_list(e):
            equip = equipament_dropdown.value

            if equip == "Condenser":
                options_list = heater_cooler_params_for_flash(c_flash_type_dropdown.value)
            elif equip == "Evaporator":
                options_list = heater_cooler_params_for_flash(e_flash_type_dropdown.value)
            elif equip == "Generator":
                options_list = heater_cooler_params_for_flash(g_flash_type_dropdown.value)
            else:
                options_list = components_param.get(equip, [])

            param_dropdown.options = [ft.DropdownOption(key=param, text=param, style=ft.ButtonStyle(color="#B7C3C1")) for param in options_list]

            param_dropdown.value = None
            empty_start = ValorUnidade3("", "start value", param_name=None)
            empty_start.expand = 1
            empty_end = ValorUnidade3("", "end value", param_name=None)
            empty_end.expand = 1
            row.controls[2] = empty_start
            row.controls[3] = empty_end
            page.update()
        
        async def check_param(e):
            current_equip = equipament_dropdown.value
            current_param = param_dropdown.value

            if not current_equip or not current_param:
                return False

            for r in independent_rows_column.controls:
                row_equip = r.controls[0].value
                row_param = r.controls[1].value

                if r.data != row_id and row_equip == current_equip and row_param == current_param:
                    await mostrar_erro(f"{current_equip}: {current_param} has already been added.")
                    param_dropdown.value = None
                    page.update()
                    return False

            return True
        
        async def update_param(e):
            is_valid = await check_param(e)
            if not is_valid:
                return
            start = ValorUnidade3("","start value", *unidades[param_dropdown.value],param_name=param_dropdown.value)
            start.expand = 4 if unidades[param_dropdown.value] != [] else 1
            row.controls[2] = start
            end = ValorUnidade3("","end value", *unidades[param_dropdown.value],param_name=param_dropdown.value)
            end.expand = 4 if unidades[param_dropdown.value] != [] else 1
            row.controls[3] = end
            page.update()
        
        def heater_cooler_params_for_flash(flash_type):
            par1, par2 = heater_cooler_flash[flash_type][0], heater_cooler_flash[flash_type][1]
            return [par1, par2, "Efficiency"]
        
        components = ["Condenser","Ejector", "Evaporator", "Generator", "Pump"]
        HeaterCooler_param =["Outlet Temperature", "Outlet Pressure", "Outlet Vapor Fraction", "Efficiency"]
        Ejector_param = ["Primary Nozzle Throat Diameter", "Primary Nozzle Exit Diameter", "Constant Area Section Diameter", "Throat Isentropic Efficiency", "Aerodynamic Throat Isentropic Efficiency", "Diffuser Isentropic Efficiency", "Mixing Loss Factor", "Expansion Coefficient"]
        Pump_param = ["Isentropic Efficiency"]
        unidades = {
            "Outlet Temperature": ["°C", "K", "°F"],
            "Outlet Pressure": ["kPa", "Pa", "MPa", "atm", "bar", "psi"],
            "Outlet Vapor Fraction": [],
            "Efficiency": [],
            "Primary Nozzle Throat Diameter": ["mm", "m", "cm"],
            "Primary Nozzle Exit Diameter" : ["mm", "m", "cm"],
            "Constant Area Section Diameter": ["mm", "m", "cm"],
            "Throat Isentropic Efficiency": [],
            "Aerodynamic Throat Isentropic Efficiency": [],
            "Diffuser Isentropic Efficiency" : [],
            "Mixing Loss Factor": [],
            "Expansion Coefficient": [],
            "Isentropic Efficiency": []
        }
        components_param = {
            "Ejector": Ejector_param,
            "Pump": Pump_param
        }
        equipament_dropdown = ft.Dropdown(
            expand=2,
            hint_text='System Component',
            hint_style=ft.TextStyle(color="#677E80", italic=True, size = 14, letter_spacing=1.5),
            options=[ft.DropdownOption(key=component, text=component, style=ft.ButtonStyle(color="#B7C3C1")) for component in components],
            color="#E3EDED",
            bgcolor="#011314",
            fill_color="#172B2B",
            filled=True,
            text_style=ft.TextStyle(font_family="Open Sans Light", size=14, letter_spacing=1.5),
            dense=True,
            editable=False,
            enable_filter=True,
            enable_search=True,
            border_radius=ft.BorderRadius(top_left=10, top_right=0, bottom_left=10, bottom_right=0),
            trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
            selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
            focused_border_color="#B6FF57",
            focused_border_width=1.5,
            border_color="#1F3835",
            tooltip=tooltip_msg("Select the system component"),
            on_select= update_param_list
        )
        param_dropdown = ft.Dropdown(
            expand=2,
            hint_text='Parameter',
            hint_style=ft.TextStyle(color="#677E80", italic=True, size = 14, letter_spacing=1.5),
            options=[],
            color="#E3EDED",
            bgcolor="#011314",
            fill_color="#172B2B",
            filled=True,
            text_style=ft.TextStyle(font_family="Open Sans Light", size=14, letter_spacing=1.5),
            dense=True,
            editable=False,
            enable_filter=True,
            enable_search=True,
            border_radius=ft.BorderRadius(top_left=0, top_right=0, bottom_left=0, bottom_right=0),
            trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_DOWN, color='#E3EDED'),
            selected_trailing_icon=ft.Icon(ft.Icons.ARROW_DROP_UP, color='#E3EDED'),
            focused_border_color="#B6FF57",
            focused_border_width=1.5,
            border_color="#1F3835",
            tooltip=tooltip_msg("Select the parameter"),
            on_select=update_param
        )
        start = ValorUnidade3("","start value",param_name=param_dropdown.value)
        start.expand = 1
        end = ValorUnidade3("","end value",param_name=param_dropdown.value)
        end.expand = 1
        n_points_field = ft.TextField(text_align=ft.TextAlign.RIGHT,value= "",hint_text="number of points",hint_style=ft.TextStyle(color="#677E80", italic=True, size=14, letter_spacing=1.5), max_lines=1,multiline=False,text_style=ft.TextStyle(color='#B7C3C1', letter_spacing=1.5),  bgcolor="#172B2B", focused_border_color="#B6FF57", border_color="#1F3835", border_radius=ft.BorderRadius.only(top_right=10, bottom_right=10), border_width=1, focused_border_width=1.5, on_change=validate_n_points)
        n_points = ft.Container(expand=2, content=n_points_field,)


        row = ft.Row(expand=True, controls=[], spacing=2)
        delete_button = ft.IconButton(
            icon=ft.Icons.DELETE_FOREVER_ROUNDED,
            icon_color="#FF5C5C",
            icon_size=30,
            on_click=remove_row,
            data=row,  # a linha vai carregar a referência a si mesma
        )
        
        
        row = ft.Row(expand=True, data=row_id, controls=[equipament_dropdown, param_dropdown, start, end, n_points, delete_button], spacing=2)
        return row

    def add_independent_row(e):
        new_row = create_independent_row()
        independent_rows_column.controls.append(new_row)
        independent_rows_column.update()

    # adiciona a primeira linha por padrão
    independent_rows_column.controls.append(create_independent_row())
    dependent_rows_column.controls.append(create_dependent_row())

    batch_progress_bar = ft.ProgressBar(height=10, width=280, color="#B6FF57", border_radius=2, bgcolor="#011314", value=0)
    batch_progress_text = ft.Text("0 / 0 simulations completed", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5))

    async def request_cancel_batch(e):
        nonlocal batch_cancel_requested
        batch_cancel_requested = True

    overlay_batch_progress = ft.Container(
        bgcolor="#50000000", expand=True, visible=False, blur=10, opacity=0,
        alignment=ft.Alignment.CENTER,
        animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT),
        content=ft.Container(
            width=320, padding=ft.Padding(24, 24, 24, 20),  expand=False,
            bgcolor="#0A2022", border_radius=10, border=ft.Border.all(1, "#50B7C3C1"),
            content=ft.Column(
                tight=True,
                alignment=ft.MainAxisAlignment.START,
                horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                spacing=0,
                controls=[
                    ft.Text("Running Batch Simulation...", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5), text_align=ft.TextAlign.CENTER),
                    ft.Container(height=16),
                    batch_progress_bar,
                    ft.Container(height=10),
                    batch_progress_text,
                    ft.Container(height=20),
                    ft.Container(
                        width=280,
                        content=ft.Button(
                            icon=ft.Icons.CANCEL_OUTLINED,
                            icon_color="#FF5C5C",
                            content="Cancel",
                            on_click=request_cancel_batch,
                            bgcolor="#172B2B",
                            color="#FF5C5C",
                            style=ft.ButtonStyle(
                                side={
                                    ft.ControlState.DEFAULT: ft.BorderSide(1.5, "#3C1F1F"),
                                    ft.ControlState.HOVERED: ft.BorderSide(1.5, "#FF5C5C"),
                                },
                                shape=ft.RoundedRectangleBorder(radius=5),
                                text_style=ft.TextStyle(letter_spacing=1.5),
                                padding=ft.Padding.symmetric(vertical=10),
                            ),
                        ),
                    ),
                ]
            )
        )
    )

    async def abre_overlay_batch(e):
        overlay_batch_progress.visible = True
        overlay_batch_progress.opacity = 0
        page.update()
        await asyncio.sleep(0.05)
        overlay_batch_progress.opacity = 1
        page.update()

    async def fecha_overlay_batch(e):
        overlay_batch_progress.opacity = 0
        page.update()
        await asyncio.sleep(0.3)
        overlay_batch_progress.visible = False
        page.update()

    def montar_linhas_tabela():
        n_indep = len(batch_independent_labels)
        colunas = []
        for i, label in enumerate(batch_independent_labels):
            unit = batch_independent_units[i]
            colunas.append(f"{label} ({unit})" if unit else label)


        for dep in batch_dependent_params:
            unit = si_unit_label(dep)
            colunas.append(f"{dep} ({unit})" if unit else dep)

        colunas.append("Status")

        total_linhas = len(batch_results)
        linhas = []
        for row_i in range(total_linhas):
            indep_vals = [batch_independent_si_columns[j][row_i] for j in range(n_indep)]

            resultado = batch_results[row_i]
            if resultado is None or "error" in resultado:
                status = resultado.get("error", "ERROR") if resultado else "ERROR"
                dep_vals = [np.nan] * len(batch_dependent_params)
            else:
                status = "OK"
                dep_vals = [resultado.get(dep, np.nan) for dep in batch_dependent_params]

            linha = indep_vals  + dep_vals + [status]
            linhas.append(linha)
        return colunas, linhas

    def montar_coluna_header(i, label):
        return ft.Text(label, style=ft.TextStyle(color="#B7C3C1", weight=ft.FontWeight.BOLD))

    def atualizar_tabela_pagina(pagina):
        total_linhas = len(batch_table_rows)
        total_paginas = max(1, -(-total_linhas // PAGE_SIZE))
        inicio = pagina * PAGE_SIZE
        fim = min(inicio + PAGE_SIZE, total_linhas)

        datatable.content = fdt.DataTable2(
            expand=True,
            min_width=len(batch_table_columns) * 150,
            column_spacing=20,
            heading_row_color="#0A2022",
            columns=[
                fdt.DataColumn2(label=montar_coluna_header(i, col))
                for i, col in enumerate(batch_table_columns)
            ],
            rows=[
                ft.DataRow(cells=[ft.DataCell(ft.Text(formatar_valor(v), style=ft.TextStyle(color="#B7C3C1"))) for v in linha])
                for linha in batch_table_rows[inicio:fim]
            ],
        )
        page_indicator_text.value = f"Page {pagina + 1} / {total_paginas}  ({inicio + 1}-{fim} of {total_linhas:,})"
        page.update()

    def pagina_anterior(e):
        nonlocal current_page
        if current_page > 0:
            current_page -= 1
            atualizar_tabela_pagina(current_page)

    def proxima_pagina(e):
        nonlocal current_page
        total_paginas = max(1, -(-len(batch_table_rows) // PAGE_SIZE))
        if current_page < total_paginas - 1:
            current_page += 1
            atualizar_tabela_pagina(current_page)

    def formatar_valor(v):
        if isinstance(v, float):
            if np.isnan(v):
                return "-"
            return f"{v:.4f}"
        return str(v)

    back_button = ft.Button(content = "Back", icon= ft.Icons.KEYBOARD_DOUBLE_ARROW_LEFT_ROUNDED, icon_color= "#B6FF57",
                  style = ft.ButtonStyle(color = "#B7C3C1", shape=ft.RoundedRectangleBorder(radius=5), bgcolor=ft.Colors.TRANSPARENT, shadow_color=ft.Colors.TRANSPARENT, text_style=ft.TextStyle(letter_spacing=1.5)), 
                  on_click= back_to_fluxogram)

    batch_sim_content =  ft.Column(spacing = 5, margin=ft.Margin.only(left=5, right=5), controls = [
        back_button,
        ft.Row(expand = 2, spacing = 5,controls=[
            ft.Container(expand = 2,
                bgcolor = "#0A2022",
                border_radius = 10,
                margin = 0,
                padding = 15,
                content = 
                ft.Column(expand=True, 
                    scroll=ft.ScrollMode.AUTO,
                    margin = 0,
                    horizontal_alignment=ft.MainAxisAlignment.CENTER,
                    controls=[
                    ft.Row(expand = False, controls=[
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2),
                        ft.Text("INDEPENDENT VARIABLES", style=ft.TextStyle(font_family="Open Sans Light", size=20, letter_spacing=2, color="#B7C3C1")),
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2),
                    ]),
                    independent_rows_column,
                    ft.IconButton(icon = ft.Icons.ADD_CIRCLE_ROUNDED, icon_color= "#B6FF57", icon_size = 32, on_click = add_independent_row, key = "independent var", tooltip= tooltip_msg("add a new independent variable", 500)),
                ])
            ),
            ft.Container(expand = 1, 
                bgcolor = "#0A2022",
                border_radius = 10,
                margin = 0,
                padding = 15,
                content = 
                ft.Column(expand=True,
                    scroll=ft.ScrollMode.AUTO,
                    margin = 0,
                    horizontal_alignment=ft.MainAxisAlignment.CENTER,
                    controls=[
                    ft.Row(expand = False, controls=[
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2),
                        ft.Text("DEPENDENT VARIABLES", style=ft.TextStyle(font_family="Open Sans Light", size=20, letter_spacing=2, color="#B7C3C1")),
                        ft.Container(expand= True, height=5, bgcolor="#B6FF57", border_radius=2),
                    ]),
                    dependent_rows_column,
                    ft.IconButton(icon = ft.Icons.ADD_CIRCLE_ROUNDED, icon_color= "#B6FF57", icon_size = 32, on_click = add_dependent_row, key = "dependent var", tooltip= tooltip_msg("add a new dependent variable", 500)),
                ])
            ),
        ]),
    ])

    datatable = ft.Container(expand=True)
    page_indicator_text = ft.Text("", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.2))

    back_button_results = ft.Button(content="Back", icon=ft.Icons.KEYBOARD_DOUBLE_ARROW_LEFT_ROUNDED, icon_color="#B6FF57",
    style=ft.ButtonStyle(color="#B7C3C1", shape=ft.RoundedRectangleBorder(radius=5), bgcolor=ft.Colors.TRANSPARENT, shadow_color=ft.Colors.TRANSPARENT, text_style=ft.TextStyle(letter_spacing=1.5)),
    on_click=back_to_fluxogram)

    overlay_exporting = ft.Container(bgcolor="#50000000", expand=True, visible=False, blur=10, opacity=0, alignment=ft.Alignment.CENTER, animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT),
        content=ft.Container(alignment=ft.Alignment.CENTER, height=150, width=300, bgcolor="#0A2022", border_radius=10, border=ft.Border.all(1, "#50B7C3C1"),
            content=ft.Column(alignment=ft.MainAxisAlignment.CENTER, horizontal_alignment=ft.CrossAxisAlignment.CENTER,
                controls=[
                    ft.Text("Exporting table...", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5)),
                    ft.ProgressBar(height=10, width=280, color="#B6FF57", border_radius=2, bgcolor="#011314"),
                ])))

    async def abre_overlay_exporting(e):
        overlay_exporting.visible = True
        overlay_exporting.opacity = 0
        page.update()
        await asyncio.sleep(0.05)
        overlay_exporting.opacity = 1
        page.update()

    async def fecha_overlay_exporting(e):
        overlay_exporting.opacity = 0
        page.update()
        await asyncio.sleep(0.3)
        overlay_exporting.visible = False
        page.update()

    async def export_to_csv(e):
        if not batch_table_columns or not batch_table_rows:
            await mostrar_erro("There are no Batch Sim results to export.")
            return

        file_path = await ft.FilePicker().save_file(
            dialog_title="Export results to CSV",
            file_name="erc_batch_results.csv",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["csv"],
        )
        if not file_path:
            return
        if not file_path.lower().endswith(".csv"):
            file_path += ".csv"

        await abre_overlay_exporting(e)
        try:
            await asyncio.get_event_loop().run_in_executor(
                None, _write_csv_file, file_path, batch_table_columns, batch_table_rows
            )
        except Exception as error:
            await fecha_overlay_exporting(e)
            await mostrar_erro(f"An error occurred: {error}")
            return
        await fecha_overlay_exporting(e)

        await mostrar_erro(
            "Table exported successfully!",
            icone=ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size=50),
        )

    async def export_to_xlsx(e):
        if not batch_table_columns or not batch_table_rows:
            await mostrar_erro("There are no Batch Sim results to export.")
            return

        file_path = await ft.FilePicker().save_file(
            dialog_title="Export results to Excel",
            file_name="erc_batch_results.xlsx",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["xlsx"],
        )
        if not file_path:
            return
        if not file_path.lower().endswith(".xlsx"):
            file_path += ".xlsx"

        await abre_overlay_exporting(e)
        try:
            await asyncio.get_event_loop().run_in_executor(
                None, _write_xlsx_file, file_path, batch_table_columns, batch_table_rows
            )
        except Exception as error:
            await fecha_overlay_exporting(e)
            await mostrar_erro(f"An error occurred: {error}")
            return
        await fecha_overlay_exporting(e)

        await mostrar_erro(
            "Table exported successfully!",
            icone=ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size=50),
        )

    async def export_to_parquet(e):
        if not batch_table_columns or not batch_table_rows:
            await mostrar_erro("There are no Batch Sim results to export.")
            return

        file_path = await ft.FilePicker().save_file(
            dialog_title="Export results to Parquet",
            file_name="erc_batch_results.parquet",
            file_type=ft.FilePickerFileType.CUSTOM,
            allowed_extensions=["parquet"],
        )
        if not file_path:
            return
        if not file_path.lower().endswith(".parquet"):
            file_path += ".parquet"

        await abre_overlay_exporting(e)
        try:
            await asyncio.get_event_loop().run_in_executor(
                None, _write_parquet_file, file_path, batch_table_columns, batch_table_rows
            )
        except Exception as error:
            await fecha_overlay_exporting(e)
            await mostrar_erro(f"An error occurred: {error}")
            return
        await fecha_overlay_exporting(e)

        await mostrar_erro(
            "Table exported successfully!",
            icone=ft.Icon(icon=ft.Icons.CHECK_ROUNDED, color="#B6FF57", size=50),
        )

    export_results_button = ft.MenuBar(
        style=ft.MenuStyle(shadow_color=ft.Colors.TRANSPARENT,
        bgcolor="#B6FF57",shape=ft.RoundedRectangleBorder(radius=10)                                             ),
            controls=[
            ft.SubmenuButton(menu_style=ft.MenuStyle(bgcolor="#011314"),content=ft.Text("Export Table", style=ft.TextStyle(font_family="Open Sans Regular", letter_spacing=1.5, color="#011314", weight=ft.FontWeight.BOLD)), controls=[
                ft.MenuItemButton(content=ft.Text("to .csv",style=ft.TextStyle(font_family="Open Sans Light", letter_spacing=1.5, color="#B7C3C1")), style=ft.ButtonStyle(bgcolor="#011314"), on_click=export_to_csv),
                ft.MenuItemButton(content=ft.Text("to .xlsx",style=ft.TextStyle(font_family="Open Sans Light", letter_spacing=1.5, color="#B7C3C1")), style=ft.ButtonStyle(bgcolor="#011314"), on_click=export_to_xlsx),
                ft.MenuItemButton(content=ft.Text("to .parquet",style=ft.TextStyle(font_family="Open Sans Light", letter_spacing=1.5, color="#B7C3C1")), style=ft.ButtonStyle(bgcolor="#011314"), on_click=export_to_parquet)])
    ])

    results_content = ft.Column(spacing=5, margin=ft.Margin.only(left=5, right=5), controls=[
        ft.Row(alignment=ft.MainAxisAlignment.SPACE_BETWEEN, expand = False, controls=[
            ft.Container(content=back_button_results),
            ft.Container(content=export_results_button)
        ]),
        ft.Row(expand=2, spacing=5, controls=[
            ft.Container(expand=2,
                bgcolor="#0A2022", border_radius=10, margin=0, padding=15,
                content=ft.Column(expand=True,
                    scroll=ft.ScrollMode.AUTO, margin=0,
                    horizontal_alignment=ft.MainAxisAlignment.CENTER,
                    controls=[
                        ft.Row(expand=False, controls=[
                            ft.Container(expand=True, height=5, bgcolor="#B6FF57", border_radius=2),
                            ft.Text("BATCH SIM RESULTS (SI UNITS)", style=ft.TextStyle(font_family="Open Sans Light", size=20, letter_spacing=2, color="#B7C3C1")),
                            ft.Container(expand=True, height=5, bgcolor="#B6FF57", border_radius=2),
                        ]),
                        datatable,
                        ft.Row(alignment=ft.MainAxisAlignment.CENTER, spacing=10, controls=[
                            ft.IconButton(icon=ft.Icons.CHEVRON_LEFT_ROUNDED, icon_color="#B6FF57", on_click=pagina_anterior),
                            page_indicator_text,
                            ft.IconButton(icon=ft.Icons.CHEVRON_RIGHT_ROUNDED, icon_color="#B6FF57", on_click=proxima_pagina),
                        ]),
                    ])
            )
        ]),
    ])

    #=================================================
    animador_de_telas = ft.AnimatedSwitcher(
        content=ft.Container(key=""), 
        transition=ft.AnimatedSwitcherTransition.FADE, # Efeito de sumir e aparecer
        duration=300, # Duração da transição em milissegundos
        switch_in_curve=ft.AnimationCurve.EASE_IN,
        switch_out_curve=ft.AnimationCurve.LINEAR,
        reverse_duration=100
    )

    # 1. Painel Esquerdo isolado com tamanho fixo e rígido
    conteudo_esquerdo = ft.Container(
        width=LARGURA_PAINEL,
        bgcolor="#0A2022",
        #bgcolor=ft.Colors.TRANSPARENT,
        border_radius=10,
        padding=10,
        content=animador_de_telas,
    )

    # Coloca o painel dentro de um container com animação de posição ativa
    painel_esquerdo = ft.Container(
        content=conteudo_esquerdo,
        left=-(LARGURA_PAINEL+20),  # Posição X inicial no Stack
        top=0,   # Força a esticar para o topo
        bottom=0, # Força a esticar para a base
        animate_position=ft.Animation(400, ft.AnimationCurve.EASE_OUT_CUBIC)
    )

    def gradiente_padrao():
        return ft.LinearGradient(colors=['#00B3F364', '#00FFFFFF'], begin=ft.Alignment.TOP_LEFT, end=ft.Alignment.BOTTOM_RIGHT, stops=[0,0.8])

    def borda_padrao():
        return ft.Border.all(0.5, color="#00B6FF57")

    # Dimensões do Stack no momento em que as coordenadas foram calibradas
    IMG_W = 751.245930599369
    IMG_H = 700.8

    # Frações proporcionais de cada overlay do stack erc
    OVERLAY_DEFS = {
        'pump':            (115/IMG_H,     20/IMG_W,  140/IMG_H,  120/IMG_W),
        'generator':       (  5/IMG_H,    284/IMG_W,  160/IMG_H,  180/IMG_W),
        'condenser':       (265/IMG_H,    162/IMG_W,  105/IMG_H,  270/IMG_W),
        'ejector':         (265/IMG_H,    480/IMG_W,  105/IMG_H,  195/IMG_W),
        'expansion_valve': (375/IMG_H,      5/IMG_W,  190/IMG_H,  165/IMG_W),
        'evaporator':      (None,        302.5/IMG_W,  165/IMG_H,  240/IMG_W),
    }
    OVERLAY_BORDER_RADIUS = {
        'pump':            10,
        'generator':       15,
        'condenser':       10,
        'ejector':         10,
        'expansion_valve': 10,
        'evaporator':      15,
    }

    def recalcular_overlays(stack_w: float, stack_h: float):
        scale = (stack_w / IMG_W + stack_h / IMG_H) / 2  # fator médio de escala
        overlay_keys = ['pump', 'generator', 'condenser', 'ejector', 'expansion_valve', 'evaporator']
        for i, key in enumerate(overlay_keys, start=1):
            top_pct, left_pct, h_pct, w_pct = OVERLAY_DEFS[key]
            c = erc.controls[i]
            c.left          = left_pct * stack_w
            c.width         = w_pct   * stack_w
            c.height        = h_pct   * stack_h
            c.border_radius = OVERLAY_BORDER_RADIUS[key] * scale  
            if key == 'evaporator':
                c.bottom = 0
                c.top    = None
            else:
                c.top    = top_pct * stack_h
                c.bottom = None
        erc.update()

    def on_erc_resize(e):
        recalcular_overlays(e.width, e.height) #será que dava pra ter usado lambda?
    
    erc = ft.Stack(
        controls=[
            ft.Image("/images/ERC.png", fit=ft.BoxFit.COVER, margin=10),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='pump',            border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Pump")),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='generator',       border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Generator")),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='condenser',       border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Condenser")),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='ejector',         border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Ejector")),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='expansion_valve', border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Expansion Valve")),
            ft.Container(gradient=gradiente_padrao(), on_click=gerenciar_clique, data='evaporator',      border=borda_padrao(), animate=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT), tooltip=tooltip_msg("Evaporator")),
        ],
        margin=10,
        on_size_change=on_erc_resize,
    )

    coluna = ft.Column(
        spacing = 5,
        controls=[
            # Barra superior (fluido, bake, run, etc)
            ft.Container(content=
                    ft.Row(
                        expand = True,
                        controls=[dropdown, bake_button, backend, batch_sim_button, results_button, run_simulation_button, settings_button],
                        alignment=ft.MainAxisAlignment.START,
                        spacing=2,
                        vertical_alignment=ft.CrossAxisAlignment.CENTER,
                        
                    ),bgcolor="#0A2022", padding=0,
            ),
            
            # ERC, tabelas etc
            ft.AnimatedSwitcher(
                transition=ft.AnimatedSwitcherTransition.FADE,
                duration=300,
                switch_in_curve=ft.AnimationCurve.EASE_IN,
                switch_out_curve=ft.AnimationCurve.LINEAR,
                reverse_duration=100,
                expand=True,
                content = ft.Container(
                    key ="erc",
                    content=erc,
                    alignment=ft.Alignment.CENTER,
                    expand=True,
                    data='out_click',
                    on_click=gerenciar_clique,
                    image= ft.DecorationImage(src="/images/pattern3.webp",repeat=ft.ImageRepeat.REPEAT,fit=ft.BoxFit.NONE, scale=20, opacity=0.5)
                )
            )
        ]   
    )
    # 2. Painel Direito usando MARGIN para abrir/fechar o espaço dinamicamente
    
    painel_direito = ft.Container(
        content=coluna,
        bgcolor="#041A1C",
        #border_radius=10,
        border_radius = ft.BorderRadius(top_left=10, top_right=10, bottom_left=10, bottom_right=10),
        alignment=ft.Alignment.CENTER,
        # Começa com margem para não ficar embaixo do painel esquerdo
        margin=ft.Margin.only(left=0),#left=LARGURA_PAINEL + SPACING
        animate=ft.Animation(400, ft.AnimationCurve.EASE_OUT_CUBIC),
        expand=True
    )


    # 3. Montagem estrutural no Stack
    stack_paineis = ft.Stack(
            expand=True,
            controls=[
                painel_direito,            # Camada de baixo (ajusta via margin)
                painel_esquerdo,   # Camada de cima (desliza via left de cima a baixo)
            ]
    )

    overlay_bake_fluid = ft.Container(bgcolor="#50000000", expand = True, visible=False, blur= 10, opacity=0, alignment= ft.Alignment.CENTER, animate_opacity=ft.Animation(300, ft.AnimationCurve.EASE_IN_OUT),
                           content=ft.Container(alignment= ft.Alignment.CENTER,height=150, width=300, bgcolor="#0A2022", border_radius=10,border =ft.Border.all(1,"#50B7C3C1"),
                                                content=ft.Column(alignment=ft.MainAxisAlignment.CENTER, horizontal_alignment=ft.CrossAxisAlignment.CENTER, 
                                                                  controls=[ft.Text("baking fluid...", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5)), ft.ProgressBar(height=10, width= 280, color="#B6FF57", border_radius=2, bgcolor="#011314"),ft.Text("(this takes around 1 min)", style=ft.TextStyle(color="#B7C3C1", letter_spacing=1.5))])))


    page.add(
        stack_paineis,
        ft.Text(f"SIMERC v.{VERSAO} - a software by GPCTM", color="#677E80", size=10, style=ft.TextStyle(font_family='Open Sans Regular'))
    )
    recalcular_overlays(IMG_W,IMG_H)
    page.overlay.append(overlay_bake_fluid)
    page.overlay.append(overlay_erro)
    page.overlay.append(overlay_batch_progress)
    page.overlay.append(overlay_settings)
    page.overlay.append(overlay_exporting)
    page.update()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    ft.run(main, assets_dir = obter_caminho_recurso("assets"))