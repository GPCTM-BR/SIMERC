import bpy
import numpy as np
import pyarrow.parquet as pq

#Rode esse código DENTRO do Blender, e certifique-se de mudar as colunas para as variaveis que você escolheu

arquivo = r"CAMINHO PARA O ARQUIVO .parquet"


#aqui estão as variáveis independentes (3 primeiras) e dependentes (resto) que eu escolhi observar. Troque pelas suas e depois reveja o resto do código (ou só usa as que eu usei)
colunas = [
    "Condenser: Outlet Temperature (K)",
    "Generator: Outlet Temperature (K)",
    "Evaporator: Outlet Temperature (K)",
    "COP",
    "EJECTOR: Entrainment Ratio",
    "EJECTOR: Pressure Lift Ratio",
    "EJECTOR: Critical Backpressure (Pa)",
    "EJECTOR: Primary Inlet Mass Flow (kg/s)",
    "EJECTOR: Secondary Inlet Mass Flow (kg/s)",
    "EJECTOR: Outlet Mass Flow (kg/s)",
    "EJECTOR: Primary Inlet Pressure (Pa)",
    "EJECTOR: Primary Inlet Specific Enthalpy (J/kg)",
    "EJECTOR: Secondary Inlet Pressure (Pa)",
    "EJECTOR: Secondary Inlet Specific Enthalpy (J/kg)",
    "EJECTOR: Outlet Pressure (Pa)",
    "EJECTOR: Outlet Specific Enthalpy (J/kg)",
    "GENERATOR: Heaty Duty (W)",
    "GENERATOR: Mass Flow (kg/s)",
    "GENERATOR: Inlet Pressure (Pa)",
    "GENERATOR: Inlet Specific Enthalpy (J/kg)",
    "GENERATOR: Outlet Pressure (Pa)",
    "GENERATOR: Outlet Specific Enthalpy (J/kg)",
    "CONDENSER: Heaty Duty (W)",
    "CONDENSER: Mass Flow (kg/s)",
    "CONDENSER: Inlet Pressure (Pa)",
    "CONDENSER: Inlet Specific Enthalpy (J/kg)",
    "CONDENSER: Outlet Pressure (Pa)",
    "CONDENSER: Outlet Specific Enthalpy (J/kg)",
    "EVAPORATOR: Heaty Duty (W)",
    "EVAPORATOR: Mass Flow (kg/s)",
    "EVAPORATOR: Inlet Pressure (Pa)",
    "EVAPORATOR: Inlet Specific Enthalpy (J/kg)",
    "EVAPORATOR: Outlet Pressure (Pa)",
    "EVAPORATOR: Outlet Specific Enthalpy (J/kg)",
    "PUMP: Work (W)",
]

chunk_size = 200000
lista_coords = []
lista_erros = []
dicionario_atributos = {col: [] for col in colunas}
erros_mapeamento = {}
codigo_erro_atual = 1

print("Iniciando leitura do arquivo parquet em blocos...")


colunas_ler = colunas + ["Status"]
parquet_file = pq.ParquetFile(arquivo)

for i, batch in enumerate(parquet_file.iter_batches(batch_size=chunk_size, columns=colunas_ler)):
    print(f"Processando bloco {i+1}...")
    chunk = batch.to_pandas()
    status_col = chunk["Status"].to_numpy()
    is_erro = status_col != "OK"
    erro_codigo_bloco = np.zeros(len(chunk), dtype=np.float32)

    if np.any(is_erro):
        for idx in np.where(is_erro)[0]:
            texto = str(status_col[idx]).strip()
            if "Current value" in texto:
                texto = texto.split("Current value")[0].strip()

            if texto not in erros_mapeamento:
                erros_mapeamento[texto] = codigo_erro_atual
                codigo_erro_atual += 1

            erro_codigo_bloco[idx] = erros_mapeamento[texto]

    lista_erros.append(erro_codigo_bloco)
    chunk_num = chunk[colunas].astype(np.float32).fillna(-1.0)

    # Extracao das coordenadas (X, Y, Z) deste bloco
    x = chunk_num.iloc[:, 0].to_numpy()
    y = chunk_num.iloc[:, 1].to_numpy()
    z = chunk_num.iloc[:, 2].to_numpy()
    lista_coords.append(np.column_stack((x, y, z)))

    # Armazena os dados numericos para os atributos do Blender
    for col in colunas:
        dicionario_atributos[col].append(chunk_num[col].to_numpy())

print("\n=== CÓDIGOS DOS ERROS ENCONTRADOS ===")
for nome, cod in erros_mapeamento.items():
    print(f"{cod}: {nome}")

print("\nConsolidando matrizes de dados...")
coords = np.vstack(lista_coords)
erros_finais = np.concatenate(lista_erros)

# Libera memória das listas temporárias
del lista_coords
del lista_erros

# Normalização das coordenadas das linhas válidas
mins = coords.min(axis=0)
maxs = coords.max(axis=0)
# Evita divisão por zero se o max for igual ao min
maxs[maxs == mins] += 1e-6
coords = (coords - mins) / (maxs - mins)

total_pontos = len(coords)
print(f"Total de vértices a criar no Blender: {total_pontos}")


# Criação da Malha (Mesh) no Blender
mesh = bpy.data.meshes.new("Scatter")
obj = bpy.data.objects.new("Scatter", mesh)
bpy.context.collection.objects.link(obj)
mesh.vertices.add(total_pontos)
mesh.vertices.foreach_set("co", coords.ravel())
mesh.update()


print("Injetando atributos de erro e COP...")
atrib_erro = mesh.attributes.new(name="Erro", type='FLOAT', domain='POINT')
atrib_erro.data.foreach_set("value", erros_finais.astype(np.float32))
del erros_finais

print("Injetando demais atributos numéricos...")
for col_nome in colunas:
    valores_coluna = np.concatenate(dicionario_atributos[col_nome]).astype(np.float32)

    atrib = mesh.attributes.new(name=col_nome, type='FLOAT', domain='POINT')
    atrib.data.foreach_set("value", valores_coluna)

    dicionario_atributos[col_nome] = None

mesh.update()
print("\nMesh criada com sucesso a partir do arquivo parquet.")
