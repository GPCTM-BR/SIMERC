# ![SIMERC - EJECTOR REFRIGERATION CYCLE SIMULATOR](docstrings/images/LOGO.png)

## SIMERC - EJECTOR REFRIGERATION CYCLE SIMULATOR
  
O SIMERC é um programa desenvolvido em Python como parte de um Trabalho de Conclusão de Curso de Engenharia Química. Ele permite simular ciclos de refrigeração com ejetores, bem como realizar simulações em série para teste de vários casos e exportar os resultados.  
O [CoolProp](https://coolprop.org/) é utilizado para os cálculos das propriedades termodinâmicas e o ejetor segue o modelo desenvolvido por [Cardemil e Colle (2012)](https://doi.org/10.1016/j.enconman.2012.05.009).
  

Visual do Programa:
![Print 01 - GUI do programa](docstrings/images/print1.png)
<div align=center>
    Aba principal
</div>
  
Entre as características desse programa destacam-se:
- Utilização do [CoolProp](https://coolprop.org/), que fornece resultados precisos para a maioria dos fluidos refrigerantes
- Dois backends do [CoolProp](https://coolprop.org/) foram implementados: O HEOS (padrão) e  [BICUBIC&HEOS](https://coolprop.org/coolprop/Tabular.html) (para cálculos mais rápidos)
- O código responsável pela parte matemática da simulação foi escrito com [Cython](https://cython.org/), aumentando muito a performance do programa
- O usuário pode realizar a simulação em batch utilizando-se de processamento paralelo, útil para casos com muitas simulações
- Os resultados do Batch Sim podem ser visualizados no próprio programa, e exportados como csv, xlsx ou parquet
- Os parâmetros $\phi_m$ e $\psi$ podem ser expressões em função de _Ar_ e _Pr_, como no trabalho de [Cardemil e Colle (2012)](https://doi.org/10.1016/j.enconman.2012.05.009).
  
Outros prits do programa:
![Print 02 - GUI do programa](docstrings/images/print2.png)
<div align=center>
    Aba Batch Sim
</div>
  
![Print 03 - GUI do programa](docstrings/images/print3.png)
<div align=center>
    Aba Batch Sim Table
</div>
  
![Print 04 - GUI do programa](docstrings/images/print4.png)
<div align=center>
    Aba de Configurações
</div>
  
  
# DOWNLOAD
