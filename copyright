import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from collections import Counter
import streamlit as st
from PIL import Image
import urllib.request
import io

# Padroniza a sequência e verifica se é DNA:


def formata_sequencia(sequencia): 
    sequencia = sequencia.upper()
        
    if sequencia[0] == ">":
        sequencia = sequencia.splitlines()
        sequencia = sequencia[1:]
        sequencia = "".join(sequencia).strip()
        sequencia = "*" + sequencia
        
    else:
        sequencia = sequencia.splitlines()
        sequencia = "".join(sequencia).strip()
        sequencia = "*" + sequencia 
        
    return sequencia


def eh_dna(seq):
    if set(seq).issubset({"A", "C", "G", "T", "*"}):
        return True
    else:
        return False


# Análise da composição de bases da sequência:

def conta_nucleotideo(seq):
    
    mapa = Counter(seq)

    adenina = mapa["A"]
    citosina = mapa["C"]
    guanina = mapa["G"]
    timina = mapa["T"]
    total = len(seq)-1
    quantidade = [adenina, citosina, guanina, timina, total]

    porc_adenina = (adenina / total) * 100
    porc_citosina = (citosina / total) * 100
    porc_guanina = (guanina / total) * 100
    porc_timina = (timina / total) * 100
    porc_total = porc_adenina + porc_citosina + porc_guanina + porc_timina
    porcentagem = [porc_adenina, porc_citosina, porc_guanina, porc_timina, porc_total]

    return quantidade, porcentagem


def df_quantidade_nucleotideos(seq1, seq2):   

    df = pd.DataFrame(columns=["Feature", "String 1", "String 2"])
  
    quantidade1, _ = conta_nucleotideo(seq1)
    quantidade2, _ = conta_nucleotideo(seq2)

    df["Feature"] = ["Adenin", "Cytosin", "Guanin", "Thymin", "Total"]
    df["String 1"] = quantidade1
    df["String 2"] = quantidade2

    return df


def df_porcentagem_nucleotideos(seq1, seq2):   

    df = pd.DataFrame(columns=["Feature", "String 1", "String 2"])
  
    _, porcentagem1 = conta_nucleotideo(seq1)
    _, porcentagem2 = conta_nucleotideo(seq2)

    df["Feature"] = ["Adenin", "Cytosin", "Guanin", "Thymin", "Total"]
    df["String 1"] = porcentagem1
    df["String 2"] = porcentagem2

    return df


def gera_grafico(df, selecao="Quantity"):

    df_melt = pd.melt(df[:4], id_vars="Feature", var_name="String", value_name=selecao)

    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x='Feature', y=selecao, hue='String', data=df_melt)
    ax.set_title(f"Nucleotides {selecao}", fontweight="bold", fontsize=14)
    ax.set_xlabel("Nucleotide", fontsize=12)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend(bbox_to_anchor=(1.35, 1))
  
    if selecao == "Percentage":
        ax.set_ylim(0, 100)
        ax.set_ylabel("Percentage (%)", fontsize=12)    
    else:
        ax.set_ylabel("Quantity", fontsize=12) 

    return fig


def conteudo_gc(seq1, seq2):
       
    df = pd.DataFrame(columns=["Feature", "String", "Percentage"])

    _, porcentagem1 = conta_nucleotideo(seq1)
    _, porcentagem2 = conta_nucleotideo(seq2)

    gc1 = porcentagem1[2] + porcentagem1[1]
    gc2 = porcentagem2[2] + porcentagem2[1]

    df["Feature"] = ["GC", "GC"]
    df["String"] = ["String 1", "String 2"]
    df["Percentage"] = [gc1, gc2]
  
    return df


def gera_grafico_gc(df):
 
    fig, ax = plt.subplots(figsize=(6,3))
    ax = sns.barplot(x="Feature", y="Percentage", hue="String", data=df)
    ax.set_title(f"GC Content", fontweight="bold", fontsize=14)
    ax.set_ylabel("Percentage (%)", fontsize=12)
    ax.set_ylim(0, 100)
    ax.set_xlabel(None)
    plt.legend(bbox_to_anchor=(1.35, 1))

    return fig


# Alinhamento com Needleman-Wunsch:

def cria_matriz_subs():

     matriz = {"col": ["A", "C", "G", "T"],
                 "A": [4, -2, -1, -2],
                 "C": [-2, 4, -2, -1],
                 "G": [-1, -2, 4, -2],
                 "T": [-2, -1, -2, 4]}

     return matriz


def calcula_score(base1, base2, matriz_subs):

    j = matriz_subs["col"].index(base1)

    for base in "ACGT":
        if base2 == base:
            score = matriz_subs[base2][j]
            return score


def valor_maximo(base1, base2, lado, cima, diagonal):

    if (base1 == base2) and (diagonal > lado) and (diagonal > cima):
        return diagonal
    elif (base1 != base2) and (diagonal > lado) and (diagonal > cima):
        return diagonal
    elif (lado > cima) and (lado > diagonal):
        return lado
    else:
        return cima


def cria_caminho(base1, base2, lado, cima, diagonal):

    if (base1 == base2) and (diagonal > lado) and (diagonal > cima):
        return "\\"
    elif (base1 != base2) and (diagonal > lado) and (diagonal > cima):
        return "\\"
    elif (lado > cima) and (lado > diagonal):
        return "-"
    else:
        return "|"


def lcs_global(seq1, seq2, matriz_subs):

    pontuacao = []
    caminho = []
    g = -3

    for i in range(0, len(seq1)):
        pontuacao.append([0] * len(seq2))
        caminho.append([""] * len(seq2))

    for i in range(0, len(seq1)):
        pontuacao[i][0] = g * i
        caminho[i][0] = "|"
    for j in range(0, len(seq2)):
        pontuacao[0][j] = g * j
        caminho[0][j] = "-"
 
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):  

            base1 = seq1[i]
            base2 = seq2[j] 
            s = calcula_score(base1, base2, matriz_subs)
        
            lado = pontuacao[i][j-1]
            cima = pontuacao[i-1][j]
            diagonal = pontuacao[i-1][j-1]
    
            pontuacao[i][j] = valor_maximo(base1, base2, lado + g, cima + g, diagonal + s)
            caminho[i][j] = cria_caminho(base1, base2, lado + g, cima + g, diagonal + s)

    return caminho


def lcs_local(seq1, seq2, matriz_subs):

    pontuacao = []
    caminho = []
    g = -3

    for i in range(0, len(seq1)):
        pontuacao.append([0] * len(seq2))
        caminho.append([""] * len(seq2))

    for i in range(0, len(seq1)):
        caminho[i][0] = "|"
    for j in range(0, len(seq2)):
        caminho[0][j] = "-"
 
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):  

            base1 = seq1[i]
            base2 = seq2[j] 
            s = calcula_score(base1, base2, matriz_subs)
        
            lado = pontuacao[i][j-1]
            cima = pontuacao[i-1][j]
            diagonal = pontuacao[i-1][j-1]
    
            pontuacao[i][j] = max(0, valor_maximo(base1, base2, lado + g, cima + g, diagonal + s))
            caminho[i][j] = cria_caminho(base1, base2, lado + g, cima + g, diagonal + s)

    return pontuacao, caminho


def alinhamento_global(seq1, seq2, matriz_caminho, matriz_subs):
    ali_seq1 = ""
    ali_seq2 = ""
    g = -3
    match = 0
    mismatch = 0
    gap = 0
    score_final = 0

    i = len(seq1)-1
    j = len(seq2)-1

    while (i != 0) or (j != 0):
        s = calcula_score(seq1[i], seq2[j], matriz_subs)

        if matriz_caminho[i][j] == "\\" and seq1[i] == seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            match += 1
            score_final += s
            i -= 1
            j -= 1

        elif matriz_caminho[i][j] == "\\" and seq1[i] != seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            mismatch += 1
            score_final += s
            i -= 1
            j -= 1    
    
        elif matriz_caminho[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            gap += 1
            score_final += g
            j -= 1
    
        elif matriz_caminho[i][j] == "|":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            gap += 1
            score_final += g
            i -= 1

    return match, mismatch, gap, score_final, ali_seq1, ali_seq2


def procura_maior_valor(seq1, matriz_pontuacao):
  
    maior_valor_matriz = 0

    for i in range(len(seq1)):
        maior_valor_da_linha = max(matriz_pontuacao[i])
        if maior_valor_da_linha > maior_valor_matriz:
            maior_valor_matriz = maior_valor_da_linha
            linha = i
  
    coluna = matriz_pontuacao[linha].index(maior_valor_matriz)
    return linha, coluna


def alinhamento_local(seq1, seq2, matriz_pontuacao, matriz_caminho, matriz_subs):
    ali_seq1 = ""
    ali_seq2 = ""
    g = -3
    match = 0
    mismatch = 0
    gap = 0
    score_final = 0

    i, j = procura_maior_valor(seq1, matriz_pontuacao)
    valor = matriz_pontuacao[i][j]

    while valor > 0:

        s = calcula_score(seq1[i], seq2[j], matriz_subs)    

        if matriz_caminho[i][j] == "\\" and seq1[i] == seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            match += 1
            score_final += s
            i -= 1
            j -= 1
            valor = matriz_pontuacao[i][j]

        elif matriz_caminho[i][j] == "\\" and seq1[i] != seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            mismatch += 1
            score_final += s
            i -= 1
            j -= 1    
            valor = matriz_pontuacao[i][j]
    
        elif matriz_caminho[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            gap += 1
            score_final += g
            j -= 1
            valor = matriz_pontuacao[i][j]
    
        elif matriz_caminho[i][j] == "|":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            gap += 1
            score_final += g
            i -= 1
            valor = matriz_pontuacao[i][j]

    return match, mismatch, gap, score_final, ali_seq1, ali_seq2


# função principal:

def main():

    url = "https://raw.githubusercontent.com/vanleiko/dna-streamlit/main/src/dna-logo2.jpg"
    with urllib.request.urlopen(url) as i:
	    byteImg = io.BytesIO(i.read())
	    imagem = Image.open(byteImg)

    default_input1 = "AAATTTACGTGGTT"
    default_input2 = "ACGTGT"
    help_text = "Free nucleotides sequence or FASTA format"

    st.image(imagem, use_column_width=True)
    st.title("Pairwise DNA sequence alignment")

    st.sidebar.markdown("""<p style='text-align: justify'>
    <b>🧬 About this web app</b><br>
    Global and local alignment between two DNA strings (pairwise alignment).<br><br> 
    <b>⚙️ How it works</b><br>
    <b>Global alignment</b> is an end-to-end alignment based on Needleman-Wunsch algorithm, 
    which searches for the optimal global alignment (highest score) between two sequences.<br> 
    <b>Local alignment</b> identify regions of similarity between two sequences based on Smith-Waterman algorithm, 
    which searches for the optimal alignment (highest score) between two subsequences.<br><br> 
    <b>🔢 Scoring</b><br> 
    The scores for match and mismatch were based on nucleotide substitution model K2P (Kimura 2-parameters)
    which allows for different rates of transition and transversion.<br>
    <i>> Match = +4</i>, for similarity<br>
    <i>> Mismatch = -1</i>, for transition<br> 
    <i>> Mismatch = -2</i>, for transversion<br> 
    <i>> Gap = -3</i>, for insertion or deletion<br><br>
    <b>📊 Nucleotides content</b><br>
    This web app also counts the amount and percentage of Adenin, Cytosin, Guanin and  
    Thymin and GC-Content from each sequence.</b><br><br> 
    <b>📚 References</b><br>
    Needleman, S.B. and Wunsch, C.D. 1970. A general method applicable to the search for similarities in 
    the amino acid sequences of two proteins. J. Mol. Bio., 48:443-453.<br>
    Kimura, M. 1980. A simple method for estimating evolutionary rates of base substitutions through comparative
    studies of nucleotide sequences. J. Mol. Evol., 16:111-120.<br>
    Smith, T.F., Waterman, M.S. 1981. Identification of common molecular subsequences. J. Mol. Biol., 147:195-197.
    </p>""", unsafe_allow_html=True)

    st.subheader("**Alignment:**")
    tipo_alinhamento = st.radio("Select:", ["Global", "Local"]) 

    st.subheader("**Enter below your DNA sequences:**")       
    input_seq1 = st.text_area(label=">>> String 1:", 
                        value=default_input1, height=100, help=help_text)
    input_seq2 = st.text_area(label=">>> String 2:", 
                        value=default_input2, height=100, help=help_text)

    if input_seq1 and input_seq2:        
        seq1 = formata_sequencia(input_seq1)
        seq2 = formata_sequencia(input_seq2)
        
        if eh_dna(seq1) and eh_dna(seq2):
            matriz_subs = cria_matriz_subs()        

            if tipo_alinhamento == "Global":
                matriz_caminho = lcs_global(seq1, seq2, matriz_subs)
                match, mismatch, gap, score_final, ali_seq1, ali_seq2 = alinhamento_global(seq1, seq2, matriz_caminho, matriz_subs)
            
            elif tipo_alinhamento == "Local":
                matriz_pontuacao, matriz_caminho = lcs_local(seq1, seq2, matriz_subs)
                match, mismatch, gap, score_final, ali_seq1, ali_seq2 = alinhamento_local(seq1, seq2, matriz_pontuacao, matriz_caminho, matriz_subs)

            df_quantidade = df_quantidade_nucleotideos(seq1, seq2)
            df_porcentagem = df_porcentagem_nucleotideos(seq1, seq2)
            df_gc = conteudo_gc(seq1, seq2)

            grafico_quantidade = gera_grafico(df_quantidade, "Quantity") 
            grafico_porcentagem = gera_grafico(df_porcentagem, "Percentage") 
            grafico_gc = gera_grafico_gc(df_gc)

            if tipo_alinhamento == "Global":            
                st.subheader("**# 1. Global alignment:**")
            elif tipo_alinhamento == "Local":
                st.subheader("**# 1. Local alignment:**")

            st.markdown(f"(1) {ali_seq1}<br>(2) {ali_seq2}", unsafe_allow_html=True)
            st.markdown(f"""<p><i>Matches: {match}<br>Mismatches: {mismatch}<br>Gaps: {gap}</i><br>
                        <i><b>Final score: {score_final}""", unsafe_allow_html=True)   
                             
            st.subheader("**# 2. Nucleotides:**")            
            selecao = st.radio("Select:", ["Quantity", "Percentage"])
            resposta = st.radio("Charts?", ["Yes", "No"], index=1) 

            if selecao == "Quantity":
                st.dataframe(df_quantidade)
                if resposta == "Yes":
                    st.pyplot(grafico_quantidade)

            elif selecao =="Percentage":
                st.dataframe(df_porcentagem)
                if resposta == "Yes":
                    st.pyplot(grafico_porcentagem)
                    
            st.subheader("**# 3. GC Content:**")
            st.dataframe(df_gc)
            if resposta == "Yes":
                st.pyplot(grafico_gc)

        else:
            st.write("*Your entered sequences are not DNA sequences =(*")
            st.write("*Please, try again*")  

        
if __name__ == "__main__":
    main()