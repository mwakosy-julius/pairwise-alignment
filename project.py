import requests
import streamlit as st
import pandas as pd
import altair as alt
from altair import *
from streamlit_lottie import st_lottie

# st.set_page_config(page_title="Final Year Project", page_icon=":tada:", layout="wide")

# def load_lottieurl(url):
#     r = requests.get(url)
#     if r.status_code != 200:
#         return None
#     return r.json()

# coding_lottie = load_lottieurl("https://assets5.lottiefiles.com/packages/lf20_fcfjwiyb.json")

def format_sequence(sequence): 
    sequence = sequence.upper()
        
    if sequence[0] == ">":
        sequence = sequence.splitlines()
        sequence = sequence[1:]
        sequence = "".join(sequence).strip()
        sequence = "*" + sequence
        
    else:
        sequence = sequence.splitlines()
        sequence = "".join(sequence).strip()
        sequence = "*" + sequence 
        
    return sequence


def is_dna(seq):
    if set(seq).issubset({"A", "C", "G", "T", "*"}):
        return True
    else:
        return False

def matrix_subs():
    matrix = {"col": ["A", "C", "G", "T"],
                "A": [4, -2, -1, -2],
                "C": [-2, 4, -2, -1],
                "G": [-1, -2, 4, -2],
                "T": [-2, -1, -2, 4]}
    return matrix


def calculate_score(base1, base2, matrix_subs):

    j = matrix_subs["col"].index(base1)

    for base in "ACGT":
        if base2 == base:
            score = matrix_subs[base2][j]
            return score

def maximum_score(base1, base2, side, top, diagonal):

    if (base1 == base2) and (diagonal > side) and (diagonal > top):
        return diagonal
    elif (base1 != base2) and (diagonal > side) and (diagonal > top):
        return diagonal
    elif (side > top) and (side > diagonal):
        return side
    else:
        return top


def create_path(base1, base2, side, top, diagonal):

    if (base1 == base2) and (diagonal > side) and (diagonal > top):
        return "\\"
    elif (base1 != base2) and (diagonal > side) and (diagonal > top):
        return "\\"
    elif (side > top) and (side > diagonal):
        return "-"
    else:
        return "|"

def lcs_global(seq1, seq2, matrix_subs):

    score = []
    path = []
    g = -3

    for i in range(0, len(seq1)):
        score.append([0] * len(seq2))
        path.append([""] * len(seq2))

    for i in range(0, len(seq1)):
        score[i][0] = g * i
        path[i][0] = "|"
    for j in range(0, len(seq2)):
        score[0][j] = g * j
        path[0][j] = "-"
 
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):  

            base1 = seq1[i]
            base2 = seq2[j] 
            s = calculate_score(base1, base2, matrix_subs)
        
            side = score[i][j-1]
            top = score[i-1][j]
            diagonal = score[i-1][j-1]
    
            score[i][j] = maximum_score(base1, base2, side + g, top + g, diagonal + s)
            path[i][j] = create_path(base1, base2, side + g, top + g, diagonal + s)

    return path


def lcs_local(seq1, seq2, matrix_subs):

    score = []
    path = []
    g = -3

    for i in range(0, len(seq1)):
        score.append([0] * len(seq2))
        path.append([""] * len(seq2))

    for i in range(0, len(seq1)):
        path[i][0] = "|"
    for j in range(0, len(seq2)):
        path[0][j] = "-"
 
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):  

            base1 = seq1[i]
            base2 = seq2[j] 
            s = calculate_score(base1, base2, matrix_subs)
        
            side = score[i][j-1]
            top = score[i-1][j]
            diagonal = score[i-1][j-1]
    
            score[i][j] = max(0, maximum_score(base1, base2, side + g, top + g, diagonal + s))
            path[i][j] = create_path(base1, base2, side + g, top + g, diagonal + s)

    return score, path

def global_alignment(seq1, seq2, matrix_path, matrix_subs):
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
        s = calculate_score(seq1[i], seq2[j], matrix_subs)

        if matrix_path[i][j] == "\\" and seq1[i] == seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            match += 1
            score_final += s
            i -= 1
            j -= 1

        elif matrix_path[i][j] == "\\" and seq1[i] != seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            mismatch += 1
            score_final += s
            i -= 1
            j -= 1    
    
        elif matrix_path[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            gap += 1
            score_final += g
            j -= 1
    
        elif matrix_path[i][j] == "|":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            gap += 1
            score_final += g
            i -= 1

    return match, mismatch, gap, score_final, ali_seq1, ali_seq2


def find_highest_value(seq1, matrix):
  
    highest_value_matrix = 0

    for i in range(len(seq1)):
        highest_value_in_line = max(matrix[i])
        if highest_value_in_line > highest_value_matrix:
            highest_value_matrix = highest_value_in_line
            line = i
  
    column = matrix[line].index(highest_value_matrix)
    return line, column


def local_alignment(seq1, seq2, score_matrix, matrix_path, matrix_subs):
    ali_seq1 = ""
    ali_seq2 = ""
    g = -3
    match = 0
    mismatch = 0
    gap = 0
    score_final = 0

    i, j = find_highest_value(seq1, score_matrix)
    value = score_matrix[i][j]

    while value > 0:

        s = calculate_score(seq1[i], seq2[j], matrix_subs)    

        if matrix_path[i][j] == "\\" and seq1[i] == seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            match += 1
            score_final += s
            i -= 1
            j -= 1
            value = score_matrix[i][j]

        elif matrix_path[i][j] == "\\" and seq1[i] != seq2[j]:
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            mismatch += 1
            score_final += s
            i -= 1
            j -= 1    
            value = score_matrix[i][j]
    
        elif matrix_path[i][j] == "-":
            ali_seq1 = " - " + ali_seq1
            ali_seq2 = seq2[j] + ali_seq2
            gap += 1
            score_final += g
            j -= 1
            value = score_matrix[i][j]
    
        elif matrix_path[i][j] == "|":
            ali_seq1 = seq1[i] + ali_seq1
            ali_seq2 = " - " + ali_seq2
            gap += 1
            score_final += g
            i -= 1
            value = score_matrix[i][j]

    return match, mismatch, gap, score_final, ali_seq1, ali_seq2

def percentage_nucleotide(sequence):
    A = sequence.count('A')
    C = sequence.count('A')
    G = sequence.count('A')
    T = sequence.count('A')
    total = A+G+C+T

    percentage_A = A/total*100
    percentage_C = C/total*100
    percentage_G = G/total*100
    percentage_T = T/total*100

    percentages = [percentage_A, percentage_C, percentage_G, percentage_T]
    return percentages

def dataframe():
    try:
        d = dict([
        ('A ', [sequence1.count('A'), sequence2.count('A')]),
        ('G ', [sequence1.count('G'), sequence2.count('G')]),
        ('C ', [sequence1.count('C'), sequence2.count('C')]),
        ('T ', [sequence1.count('T'), sequence2.count('T')])
        ])
    except:
        ZeroDivisionError
    return d

def table():
    X = dataframe()
    X_label = list(X)
    X_values = list(X.values())

    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'sequence1', 1: 'sequence2'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index': 'nucleotide'})
    return df

with st.container():
    left_column, right_column = st.columns((2,1))

    with left_column:
        st.sidebar.markdown("""<p style='text-align: justify'>
        <b>🧬 About this web app</b><br>
        Global and local alignment between two DNA sequences (pairwise alignment).<br><br> 
        <b>⚙️ How to use</b><br>
        Choose the type of alignment to use; Global or Local alignment.<br>
        Copy and paste your sequences one and two into the text entry areas and place Enter.<br>
        <br> 
        <b>🔢 Scoring</b><br> 
        The scores for match and mismatch were based on nucleotide substitution model K2P (Kimura 2-parameters)
        which allows for different rates of transition and transversion.<br>
        <i>> Match = +4</i>, for similarity<br>
        <i>> Mismatch = -1</i>, for transition<br> 
        <i>> Mismatch = -2</i>, for transversion<br> 
        <i>> Gap = -3</i>, for insertion or deletion<br><br>
        <b>📊 Nucleotides content</b><br>
        This web app also counts the amount of Adenine, Cytosine, Guanine and  
        Thymine and plots a bar chart for the sequences' count.</b><br><br> 
        </p>""", unsafe_allow_html=True)

        st.header("Pairwise Alignment Tool")
        
        type_of_alignment = st.radio("Select:", ["Global_Alignment", "Local_Alignment"])

        sequence1 = st.text_area(label="sequence1", height=20, key=1)
        sequence2 = st.text_area(label="sequence2", height=20, key=2)

        if sequence1 and sequence2:        
            seq1 = format_sequence(sequence1)
            seq2 = format_sequence(sequence2)

            if is_dna(seq1):
                if is_dna(seq2):
                    matrix_subs = matrix_subs()        

                    if type_of_alignment == "Global_Alignment":
                        matrix_score = lcs_global(seq1, seq2, matrix_subs)
                        match, mismatch, gap, score_final, ali_seq1, ali_seq2 = global_alignment(seq1, seq2, matrix_score, matrix_subs)
                    
                    elif type_of_alignment == "Local_Alignment":
                        matrix_path, matrix_score = lcs_local(seq1, seq2, matrix_subs)
                        match, mismatch, gap, score_final, ali_seq1, ali_seq2 = local_alignment(seq1, seq2, matrix_path, matrix_score, matrix_subs)
                    
                    if type_of_alignment == "Global_Alignment":            
                        st.subheader("**# 1. Global alignment:**")
                    elif type_of_alignment == "Local_Alignment":
                        st.subheader("**# 1. Local alignment:**")
                    
                    similarity = round(match/(match+mismatch+gap)*100)

                    st.markdown(f"(1) {ali_seq1}<br>(2) {ali_seq2}", unsafe_allow_html=True)
                    st.markdown(f"""<p><i>
                    Matches: {match}<br>
                    Mismatches: {mismatch}<br>
                    Gaps: {gap}</i><br>
                    <i>Final score: {score_final}</i><br>
                    <i><b>Similarity: {similarity}%</i><br>
                    """, unsafe_allow_html=True) 
                else:
                    st.error("Incorrect input in sequence2!")
            else: 
                st.error("Incorrect input in sequence1!")

    with right_column:
        df = table()
#         st_lottie(coding_lottie, height=300, key="coding")
        st.dataframe(df, use_container_width = True)
        try:
            ds = pd.DataFrame([
                ['A', sequence1.count('A'), 'seq1'],
                ['C', sequence1.count('C'), 'seq1'],
                ['G', sequence1.count('G'), 'seq1'],
                ['T', sequence1.count('T'), 'seq1'],
                ['A', sequence2.count('A'), 'seq2'],
                ['C', sequence2.count('C'), 'seq2'],
                ['G', sequence2.count('G'), 'seq2'],
                ['T', sequence2.count('T'), 'seq2']],
                columns = ['Nucleotide', 'Percentage_Count', 'Sequences']
            )

        except:
            ZeroDivisionError

        bar_chart = alt.Chart(ds).mark_bar().encode(
            column = Column('Nucleotide'),
            x=X('Sequences'),
            y=Y('Percentage_Count'),
            color = Color('Sequences', scale=Scale(range=['#EA98D2', '#659CCA']))
            ).configure_view(
            strokeWidth=1.0,
            height=200,
            width=80
                )
 
        st.altair_chart(bar_chart, use_container_width=False)
        #st.bar_chart(df, x=df.columns[0], use_container_width=True)
