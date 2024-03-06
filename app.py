import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

#make page layout wide
st.set_page_config(layout="wide")

def a_simulator(number_of_generations, popsize, simulations, WAA, WAa, Waa, frequency_A):
    """Callable function 'a_simulator' which takes as inputs the number of generations, population size, number of simulationos, 
    fittnes and frequency and returns a graph showing the allele frequency at each generation, and a histogram with the final values 
    of the allele frequency."""
    simulation_results = []
    for sim in range(simulations):
        allele_freqs = [frequency_A] 
        for generation in range(1, number_of_generations):
            FA = allele_freqs[-1]
            Fa = 1 - FA

            # Frequency of each genotype in the next gen
            M = WAA*(FA**2) + WAa*(2 * FA * Fa) + Waa*(Fa**2)
            
            # make sure that the probability ends up inbetween 0/1 - dealing with small fitness values
            pA = (((WAA * (FA**2)) / M))
            pA = max(0, min(pA, 1))  # ensure that pA is always between 0 and 1, if less than 0 make 0, more than 1 make 1
            pAa = 0.5 * ((WAa * (2 * FA * Fa)) / M)
            pAa = max(0, min(pAa, 1))
            
            # Use binomial distribution to account for random fluctuations
            fA_freq = (np.random.binomial(popsize, pA) + np.random.binomial(popsize, pAa)) / popsize

            # append to list
            allele_freqs.append(fA_freq)
        
        #append list to list :)
        simulation_results.append(allele_freqs)

    ##
    # FIGURES
    ##
            
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Simulation Results", "Final Frequencies Histogram"))

    # Frequency Scatter
    for i, freqs in enumerate(simulation_results):
        fig.add_trace(
            go.Scatter(x=list(range(number_of_generations)), y=freqs, mode='lines', name=f"Simulation {i+1}"),
            row=1, col=1
        )

    # Histogram
    final_frequencies = [freqs[-1] for freqs in simulation_results]
    fig.add_trace(
        go.Histogram(x=final_frequencies, nbinsx=30, name='A Frequency'),
        row=2, col=1
    )

    fig.update_layout(plot_bgcolor='white', width=700, height=700)
    fig.update_xaxes(range=[0, 1], row=2, col=1)

    return fig

left_col, spacer, right_col = st.columns([1, 0.5, 2])

with left_col:
    st.header("User Inputs")
    generations = st.slider("Enter a number for the number of generations:", 1, 500, value=100, key='gens')
    pop = st.slider("Enter a number for the population size:", 1, 5000, value=500, key='pop')
    sims = st.slider("Enter a number for the # of simulations:", 1, 500, value=80, key='sim')
    AAfit = st.slider("Fitness of AA allele:", min_value=0.0, max_value=1.0, step=0.05, value=1.0, key='WAA')
    Aafit = st.slider("Fitness of Aa allele:", min_value=0.0, max_value=1.0, step=0.05, value=1.0, key='WAa')
    aafit = st.slider("Fitness of aa allele:", min_value=0.0, max_value=1.0, step=0.05, value=1.0, key='Waa')
    a_fre = st.slider("Frequency of A allele:", min_value=0.0, max_value=1.0, step=0.05, value=0.5, key='frequency_A')

with right_col:
    st.header("Pop Gen Sim & Histogram")
    graph = a_simulator(generations, pop, sims, AAfit, Aafit, aafit, a_fre)
    st.plotly_chart(graph)