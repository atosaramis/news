import streamlit as st
from datetime import datetime
import requests
from Bio import Entrez

def fetch_news(api_key, query, from_date, to_date, num_results):
    # ... rest of your fetch_news code ...

def fetch_pubmed(query, from_date, to_date, num_results):
    Entrez.email = st.session_state['email']  # Use email from session state
    # ... rest of your fetch_pubmed code ...

menu = st.sidebar.selectbox("Choose a section", ["PubMed Searcher", "News API Searcher"])

if menu == "PubMed Searcher":
    st.title("PubMed Searcher")
    if 'email' not in st.session_state:
        st.session_state['email'] = st.text_input("Enter your email:")

    if st.session_state['email']:
        query = st.text_input("Search Term:", "COVID-19")
        from_date = st.date_input("From Date:", datetime.now().date())
        to_date = st.date_input("To Date:", datetime.now().date())
        num_results = st.number_input("Number of Results to Return:", 1, 100, 10)

        if st.button("Fetch Articles"):
            # ... rest of your PubMed code ...

elif menu == "News API Searcher":
    st.title("News API Searcher")
    api_key = st.text_input("Enter your NewsAPI key:", type="password")

    if api_key:
        query = st.text_input("Search Term:", "Python")
        date_range = st.date_input("Date Range:", [datetime.now().date(), datetime.now().date()])
        num_results = st.number_input("Number of Results to Return:", 1, 100, 10)

        if st.button("Fetch News"):
            from_date = date_range[0].strftime('%Y-%m-%d')
            to_date = date_range[1].strftime('%Y-%m-%d')
            articles = fetch_news(api_key, query, from_date, to_date, num_results)

            for article in articles:
                st.write(f"### [{article['title']}]({article['url']})")
                st.write(f"- Source: {article['source']['name']}")
                st.write(f"- Published At: {article['publishedAt']}")
                st.write(f"- Description: {article['description']}")
                st.write("---")
