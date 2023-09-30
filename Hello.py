import streamlit as st
import requests
from datetime import datetime

def fetch_news(api_key, query, from_date, to_date, num_results):
    url = "https://newsapi.org/v2/everything"
    params = {
        "q": query,
        "from": from_date,
        "to": to_date,
        "apiKey": api_key,
        "pageSize": num_results
    }
    response = requests.get(url, params=params)
    if response.status_code != 200:
        st.error(f"Failed to retrieve news: {response.text}")
        return []
    else:
        return response.json()['articles']

def main():
    st.title("News API Searcher")
    
    # Initialize session state variables
    if 'api_key' not in st.session_state:
        st.session_state.api_key = ""
    if 'query' not in st.session_state:
        st.session_state.query = "Python"
    if 'date_range' not in st.session_state:
        st.session_state.date_range = [datetime.now().date(), datetime.now().date()]
    if 'num_results' not in st.session_state:
        st.session_state.num_results = 10
    
    st.session_state.api_key = st.text_input("Enter your NewsAPI key:", value=st.session_state.api_key, type="password")
    st.session_state.query = st.text_input("Search Term:", value=st.session_state.query)
    st.session_state.date_range = st.date_input("Date Range:", value=st.session_state.date_range)
    st.session_state.num_results = st.number_input("Number of Results to Return:", 1, 100, value=st.session_state.num_results)
    
    if st.button("Fetch News"):
        from_date = st.session_state.date_range[0].strftime('%Y-%m-%d')
        to_date = st.session_state.date_range[1].strftime('%Y-%m-%d')
        
        articles = fetch_news(st.session_state.api_key, st.session_state.query, from_date, to_date, st.session_state.num_results)
        
        for article in articles:
            st.write(f"### [{article['title']}]({article['url']})")
            st.write(f"- Source: {article['source']['name']}")
            st.write(f"- Published At: {article['publishedAt']}")
            st.write(f"- Description: {article['description']}")
            st.write("---")

if __name__ == "__main__":
    main()
