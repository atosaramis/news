import streamlit as st
import requests
from datetime import datetime
import pandas as pd
import base64

# Function to get the news articles
def get_articles(api_key, keyword, from_date, to_date, num_results):
    url = "https://newsapi.org/v2/everything"
    params = {
        "q": keyword,
        "from": from_date,
        "to": to_date,
        "pageSize": num_results,
        "apiKey": api_key
    }
    response = requests.get(url, params=params)
    return response.json()

# Streamlit app
st.title('News Search')

api_key = st.text_input('Enter your NewsAPI API key', type="password")  # Hide the API key

if api_key:
    keyword = st.text_input('Enter a keyword to search for')
    from_date = st.date_input('Choose a start date')
    to_date = st.date_input('Choose an end date')
    num_results = st.number_input('Enter the number of results to return', min_value=1, max_value=100)

    if keyword and from_date and to_date and num_results:
        articles = get_articles(api_key, keyword, from_date.strftime('%Y-%m-%d'), to_date.strftime('%Y-%m-%d'), num_results)
        df = pd.DataFrame(articles['articles'])
        st.markdown(df.to_markdown())
