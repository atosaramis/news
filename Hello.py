import streamlit as st
import requests
from datetime import datetime

# Requirements.txt file:
# streamlit
# requests


def fetch_news(api_key, query, from_date, to_date, num_results):
  """Fetches news articles from the NewsAPI.

  Args:
    api_key: The NewsAPI key.
    query: The search query.
    from_date: The start date for the search.
    to_date: The end date for the search.
    num_results: The number of results to return.

  Returns:
    A list of news articles.
  """

  url = "https://newsapi.org/v2/everything?q={}&from={}&to={}&apiKey={}&pageSize={}"
  response = requests.get(url.format(query, from_date, to_date, api_key, num_results))
  articles = response.json()['articles']
  return articles


menu = st.sidebar.selectbox("Choose a section", ["News API Searcher"])

if menu == "News API Searcher":
  st.title("News API Searcher")

  # Validate the user input.
  api_key = st.text_input("Enter your NewsAPI key:")
  if not api_key:
    st.error("Please enter a valid NewsAPI key.")
    st.stop()

  # Display the search form.
  query = st.text_input("Search Term:", "COVID-19")
  from_date = st.date_input("From Date:", datetime.now().date())
  to_date = st.date_input("To Date:", datetime.now().date())
  num_results = st.number_input("Number of Results to Return:", 1, 100, 10)

  # Fetch the search results.
  if st.button("Fetch Articles"):
    st.write("Fetching articles...")
    articles = fetch_news(api_key, query, from_date, to_date, num_results)

    # Display the search results.
    if articles:
      st.write("Search results:")
      for article in articles:
        st.write(f"### [{article['title']}]({article['url']})")
        st.write(f"- Source: {article['source']['name']}")
        st.write(f"- Published At: {article['publishedAt']}")
