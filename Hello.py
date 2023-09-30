import streamlit as st
from datetime import datetime
import requests
from Bio import Entrez

# Requirements.txt file:
# streamlit
# requests
# biopython

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

def fetch_pubmed(query, from_date, to_date, num_results, email):
  """Fetches PubMed articles.

  Args:
    query: The search query.
    from_date: The start date for the search.
    to_date: The end date for the search.
    num_results: The number of results to return.
    email: The user's email address.

  Returns:
    A list of PubMed articles.
  """

  Entrez.email = email
  handler = Entrez.esearch(db="pubmed", term=query, datetype="pdat", mindate=from_date, maxdate=to_date, retmax=num_results)
  record = Entrez.read(handler)
  pmids = record["IdList"]

  articles = []
  for pmid in pmids:
    article = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
    article_record = Entrez.read(article)[0]

    articles.append({
      "title": article_record["Article"]["ArticleTitle"],
      "abstract": article_record["Article"]["Abstract"]["AbstractText"][0],
      "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}",
      "authors": [author["LastName"] + ", " + author["FirstName"] for author in article_record["Article"]["AuthorList"]],
      "journal": article_record["Article"]["Journal"],
      "publication_date": article_record["Article"]["PublicationDate"],
    })

  return articles

menu = st.sidebar.selectbox("Choose a section", ["PubMed Searcher", "News API Searcher"])

if menu == "PubMed Searcher":
  st.title("PubMed Searcher")

  # Validate the user input.
  email = st.text_input("Enter your email:", type="email")
  if not email:
    st.error("Please enter a valid email address.")
    st.stop()

  # Display the search form.
  query = st.text_input("Search Term:", "COVID-19
