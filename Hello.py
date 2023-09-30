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
    data = response.json()
    
    if response.status_code != 200:
        st.error(f"Failed to retrieve news: {response.text}")
        return []
    else:
        return data['articles']

def main():
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
            
            if not articles:
                st.write("No articles found.")
            
            for article in articles:
                article_md = f"""
### [{article['title']}]({article['url']})
- **Source:** {article['source']['name']}
- **Published At:** {article['publishedAt']}
- **Description:** {article['description']}
---
"""
                st.markdown(article_md)

if __name__ == "__main__":
    main()
