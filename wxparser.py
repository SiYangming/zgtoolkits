#!/usr/bin/python3

import requests
from bs4 import BeautifulSoup
import re

def getSoup(url, headers):
    response = requests.get(url=url, headers=headers)
    soup = BeautifulSoup(response.text, "lxml")
    return soup

def parsePage(soup):
    title= soup.select_one("#activity-name").get_text(strip=True)
    date = soup.select_one("#post-date").get_text(strip=True)
    author = soup.find_all("em")[1].get_text(strip="True")
    post_name = soup.select_one("#post-user").get_text(strip=True)
    article = [p.get_text(strip=True) for p in soup.find_all("p")]
    print(title)
    return [title, date, author, post_name, '\n'.join(article)]


def login_db(dbfile):
    import sqlite3
    conn = sqlite3.connect(dbfile)
    curs = conn.cursor()
    return conn, curs

def load_db(curs, table, data, conn=None):
    dataset = data
    curs.execute('create table IF NOT EXISTS ' + table + \
        ' (title char(100), date char(10), author char(20), poster char(100), article char(10000))')
    curs.execute('insert into ' + table + ' values(?, ?, ?, ?, ?)', dataset)

def logout_db(conn):
    conn.commit()
    conn.close()


if __name__ == '__main__':
    import sys
    headers= {
    "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8",
    "Host":"mp.weixin.qq.com",
    "User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/59.0.3071.115 Safari/537.36"
    }
    url = sys.argv[1]
    print(url)
    soup = getSoup(url, headers)
    content = parsePage(soup)
    conn, curs = login_db("bioinfo")
    load_db(curs, table = "bioinfo", data = content, conn = conn)
    logout_db(conn)
