#!/usr/bin/python3

import requests
from bs4 import BeautifulSoup
import re

def parse_url(url,headers=None):
    '''
    it receive a url as parameter
    it will parse the information
    it return a dict including: title, post_date, author, url, contents
    '''
    rs = requests.get(url.strip(), headers=headers)
    soup = BeautifulSoup(rs.content, "lxml")
    title = soup.title.get_text()
    # extract the post date
    if soup.find(id='post-date') is not None:
        post_date = soup.find(id="post-date").get_text().strip()
    else:
        post_date = soup.find(id='publish_time').get_text().strip()
    # extract the author
    if soup.find(id="meta_content") is not None:
        author = soup.find(id="meta_content").find_all('em')[1].get_text().strip()
    else:
        author = "分享"
    url = url.strip()
    page_content = [p.get_text().strip() for p in soup.find_all("p")]
    info = dict(title = title,
                time = post_date,
                author = author,
                url = url,
                contents = ''.join(page_content)
               )
    return info


def load_db_mysql(dicts, host, user, passwd, dbname):
    '''
    it receive a dicts and database connection as parameters
    it will insert the data to the databases
    the table of the db must include time, author, url, title, contents
    '''
    import pymysql
    connection = pymysql.connect(
    host = host,
    user = user,
    password = passwd,
    db = dbname,
    charset = 'utf8mb4',
    cursorclass=pymysql.cursors.DictCursor
    )
    try:
        with connection.cursor() as cursor:
            # Create a new record
            sql = '''
            INSERT INTO `articles` (`time`, `author`, `url`, `title`, `content`)
            VALUES (%s, %s, %s, %s, %s)
            '''
            cursor.execute(sql, (dicts['time'], dicts['author'],dicts['url'],
                                dicts['title'], dicts['contents']))

        # connection is not autocommit by default. So you must commit to save
        # your changes.
        connection.commit()
    finally:
        connection.close()

def load_db_sqlite(dicts, dbname, table):
    '''
    dicts: a dictionary
    dbname: sqlite3 database file name
    '''
    import sqlite3
    conn = sqlite3.connect(dbname)
    curs = conn.cursor()
    curs.execute('''
        create table IF NOT EXISTS ?
        (time char(100), author char(10), url char(100), title char(150), content char(10000))
        ''', table)
    curs.execute('''
    INSERT INTO `articles` (`time`, `author`, `url`, `title`, `content`)
    VALUES (?,?,?,?,?)
    ''', dicts)
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
    sqltype = sys.argv[2]
    dbname = sys.argv[3]
    print(url)
    info = parse_url(url=url, headers=headers)
    if sqltype.lower() == 'mysql':
        import os
        host =  os.getenv('host')
        if not os.getenv('user'):
            user = 'root'
        else:
            user = os.getenv('user')
        password = os.getenv('password')
        load_db_mysql(info ,host=host,user=user,passwd=password, dbname=dbname)
    else:
        load_db_sqlite(info, dbname, table='articles')
