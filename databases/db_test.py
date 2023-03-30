# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 19:18:07 2022

@author: sequo
"""
import pymysql

# database connection
connection = pymysql.connect(host="localhost", port=3306, user="root", passwd="M3murL3mur!", database="vet_clinic")
cursor = connection.cursor()
query = "select * from pets"
cursor.execute(query)
output = cursor.fetchall()
print(output)
# some other statements  with the help of cursor
connection.close()
