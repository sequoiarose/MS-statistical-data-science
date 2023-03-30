# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 19:38:47 2022

@author: sequo
"""
import os

from flask import Flask, render_template, url_for, request
import pymysql

# create and configure the app
app = Flask(__name__)
app.config.update(
    SECRET_KEY='192b9bdd22ab9ed4d12e236c78afcb9a393ec15f71bbf5dc987d54727823bcbf'
)

def setup_cursor():
    connection = pymysql.connect(host="localhost", port=3306, user="root", passwd="", database="vet_clinic")
    cursor = connection.cursor()
    return cursor, connection

# a simple page that says hello
@app.route('/')
def home():
    return render_template("home.html")

@app.route('/data_records')
def records():
    return render_template("home.html")

@app.route('/data_records/pets')
def records_pets():
    # database connection
    cursor, connection = setup_cursor()
    query = "select * from pets"
    cursor.execute(query)
    data = cursor.fetchall()
    # some other statements  with the help of cursor
    connection.close()
    return render_template("data.html", data=data)

@app.route('/data_records/owners')
def records_owners():
    # database connection
    cursor, connection = setup_cursor()
    query = "select * from owners"
    cursor.execute(query)
    data = cursor.fetchall()
    # some other statements  with the help of cursor
    connection.close()
    return render_template("data_owners.html", data=data)

@app.route('/data_records/vets')
def records_vets():
    # database connection
    cursor, connection = setup_cursor()
    query = "select * from vets"
    cursor.execute(query)
    data = cursor.fetchall()
    # some other statements  with the help of cursor
    connection.close()
    return render_template("data_vets.html", data=data)


@app.route('/data_records/appointments')
def records_appts():
    # database connection
    cursor, connection = setup_cursor()
    query = "select a.AID, a.appt_date, a.appt_time, a.exam_room, pa.PID, va.EID from appointments a, pet_appointments pa, vet_appointments va \
        where a.AID=pa.AID and a.AID=va.AID"
    cursor.execute(query)
    data = cursor.fetchall()
    # some other statements  with the help of cursor
    connection.close()
    return render_template("data_appointments.html", data=data)

@app.route('/data_records/vax')
def records_vax():
    # database connection
    cursor, connection = setup_cursor()
    query = "select v.VID, pv.vax_date, v.vax_type, pv.PID from vaccines v, pet_vaccines pv \
        where v.VID=pv.VID"
    cursor.execute(query)
    data = cursor.fetchall()
    # some other statements  with the help of cursor
    connection.close()
    return render_template("data_vax.html", data=data)

@app.route('/data_input')
def data_in():
    return render_template("home.html")

@app.route('/data_input/pets', methods=["POST", "GET"])
def data_input_pets():
    cursor, connection = setup_cursor()
    pets_sql = "SELECT OID FROM Owners"
    cursor.execute(pets_sql)
    owners = cursor.fetchall()
    connection.close()
    if (request.method=='POST'):
        cursor, connection = setup_cursor()
        PID_query = "SELECT count(PID) FROM PETS"
        cursor.execute(PID_query)
        PID = cursor.fetchone()
        PID = str(int(PID[0])+1)
        first_name = request.form['first_name']
        last_name = request.form['last_name']
        gender = request.form['gender']
        weight = float(request.form['weight'])
        age = int(request.form['age'])
        species = request.form['species']
        breed = request.form['breed']
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO pets (PID, first_name, last_name, gender, weight, age, species, breed) VALUES
                       ("%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s")""" 
        % (PID, first_name, last_name, gender, weight, age, species, breed))
        connection.commit()
        connection.close()
        
        OID = request.form['owner']
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO Pet_Owners(PID, OID) VALUES
                       ("%s", "%s")""" 
        % (PID, OID))
        connection.commit()
        connection.close()
    return render_template("data_input_pets.html", owners=owners)

@app.route('/data_input/owners', methods=["POST", "GET"])
def data_input_owners():
    if (request.method=='POST'):
        cursor, connection = setup_cursor()
        OID_query = "SELECT count(OID) FROM owners"
        cursor.execute(OID_query)
        OID = cursor.fetchone()
        OID = str(int(OID[0])+1)
        first_name = request.form['first_name']
        last_name = request.form['last_name']
        address = request.form['address']
        phone_number = float(request.form['phone_number'])
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO owners (OID, first_name, last_name, address, phone_number) VALUES
                       ("%s", "%s", "%s", "%s", "%s")""" 
        % (OID, first_name, last_name, address, phone_number))
        connection.commit()
        connection.close()
        
    return render_template("data_input_owners.html")

@app.route('/data_input/vets', methods=["POST", "GET"])
def data_input_vets():
    if (request.method=='POST'):
        cursor, connection = setup_cursor()
        EID_query = "SELECT max(EID) FROM vets"
        cursor.execute(EID_query)
        EID = cursor.fetchone()
        EID = str(int(EID[0])+100)
        first_name = request.form['first_name']
        last_name = request.form['last_name']
        address = request.form['address']
        phone_number = float(request.form['phone_number'])
        liscense = request.form['liscense']
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO vets (EID, first_name, last_name, address, phone_number, license_number) VALUES
                       ("%s", "%s", "%s", "%s", "%s", "%s")""" 
        % (EID, first_name, last_name, address, phone_number, liscense))
        connection.commit()
        connection.close()
        
    return render_template("data_input_vets.html")

@app.route('/data_input/appointments', methods=["POST", "GET"])
def data_input_appts():
    cursor, connection = setup_cursor()
    pets_sql = "SELECT PID FROM pets"
    cursor.execute(pets_sql)
    pets = cursor.fetchall()
    connection.close()
    cursor, connection = setup_cursor()
    vets_sql = "SELECT EID FROM vets"
    cursor.execute(vets_sql)
    vets = cursor.fetchall()
    connection.close()
    if (request.method=='POST'):
        cursor, connection = setup_cursor()
        AID_query = "SELECT count(AID) FROM appointments"
        cursor.execute(AID_query)
        AID = cursor.fetchone()
        AID = str(int(AID[0])+1)
        year = request.form['appt_year']
        month = request.form['appt_month']
        day = request.form['appt_day']
        appt_date = year+"-"+month+"-"+day
        appt_time = request.form['appt_time']
        exam_room = request.form['exam_room']
        PID = request.form['pet']
        EID = request.form['vet']
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO appointments (AID, exam_room, appt_date, appt_time) VALUES
                       ("%s", "%s", "%s", "%s")""" 
        % (AID, exam_room, appt_date, appt_time))
        connection.commit()
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO Vet_Appointments(EID, AID) VALUES
                       ("%s", "%s")""" 
        % (EID, AID))
        connection.commit()
        connection.close()
        cursor, connection = setup_cursor()
        cursor.execute("""INSERT INTO Pet_Appointments(AID, PID) VALUES
                       ("%s", "%s")""" 
        % (AID, PID))
        connection.commit()
        connection.close()
        
        
    return render_template("data_input_appointments.html", pets=pets, vets=vets)

@app.route('/data_input/pet_vax', methods=["POST", "GET"])
def data_input_pet_vax():
    if (request.method=='POST'):
        print(request.form.keys())
        cols = ", ".join(request.form.keys())
        cols = "( "+cols+")"
        print(cols)
        vals = ", ".join(request.form.values())
        vals = "( "+vals+")"
        print(vals)
        cursor, connection = setup_cursor()
        
    return render_template("data_input_pet_vax.html")


if (__name__ == "__main__"):
    app.run(debug=True)