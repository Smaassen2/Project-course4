from flask import Flask, render_template, request
import mysql.connector
app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def website():
    """De zoekterm die op de website is meegegeven wordt in de database
    gezocht en in een lijst gezet en deze lijst wordt gereturned naar
    website.html.

    :return: website.html - html
    :return: zoekterm - string
    :return: filter - string
    :return: results - list
    :return: len_results - integer
    """
    # Er wordt gekeken of er gebruik wordt gemaakt van de POST methode
    try:
        if request.method == "POST":
            zoekterm = request.form.get("zoekterm", "")
            filter = request.form.get("filter", '')
            # Connecten met de database
            conn = mysql.connector.connect(host="145.74.104.145",
                                           user="oppvr",
                                           password="A1b2c3d4!",
                                           auth_plugin=
                                           'mysql_native_password',
                                           db="oppvr")
            # Openen van een cursor
            cursor = conn.cursor()
            # Uitvoeren van een query
            # De data die deze zoekterm bevat wordt in de
            # De_blast_resultaten tabel in de kolom description_protein
            # gezocht in de rijen die bij de forward reads horen
            # Gaat over de 72-karakterlijn heen omdat de sql-code anders
            # niet werkt
            if filter == "forward":
                cursor.execute("select * from De_blast_resultaten where description_protein like '%" + zoekterm + "%' and De_sequentie_data_seq_id % 2 = 1")
            # De data die deze zoekterm bevat wordt in de
            # De_blast_resultaten tabel in de kolom description_protein
            # gezocht in de rijen die bij de reverse reads horen
            elif filter == "reverse":
                cursor.execute("select * from De_blast_resultaten where description_protein like '%" + zoekterm + "%' and De_sequentie_data_seq_id % 2 = 0")
            # De data die deze zoekterm bevat wordt in de
            # De_blast_resultaten tabel in de kolom description_protein
            # gezocht in zowel de forward als reverse reads
            elif filter == "":
                cursor.execute("select * from De_blast_resultaten where description_protein like '%" + zoekterm + "%'")
            # De gefilterde data wordt opgehaald
            alle_rijen = cursor.fetchall()
            results = []
            # De opgehaalde data wordt in een lijst gezet
            for data in alle_rijen:
                results.append(data)
            # Sluiten van de cursor en connectie
            cursor.close()
            conn.close()
            # De relevante data wordt gereturned naar website.html
            return render_template("website.html", zoekterm=zoekterm,
                                   filter=filter, results=results,
                                   len_results=len(results))
    except IndexError:
        print("Er is een invalide index gevonden in de functie"
              "website")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "website")
    except TypeError:
        print("Een of meer variabelen in de functie website"
              "heeft/hebben niet de correcte variabele type(n).")

    # Alle variabelen worden leeg gereturnd
    else:
        return render_template("website.html", zoekterm="", filter="",
                               results=[], len_results=0)

@app.route('/image.html')
def image():
    # Image.html wordt gereturned zodat deze op de website.html geopend
    # kan worden
    return render_template('image.html')


if __name__ == '__main__':
    app.run()
