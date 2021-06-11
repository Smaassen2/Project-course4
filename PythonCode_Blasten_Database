from Bio.Blast import NCBIXML, NCBIWWW
import mysql.connector
import re


def lijstfilteren(bestandnaam):
    """Opent een TSV bestand. Vanuit hier wordt iedere regel gestript en
    gesplitst op een tab. De sequentie en headers worden aan aparte
    nieuwe lijsten toegevoegd
    :param bestandnaam: TSV File
    :return: l_seq, l_headers
    """
    try:
        # Opent een bestand
        document = open(bestandnaam, "r")
        l_seq = []
        l_headers = []
        # Kijkt naar iedere regel van het bestand
        for line in document:
            if line != "    ":
                l_seq.append(line.strip().split("	")[1])
                l_headers.append(line.strip().split("	")[0])
                l_seq.append(line.strip().split("	")[4])
                l_headers.append(line.strip().split("	")[3])
        # Sluit een bestand
        document.close()
        return l_seq, l_headers
    except FileNotFoundError:
        print("Het bestand is niet gevonden in de functie "
              "lijstfilteren")
    except IOError:
        print("Het bestand is niet leesbaar in de functie "
              "lijstfilteren")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "lijstfilteren")
    except TypeError:
        print("Een of meer variabelen in de functie lijstfilteren "
              "heeft / hebben niet de correcte variabele type(n).")


def blasten(l_seq, xmlbestand):
    """Iedere sequentie in de lijst wordt geblast door middel van
    qblast. De gebruikte blast is Blastx. De gebruikte database is
    non-redundant database. De gebruikte filters zijn een e_value filter
    en hoeveelheid hits. Alle hits zijn hierna toegevoegd aan een XML
    bestand.
    :param xmlbestand: file
    :param l_seq: list
    :return: Alle data van de sequenties in een XML bestand
    """
    try:
        index = 0
        # Loopt door iedere sequentie in de sequentie lijst heen
        for sequence in l_seq:
            print("Start", index)
            print("Start BLAST...")
            # Blasten van iedere sequentie tegen de non-redunant
            # database
            result_handle = NCBIWWW.qblast("blastx", "nr", sequence,
                                           expect=0.01, hitlist_size=10)
            print("BLAST resultaat in variabele")
            # Append de data in een XML bestand
            with open(xmlbestand, "a") as out_handle:
                out_handle.write(result_handle.read())
            print("eind", index)
            index += 1
    except MemoryError:
        print("De operator is out of memory in de functie blasten")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "blasten")
    except TypeError:
        print("Een of meer variabelen in de functie blasten "
              "heeft/hebben niet de correcte variabele type(n).")
    except SyntaxError:
        print("Er is een invalide syntax gevonden in de functie "
              "blasten")


def header_seq(l_seq, l_headers):
    """Het inserten van de sequenties en de headers in de tabel
    De_sequentie_data in de database
    :param l_seq: list
    :param l_headers: list
    :return: De_sequentie_data kolom gevuld
    """
    try:
        count = 1
        # Verbinding aangaan met de database
        conn = mysql.connector.connect(host="145.74.104.145",
                                       user="oppvr",
                                       password="A1b2c3d4!",
                                       auth_plugin=
                                       'mysql_native_password',
                                       db="oppvr")
        cursor = conn.cursor()
        # Tegelijk door de header lijst en sequentie lijst heen loopen
        # door middel van een for loop
        for header, sequentie in zip(l_headers, l_seq):
            # Inserten van de data in de database
            cursor.execute(
                # Gaat over de 72-karakterlijn heen omdat de sql-code
                # anders niet werkt
                "insert into De_sequentie_data (seq_id, name_header, `read`) values "
                "('" + str(
                    count) + "','" + header + "','" + sequentie + "')")
            # Uitvoeren van de execute
            conn.commit()
            count += 1
    except MemoryError:
        print("De operator is out of memory in de functie header_seq")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "header_seq")
    except TypeError:
        print("Een of meer variabelen in de functie header_seq heeft/"
              "hebben niet de correcte variabele type(n).")
    except SyntaxError:
        print("Er is een invalide syntax gevonden in de functie "
              "header_seq")


def parsen(bestandnaam):
    """Opent een XML bestand. Deze wordt geparsd met behulp van de
    NCBIXML parser. Hier worden de values: discriptie, e_value,
    accession, bitscore uitgefilterd. Hierna wordt de scientific_name
    gefilterd uit de discriptie. Tot slot worden al deze values
    geïnsert in de kolom De_blast_resultaten in de database.
    :param bestandnaam: XML File
    :return: De_blast_resultaten kolom gevuld
    """
    try:
        index = 0
        counter = 0
        # Verbinding aangaan met de database
        conn = mysql.connector.connect(host="145.74.104.145",
                                       user="oppvr",
                                       password="A1b2c3d4!",
                                       auth_plugin=
                                       'mysql_native_password',
                                       db="oppvr")
        cursor = conn.cursor()
        # Opent een bestand
        with open(bestandnaam) as out_handle:
            blast_records = NCBIXML.parse(out_handle)
            try:
                # Loopt door iedere read heen
                for blast_record in blast_records:
                    try:
                        index += 1
                        # Loopt door iedere hit heen
                        for alignment in blast_record.alignments:
                            try:
                                # Loopt door iedere hsp heen
                                for hsp in alignment.hsps:
                                    description = alignment.title
                                    e_value = hsp.expect
                                    accession = alignment.accession
                                    bitscore = hsp.bits
                                    counter += 1
                                    scientific_name = re.search(
                                        "\[[^\]]*\]", description)
                                    scientific_name = \
                                        scientific_name.group()
                                    scientific_name = \
                                        scientific_name.strip("[]")
                                    # Inserten van de data in de
                                    # database
                                    # Gaat over de 72-karakterlijn heen
                                    # omdat de sql-code anders niet
                                    # werkt
                                    cursor.execute(
                                        "insert into De_blast_resultaten (blast_id, accession_code, description_protein, scientific_name, bit_score, e_value, De_sequentie_data_seq_id) values "
                                        "('" + str(counter) + "','" +
                                        accession + "','" + description
                                        + "','" + scientific_name
                                        + "','" + str(bitscore) + "','"
                                        + str(e_value) + "','"
                                        + str(index) + "')")
                                    # Uitvoeren van de execute
                                    conn.commit()
                            except NameError:
                                print("Functienaam of variabelenaam "
                                      "bestaat niet in de functie"
                                      "parsen in de hsp loop")
                            except TypeError:
                                print("Een of meer variabelen in de "
                                      "functie parsen in de hsp loop "
                                      "heeft/hebben niet de correcte "
                                      "variabele type(n).")
                    except ValueError:
                        print("Er is geen hit gevonden na ", accession,
                              "in de functie parsen in de alignment"
                              "loop")
            except NameError:
                print("Functienaam of variabelenaam bestaat niet in de "
                      "functie parsen in de blast_record loop")
    except FileNotFoundError:
        print("Het bestand is niet gevonden in de functie parsen")
    except IOError:
        print("Het bestand is niet leesbaar in de functie parsen")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "parsen")
    except IndexError:
        print("Er is een invalide index gevonden in de functie parsen")


def main():
    l_seq, l_headers = lijstfilteren("B4_seq_app.txt")
    blasten(l_seq, "De_definitieve_blast.xml")
    header_seq(l_seq, l_headers)
    parsen("De_definitieve_blast.xml")


main()
