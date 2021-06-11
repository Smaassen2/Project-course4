from Bio.Blast import NCBIXML
import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def scientific_name_filter(bestandnaam):
    """Opent een XML bestand. Deze wordt geparset met behulp van de
    NCBIXML parser. Hier worden de value discriptie uitgefilterd.
    Hieruit wordt de scientific_name gefilterd. De scientific_name wordt
    in een lijst gezet.
    :param bestandnaam: XML File
    :return: l_scientific_name
    """
    try:
        l_scientific_name = []
        # Opent een bestand
        with open(bestandnaam) as out_handle:
            blast_records = NCBIXML.parse(out_handle)
            try:
                # Loopt door iedere read heen
                for blast_record in blast_records:
                    try:
                        # Loopt door iedere hit heen
                        for alignment in blast_record.alignments:
                            try:
                                for hsp in alignment.hsps:
                                    description = alignment.title
                                    scientific_name = re.search(
                                        "\[[^\]]*\]", description)
                                    scientific_name = scientific_name.\
                                        group()
                                    scientific_name = scientific_name.\
                                        strip("[]")
                                    l_scientific_name.append(
                                        scientific_name)
                            except IndexError:
                                print("Er is een invalide index "
                                      "gevonden in de functie "
                                      "scientific_name_filter in de"
                                      "hsp loop")
                            except NameError:
                                print("Functienaam of variabelenaam "
                                      "bestaat niet in de functie "
                                      "scientific_name_filter in de"
                                      "hsp loop")
                    except IndexError:
                        print("Er is een invalide index gevonden in de "
                              "functie scientific_name_filter in de"
                              "alignment loop")
                    except NameError:
                        print("Functienaam of variabelenaam bestaat "
                              "niet in de functie "
                              "scientific_name_filter in de "
                              "alignment loop")
            except IndexError:
                print("Er is een invalide index gevonden in de functie "
                      "scientific_name_filter in de blast_record loop")
            except NameError:
                print("Functienaam of variabelenaam bestaat niet in de "
                      "functie scientific_name_filter in de"
                      "blast_record loop")
        return l_scientific_name
    except FileNotFoundError:
        print("Het bestand is niet gevonden in de functie "
              "scientific_name_filter")
    except IOError:
        print("Het bestand is niet leesbaar in de functie "
              "scientific_name_filter")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "scientific_name_filter")


def scientific_name_frequency(l_scientific_name):
    """
    Maakt van een lijst van organismen een lijst en een dictionary die
    aangeven hoevaak die organismen voorkomen. 1 lijst met absoluten
    waarden en 1 lijst met relatieve waarden.
    :param l_scientific_name: Lijst met organismen namen
    :return: 1 lijst en 1 dictionary die per organisme een waarde bevat
    hoevaak die
    voorkomt. 1 dictionary met absolute waarden (freq) en 1 lijst
    relatieve waarden (frequency)
    """
    try:
        freq = {}
        counter = 1
        frequency = []
        # Maakt een dictionary met als opbouw;
        # (organisme, hoevaak komt het organisme voor in absolute
        # waarde)
        for sn in l_scientific_name:
            freq[sn] = freq.get(sn, 0) + 1
            counter += 1
        # Maakt een lijst met als opbouw;
        # (organisme, hoevaak komt het organisme voor in relatieve
        # waarde)
        for key, value in freq.items():
            f = (key, value / counter * 100)
            frequency.append(f)
        return frequency, freq
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "scientific_name_frequency")
    except KeyError:
        print("Er is een invalide key gevonden in de functie "
              "scientific_name_frequency")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "scientific_name_frequency")


def sorting_one(frequency):
    """ Zet de lijst frequency om in 2 lijsten waarin 1 de top 10
    voorkomende organismen zitten en 1 waarin de bijbehorende relatieve
    frequenties zitten. Bij beide lijsten wordt ook de naam other en de
    bijbehorende waarde toegevoegd
    :param frequency: Lijst die per organisme een relatieve waarde bevat
     over hoevaak die voorkomt
    :return: 2 lijsten waarin 1 de top 10 voorkomende organismen zitten
    (organism) en 1 waarin de bijbehorende relatieve frequentie zit
    (frequency)
    """
    try:
        # Zet per tuple die in elk lijst item zit in de frequency lijst,
        # de eerste waarde in org en de tweede waarde in frq
        org = []
        frq = []
        for item in frequency:
            org.append(item[0])
            frq.append(item[1])
        # Sorteert de frequenties van hoog naar laag waarbij de
        # bijbehorende organismen meegenomen wordt. Daarbij wordt
        # alleen de top 10 opgeslagen
        mx = sorted(zip(frq, org), reverse=True)[:10]
        # Zet per tuple die in elk lijst item zit in de mx lijst,
        # de eerste waarde in frequency en de tweede waarde in organism
        organism = []
        frequency = []
        for item in mx:
            organism.append(item[1])
            frequency.append(item[0])
        # Berekent het percentage van de niet top 10 voorkomende
        # organismen
        other = 100
        for f in frequency:
            other -= f
        # Voegt de naam "Other" toe aan lijst organism en het
        # bijbehorende percentage bij frequency
        organism.append("Other")
        frequency.append(other)
        return frequency, organism
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "sorting_one")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "sorting_one")


def sorting_two(freq):
    """ Zet de dictionary freq om in 2 lijsten waarin 1 de top 10
    voorkomende organismen zitten en 1 waarin de bijbehorende absolute
    frequenties zitten
    :param freq: Dictionary die per organisme een absolute waarde bevat
    over hoevaak die voorkomt
    :return: 2 lijsten waarin 1 de top 10 voorkomende organismen zitten
    (orga) en 1 waarin de bijbehorende absolute frequentie zit (freq)
    """
    try:
        # Zet per dictionary item die in de freq dictionary zit,
        # de key organism in org om en de value in frq
        org = []
        frq = []
        for organism, frequency in freq.items():
            org.append(organism)
            frq.append(frequency)
        # Sorteert de frequenties van hoog naar laag waarbij de
        # bijbehorende organismen meegenomen wordt. Daarbij wordt
        # alleen de top 10 opgeslagen
        mx = sorted(zip(frq, org), reverse=True)[:10]
        # Zet per tuple die in elk lijst item zit in de mx lijst,
        # de eerste waarde in frequency en de tweede waarde in organism
        orga = []
        freq = []
        for item in mx:
            orga.append(item[1])
            freq.append(item[0])
        return freq, orga
    except KeyError:
        print("Er is een invalide key gevonden in de functie "
              "sorting_two")
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "sorting_two")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "sorting_two")


def top_ten_one(frequency, organism):
    """ Maakt een taartdiagram over de top 10 meest voorkomende
    organismen tegenover de rest
    :param frequency: Lijst met de relatieve waardes over de top 10
    voorkomende organismen + de waarde van other
    :param organism: Lijst met de namen van de top 10 organismen
    + "other"
    :return: Een taartdiagram
    """
    try:
        plt.pie(frequency, labels=organism, startangle=90,
                shadow=True, textprops={"fontsize": 6.9},
                colors=["Red", "Orange", "Yellow", "Lime", "Green",
                        "Teal", "Dodgerblue", "Blue", "Blueviolet",
                        "Magenta", "Grey"],
                autopct='%.2f%%', pctdistance=0.85)
        plt.title("Top 10 organism appearances in percentages")
        plt.show()
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "top_ten_one")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "top_ten_one")


def addlabels(orga, freq):
    """Zorgt voor de aanmaak van waarde labels bij de barplot
    :param orga: Lijst met de namen van de top 10 organismen
    :param freq: Lijst met de absolute waardes over de top 10
    voorkomende organismen
    :return: n.v.t.
    """
    try:
        for i in range(len(orga)):
            plt.text(i, freq[i]//2, freq[i], ha="center")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "addlabels")
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "addlabels")


def top_ten_two(freq, orga):
    """ Maakt een barplot over de top 10 meest voorkomende organismen
    :param orga: Lijst met de namen van de top 10 organismen
    :param freq: Lijst met de absolute waardes over de top 10
    :return: Een barplot
    """
    try:
        # Labels x-as
        label_org = ["Gam. b.", "Pro. b.", "Fir. b.", "Xan. b.",
                     "Del. b.", "Gem. b.", "Chi. b.", "Chl. b.",
                     "Aci. b.", "Sin. a."]
        # Kleuren van de bars
        colors = ["Red", "Orange", "Yellow", "Lime", "Green", "Teal",
                  "Dodgerblue", "Blue", "Blueviolet", "Magenta",
                  "Black"]
        plt.bar(orga, freq, color=colors)
        plt.xlabel("Species")
        plt.ylabel("Number of appearences in absolute values")
        # Aanpassingen x-labels
        plt.xticks(orga, label_org, rotation='vertical')
        # Maakt witte horizontale rasterlijnen aan
        plt.grid(True, color="w", axis="y")
        # Maakt een custom legenda aan
        legend = [Line2D([0], [0], color=colors[0], label=orga[0]),
                  Line2D([0], [0], color=colors[1], label=orga[1]),
                  Line2D([0], [0], color=colors[2], label=orga[2]),
                  Line2D([0], [0], color=colors[3], label=orga[3]),
                  Line2D([0], [0], color=colors[4], label=orga[4]),
                  Line2D([0], [0], color=colors[5], label=orga[5]),
                  Line2D([0], [0], color=colors[6], label=orga[6]),
                  Line2D([0], [0], color=colors[7], label=orga[7]),
                  Line2D([0], [0], color=colors[8], label=orga[8]),
                  Line2D([0], [0], color=colors[9], label=orga[9])]
        plt.legend(handles=legend, loc="upper right")
        # Voegt waarde labels toe
        addlabels(orga, freq)
        plt.title("Top 10 organism appearances in absolute values")
        plt.show()
    except IndexError:
        print("Er is een invalide index gevonden in de functie "
              "top_ten_two")
    except NameError:
        print("Functienaam of variabelenaam bestaat niet in de functie"
              "top_ten_two")


def main():
    l_scientific_name = scientific_name_filter("De_definitieve_blast."
                                               "xml")
    frequency, freq = scientific_name_frequency(l_scientific_name)
    frequency, organism = sorting_one(frequency)
    freq, orga = sorting_two(freq)
    top_ten_one(frequency, organism)
    top_ten_two(freq, orga)


main()
