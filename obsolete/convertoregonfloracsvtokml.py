#!/usr/bin/env python2.3
#
#  oregon_flora_csv_to_kml.py
#  
#  Created by Cody Hinchliff on 3/31/11.
#

import csv

name = "Oregon Flora Apiaceae"

infile = file("results.csv","rU")
data = csv.DictReader(infile)

outfile = file("data.kml","wb")

indent = 0

outfile.write("<?xml version='1.0' encoding='ISO-8859-1'?>" \
"<kml xmlns='http://www.opengis.net/kml/2.2'>\n" \
"\t<Document>\n" \
"\t\t<Folder>\n" \
"\t\t\t<name>{0}</name>\n".format(name))

indent += 3

for record in data:
	
	outfile.write("{0}<Placemark>\n".format(indent * "\t"))
	indent += 1
	
	outfile.write("{0}<name>{1}</name>\n".format(indent * "\t", record["taxon"]))

	try:
		taxon_auths = record["taxonauths"].split(record["taxon"])[1]
	except IndexError:
		taxon_auths = ""

	description = "<b>" + record["taxon"] + "</b>" + taxon_auths + "\n\n " + record["county"] + \
		". " + record["location"] + " " + record["habitat"] + " " + record["phenology"] + ".\n\n"
	
	if record["data_type"] == "obs":
		description += "Observation. " + record["collector_observer"] + " on " + record["date"] + ", " + \
			record["year"] + "."

	elif record["data_type"] == "voucher, not at OSU":
		description += "Observation? " + record["collector_observer"] + " on " + record["date"] + ", " + \
			record["year"] + ". Voucher not at OSU; voucher number unknown."

	else:
		description += "Vouchered collection. " + record["collector_observer"] + " " + record["collectors_num"] + \
			". " + record["date"] + ", " + record["year"] + "."	

	outfile.write("{0}<description><![CDATA[<p>{1}]]></description>\n".format(indent * "\t", description))
	
#	if record['taxon'].split(" ")[0] == "Lomatium":
	outfile.write("{0}<styleUrl>#Lomatium</styleUrl>\n".format(indent * "\t"))

	outfile.write("{0}<Point><coordinates>{1},{2},{3}</coordinates></Point>\n".format(indent * "\t", \
		record["long"], record["lat"], record["elevation (m)"]))
	
	indent -= 1

	outfile.write("{0}</Placemark>\n".format(indent * "\t"))

outfile.write("\t\t</Folder>\n" \
"\t</Document>\n" \
"</kml>")

outfile.close()
